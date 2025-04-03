import re
import arviz as az
import pandas as pd
from glob import glob
from pathlib import Path
from multiprocessing.pool import ThreadPool
#from src._utils import _create_folder_without_clear
from ._utils import _create_folder_without_clear
from .logger import setup_loggers


def _process_dataframe(df, feat_id, suffix=""):
    df = df.copy()
    df.reset_index(inplace=True, drop=True)
    df.columns.name = ""
    df.index = [feat_id]
    df.columns = [x + suffix for x in df.columns]
    return df


def _reformat_multiindex(df, feat_id, suffix="", logger=None):
    df = df.copy().reset_index()
    new_df = pd.DataFrame(columns=df.covariate.unique(), index=[feat_id])
    for c in new_df.columns:
        logger.info(f"DEBUG: Processing covariate {c} for feature {feat_id}")
        logger.info(f"DEBUG: DataFrame shape: {df.shape}, columns: {df.columns.tolist()}")
        logger.info(f"DEBUG: Filtering for covariate {c} and hdi='lower'")
        lower_df = df.loc[(df["covariate"] == c) & (df["hdi"] == "lower")]
        logger.info(f"DEBUG: Lower dataframe shape: {lower_df.shape}")
        if lower_df.empty:
            logger.info(f"DEBUG: No lower values found for covariate {c}")
            continue
        logger.info(f"DEBUG: Accessing 'beta_var' column from lower dataframe")
        lower = lower_df["beta_var"].values[0]
        
        logger.info(f"DEBUG: Filtering for covariate {c} and hdi='higher'")
        higher_df = df.loc[(df["covariate"] == c) & (df["hdi"] == "higher")]
        logger.info(f"DEBUG: Higher dataframe shape: {higher_df.shape}")
        if higher_df.empty:
            logger.info(f"DEBUG: No higher values found for covariate {c}")
            continue
        logger.info(f"DEBUG: Accessing 'beta_var' column from higher dataframe")
        higher = higher_df["beta_var"].values[0]
        
        new_df[c][feat_id] = (lower, higher)
    new_df.columns = [c + suffix for c in new_df.columns]
    return new_df


def _parallel(threads, unit_func, arg_list):
    p = ThreadPool(processes=threads)
    results = p.map(unit_func, arg_list)
    p.close()
    p.join()
    return results

def convert_types(df):
    for col in df.columns:
        try:
            df[col] = df[col].astype(float)
        except ValueError:
            df[col] = df[col].astype(str)
    return df

def summarize_inferences_single_file(inf_file, logger=None):
    FEAT_REGEX = re.compile("F\d{4}_(.*).nc")
    try:
        logger.info(f"DEBUG: Processing file: {inf_file}")
        match = FEAT_REGEX.search(inf_file)
        if match is None:
            logger.info(f"DEBUG: Regex pattern did not match filename: {inf_file}")
            return None
        logger.info(f"DEBUG: Regex match groups: {match.groups()}")
        this_feat_id = match.groups()[0]
        logger.info(f"DEBUG: Extracted feature ID: {this_feat_id}")
        
        logger.info(f"DEBUG: Loading inference data from {inf_file}")
        this_feat_diff = az.from_netcdf(inf_file).posterior["beta_var"]
        logger.info(f"DEBUG: Inference data shape: {this_feat_diff.shape}")
        
        logger.info(f"DEBUG: Calculating mean")
        this_feat_diff_mean = this_feat_diff.mean(["chain", "draw"]).to_dataframe().T
        logger.info(f"DEBUG: Mean dataframe shape: {this_feat_diff_mean.shape}")
        
        logger.info(f"DEBUG: Calculating std")
        this_feat_diff_std = this_feat_diff.std(["chain", "draw"]).to_dataframe().T
        logger.info(f"DEBUG: Std dataframe shape: {this_feat_diff_std.shape}")
        
        logger.info(f"DEBUG: Calculating HDI")
        this_feat_diff_hdi = az.hdi(this_feat_diff).to_dataframe()
        logger.info(f"DEBUG: HDI dataframe shape: {this_feat_diff_hdi.shape}")
        
        logger.info(f"DEBUG: Processing mean dataframe")
        this_feat_diff_mean = _process_dataframe(
            this_feat_diff_mean, this_feat_id, suffix="_mean"
        )
        
        logger.info(f"DEBUG: Processing std dataframe")
        this_feat_diff_std = _process_dataframe(
            this_feat_diff_std, this_feat_id, suffix="_std"
        )
        
        logger.info(f"DEBUG: Reformatting multiindex")
        this_feat_diff_hdis = _reformat_multiindex(
            this_feat_diff_hdi, this_feat_id, suffix="_hdi", logger=logger
        )
        
        logger.info(f"DEBUG: Concatenating dataframes")
        result = pd.concat(
            [this_feat_diff_mean, this_feat_diff_std, this_feat_diff_hdis], axis=1
        )
        logger.info(f"DEBUG: Final result shape: {result.shape}")
        return result
    except Exception as e:
        logger.error(f"Error processing file {inf_file}: {str(e)}")
        import traceback
        logger.error(f"DEBUG: Traceback: {traceback.format_exc()}")
        return None


def summarize_inferences(input_dir, threads=1, logfile=None):
    #_create_folder_without_clear(output_dir)
    logger = setup_loggers(logfile) if logfile else None
    
    all_inf_files = glob(f"{input_dir}/inferences/*.nc")
    logger.info(f"DEBUG: Found {len(all_inf_files)} inference files")

    results = _parallel(threads, lambda x: summarize_inferences_single_file(x, logger), all_inf_files)
    feat_diff_df_list = [df for df in results if df is not None]
    logger.info(f"DEBUG: Processed {len(feat_diff_df_list)} inference files successfully")

    if feat_diff_df_list:
        all_feat_diffs_df = pd.concat(feat_diff_df_list, axis=0)
        all_feat_diffs_df.index.name = "feature id"
        all_feat_diffs_df.to_csv(
            f"{input_dir}/results/beta_var.tsv", sep="\t", index=True
        )
        return convert_types(all_feat_diffs_df)
    else:
        logger.warning("No available feat_diff_dfs...")  # TODO: chaneg this to log
        return None
