import re
import arviz as az
import pandas as pd
from glob import glob
from pathlib import Path
from multiprocessing.pool import ThreadPool
#from src._utils import _create_folder_without_clear
from ._utils import _create_folder_without_clear


def _process_dataframe(df, feat_id, suffix=""):
    df = df.copy()
    df.reset_index(inplace=True, drop=True)
    df.columns.name = ""
    df.index = [feat_id]
    df.columns = [x + suffix for x in df.columns]
    return df


def _reformat_multiindex(df, feat_id, suffix=""):
    df = df.copy().reset_index()
    new_df = pd.DataFrame(columns=df.covariate.unique(), index=[feat_id])
    for c in new_df.columns:
        print(f"DEBUG: Processing covariate {c} for feature {feat_id}")
        print(f"DEBUG: DataFrame shape: {df.shape}, columns: {df.columns.tolist()}")
        print(f"DEBUG: Filtering for covariate {c} and hdi='lower'")
        lower_df = df.loc[(df["covariate"] == c) & (df["hdi"] == "lower")]
        print(f"DEBUG: Lower dataframe shape: {lower_df.shape}")
        if lower_df.empty:
            print(f"DEBUG: No lower values found for covariate {c}")
            continue
        print(f"DEBUG: Accessing 'beta_var' column from lower dataframe")
        lower = lower_df["beta_var"].values[0]
        
        print(f"DEBUG: Filtering for covariate {c} and hdi='higher'")
        higher_df = df.loc[(df["covariate"] == c) & (df["hdi"] == "higher")]
        print(f"DEBUG: Higher dataframe shape: {higher_df.shape}")
        if higher_df.empty:
            print(f"DEBUG: No higher values found for covariate {c}")
            continue
        print(f"DEBUG: Accessing 'beta_var' column from higher dataframe")
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

def summarize_inferences_single_file(inf_file):
    FEAT_REGEX = re.compile("F\d{4}_(.*).nc")
    try:
        print(f"DEBUG: Processing file: {inf_file}")
        match = FEAT_REGEX.search(inf_file)
        if match is None:
            print(f"DEBUG: Regex pattern did not match filename: {inf_file}")
            return None
        print(f"DEBUG: Regex match groups: {match.groups()}")
        this_feat_id = match.groups()[0]
        print(f"DEBUG: Extracted feature ID: {this_feat_id}")
        
        print(f"DEBUG: Loading inference data from {inf_file}")
        this_feat_diff = az.from_netcdf(inf_file).posterior["beta_var"]
        print(f"DEBUG: Inference data shape: {this_feat_diff.shape}")
        
        print(f"DEBUG: Calculating mean")
        this_feat_diff_mean = this_feat_diff.mean(["chain", "draw"]).to_dataframe().T
        print(f"DEBUG: Mean dataframe shape: {this_feat_diff_mean.shape}")
        
        print(f"DEBUG: Calculating std")
        this_feat_diff_std = this_feat_diff.std(["chain", "draw"]).to_dataframe().T
        print(f"DEBUG: Std dataframe shape: {this_feat_diff_std.shape}")
        
        print(f"DEBUG: Calculating HDI")
        this_feat_diff_hdi = az.hdi(this_feat_diff).to_dataframe()
        print(f"DEBUG: HDI dataframe shape: {this_feat_diff_hdi.shape}")
        
        print(f"DEBUG: Processing mean dataframe")
        this_feat_diff_mean = _process_dataframe(
            this_feat_diff_mean, this_feat_id, suffix="_mean"
        )
        
        print(f"DEBUG: Processing std dataframe")
        this_feat_diff_std = _process_dataframe(
            this_feat_diff_std, this_feat_id, suffix="_std"
        )
        
        print(f"DEBUG: Reformatting multiindex")
        this_feat_diff_hdis = _reformat_multiindex(
            this_feat_diff_hdi, this_feat_id, suffix="_hdi"
        )
        
        print(f"DEBUG: Concatenating dataframes")
        result = pd.concat(
            [this_feat_diff_mean, this_feat_diff_std, this_feat_diff_hdis], axis=1
        )
        print(f"DEBUG: Final result shape: {result.shape}")
        return result
    except Exception as e:
        print(f"Error processing file {inf_file}: {str(e)}")
        import traceback
        print(f"DEBUG: Traceback: {traceback.format_exc()}")
        return None


def summarize_inferences(input_dir, threads=1):
    #_create_folder_without_clear(output_dir)
    all_inf_files = glob(f"{input_dir}/inferences/*.nc")

    results = _parallel(threads, summarize_inferences_single_file, all_inf_files)
    feat_diff_df_list = [df for df in results if df is not None]

    if feat_diff_df_list:
        all_feat_diffs_df = pd.concat(feat_diff_df_list, axis=0)
        all_feat_diffs_df.index.name = "feature id"
        all_feat_diffs_df.to_csv(
            f"{input_dir}/results/beta_var.tsv", sep="\t", index=True
        )
        return convert_types(all_feat_diffs_df)
    else:
        print("No available feat_diff_dfs...")  # TODO: chaneg this to log
