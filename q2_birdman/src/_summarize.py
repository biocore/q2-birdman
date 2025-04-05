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
        lower = df.loc[(df["covariate"] == c) & (df["hdi"] == "lower")][
            "beta_var"
        ].values[0]
        higher = df.loc[(df["covariate"] == c) & (df["hdi"] == "higher")][
            "beta_var"
        ].values[0]
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
        this_feat_id = FEAT_REGEX.search(inf_file).groups()[0]
        this_feat_diff = az.from_netcdf(inf_file).posterior["beta_var"]
        this_feat_diff_mean = this_feat_diff.mean(["chain", "draw"]).to_dataframe().T
        this_feat_diff_std = this_feat_diff.std(["chain", "draw"]).to_dataframe().T
        this_feat_diff_hdi = az.hdi(this_feat_diff).to_dataframe()

        this_feat_diff_mean = _process_dataframe(
            this_feat_diff_mean, this_feat_id, suffix="_mean"
        )
        this_feat_diff_std = _process_dataframe(
            this_feat_diff_std, this_feat_id, suffix="_std"
        )
        this_feat_diff_hdis = _reformat_multiindex(
            this_feat_diff_hdi, this_feat_id, suffix="_hdi"
        )
        return pd.concat(
            [this_feat_diff_mean, this_feat_diff_std, this_feat_diff_hdis], axis=1
        )
    except Exception as e:
        print(f"Error processing file {inf_file}: {str(e)}")  # TODO: chaneg this to log
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
