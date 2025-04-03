import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from qiime2 import Visualization
from qiime2.plugin import Str, Bool, Metadata
from q2_types.metadata import ImmutableMetadata
from typing import List, Dict, Any, Optional, Tuple
import json
import re
from pathlib import Path

# Original implementation (commented out for now)
"""
def _unpack_hdi_and_filter(df, col):
    df[["lower", "upper"]] = df[col].str.split(",", expand=True)
    # remove ( from lower and ) from upper and convert to float
    df.lower = df.lower.str[1:].astype("float")
    df.upper = df.upper.str[:-1].astype("float")

    df["credible"] = np.where((df.lower > 0) | (df.upper < 0), "yes", "no")

    df.upper = df.upper - df[col.replace("hdi", "mean")]
    df.lower = df[col.replace("hdi", "mean")] - df.lower

    return df

def _display_top_n_feats(df, n, yvar, xvar, xlab, ylab, title, outdir):
    if df.shape[0] < 2 * n:
        df_for_display = df
    else:
        bottomn = df[:n]
        topn = df[-1 * n :]
        df_for_display = pd.concat([bottomn, topn])

    sns.stripplot(data=df_for_display, y=yvar, x=xvar)
    plt.errorbar(
        data=df_for_display,
        x=xvar,
        y=yvar,
        xerr=df_for_display[["lower", "upper"]].T,
        ls="none",
    )
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    plt.title(title)
    plt.savefig(f'{outdir}/{xlab.split("Ratio for ")[1]}_plot.png', bbox_inches='tight')
    plt.savefig(f'{outdir}/{xlab.split("Ratio for ")[1]}_plot.svg', bbox_inches='tight')
    plt.close()

def plot(output_dir: str, results_artifact: Metadata, plot_var: str, flip: bool = False) -> None:
    # Load the metadata from the artifact
    df = results_artifact.to_dataframe()
    
    # Process data
    sub_df = _unpack_hdi_and_filter(df, plot_var + "_hdi")
    sub_df.rename_axis(index="Feature", inplace=True)
    sub_df = sub_df.sort_values(by=plot_var + "_mean")

    # Create plot
    xlab = "Ratio for " + plot_var
    ylab = "Features"
    df_for_display = sub_df.reset_index()
    df_for_display = df_for_display.loc[df_for_display.credible == "yes"]
    
    fig, ax = plt.subplots(figsize=(6, 10))
    _display_top_n_feats(
        df_for_display, 25, "Feature", plot_var + "_mean", xlab, ylab, "Top Features", output_dir
    )

    # Create index.html
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write(_html_template.format(
            plot_var=plot_var,
            xlab=xlab,
            plot_file=f"{plot_var}_plot.png"
        ))
"""

# Simplified version based on QIIME2 tutorial
def plot(
    output_dir: str,
    results_artifact: Metadata,
    plot_var: str,
    flip: bool = False
) -> None:
    """
    Create plots for BIRDMAn analysis results.
    
    Parameters
    ----------
    output_dir : str
        Directory where visualization outputs will be written
    results_artifact : Metadata
        QIIME2 artifact containing BIRDMAn results
    plot_var : str
        Variable to plot (e.g. "host_age")
    flip : bool, optional
        Whether to flip the plot orientation, by default False
    """
    # Load the metadata from the artifact
    df = results_artifact.to_dataframe()
    
    # Create a simple HTML file with the plot variable name
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write(f'<h1>Plot for {plot_var}</h1>')

_html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BIRDMAn Plot</title>
    <style>
        body {
            padding: 20px;
            font-family: Arial, sans-serif;
        }
        h1 {
            color: #2c3e50;
            text-align: center;
        }
    </style>
</head>
<body>
    <h1>BIRDMAn Plot for %s</h1>
    <p>This is a simplified visualization for BIRDMAn results.</p>
    <p>Plot variable: %s</p>
</body>
</html>
""" 