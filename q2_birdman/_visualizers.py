import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from qiime2 import Visualization, Artifact, Metadata
from q2_types.metadata import ImmutableMetadata
from q2_types.feature_data import FeatureData, Taxonomy
from typing import List, Dict, Any, Optional, Tuple
import json
import re
from pathlib import Path

_html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BIRDMAn Plot</title>
    <style>
        body {{
            padding: 20px;
            font-family: Arial, sans-serif;
        }}
        h1 {{
            color: #2c3e50;
            text-align: center;
        }}
    </style>
</head>
<body>
    <img src="{plot_file}" alt="Plot for {plot_var}">
</body>
</html>
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

def _parse_taxon(taxon):
    levels = taxon.split(';')
    genus = levels[-2].strip() if len(levels) > 1 else ""
    species = levels[-1].strip() if len(levels) > 0 else ""
    
    # Replace placeholder values (e.g., 'g__', 's__') with an empty string
    genus = genus if not genus.endswith("__") else ""
    species = species if not species.endswith("__") else ""
    
    return genus, species

def _display_top_n_feats(df, n, yvar, xvar, xlab, ylab, title, outdir, invert=False):
    if df.shape[0] < 2 * n:
        df_for_display = df
    else:
        bottomn = df[:n]
        topn = df[-1 * n :]
        df_for_display = pd.concat([bottomn, topn])

    if invert:
        df_for_display = df_for_display.iloc[::-1]
        df_for_display[xvar] = df_for_display[xvar] * -1
        df_for_display[['lower', 'upper']] = df_for_display[['upper', 'lower']]

    # Check for duplicated labels
    df_for_display[yvar] = df_for_display[yvar].astype(str)
    duplicates = df_for_display.duplicated(subset=[yvar], keep=False)
    df_for_display.loc[duplicates, yvar] += df_for_display[duplicates].groupby(yvar).cumcount().add(1).astype(str).radd(' #')

    fig, ax = plt.subplots(figsize=(6, 10))
    sns.stripplot(data=df_for_display, y=yvar, x=xvar, size=5, edgecolor='gray', linewidth=0.5, palette=['black'])
    plt.errorbar(
        x=df_for_display[xvar],
        y=np.arange(df_for_display.shape[0]),
        xerr=df_for_display[['lower', 'upper']].T.values,
        ls='none',
        ecolor='gray'
    )
    
    ax.set_ylabel(ylab, fontsize=14, fontweight='normal')
    ax.set_xlabel(xlab, fontsize=14, fontweight='normal')
    ax.set_title(title, fontsize=16, fontweight='bold')
    
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    
    plt.tight_layout()
    plt.savefig(f'{outdir}/{xlab.split("Ratio for ")[1]}_plot.png', bbox_inches='tight', dpi=300)
    plt.savefig(f'{outdir}/{xlab.split("Ratio for ")[1]}_plot.svg', bbox_inches='tight', dpi=300)
    plt.close()

def plot(output_dir: str, results_artifact: Metadata, plot_var: str, flip: bool = False, taxonomy: pd.DataFrame = None) -> None:
    """
    Create plots for BIRDMAn analysis results.
    
    Parameters
    ----------
    output_dir : str
        Directory where visualization outputs will be written
    results_artifact : ImmutableMetadata
        QIIME2 artifact containing BIRDMAn results
    plot_var : str
        Variable to plot (e.g. "host_age")
    flip : bool, optional
        Whether to flip the plot orientation, by default False
    taxonomy : pd.DataFrame, optional
        Optional taxonomy information to annotate features
    """
    # Load the metadata from the artifact
    df = results_artifact.to_dataframe()
    
    # Process data
    sub_df = _unpack_hdi_and_filter(df, plot_var + "_hdi")
    sub_df.rename_axis(index="Feature", inplace=True)
    sub_df = sub_df.sort_values(by=plot_var + "_mean")

    # Add taxonomy information if provided
    if taxonomy is not None:
        sub_df['taxon'] = taxonomy.loc[sub_df.index]['Taxon']
        sub_df['Genus'], sub_df['Species'] = zip(*sub_df['taxon'].apply(_parse_taxon))
        yvar = 'Species'  # Use species as the y-axis label
    else:
        yvar = 'Feature'  # Use feature ID as the y-axis label

    # Create plot
    xlab = "Ratio for " + plot_var
    ylab = f"{plot_var} Feature"
    df_for_display = sub_df.reset_index()
    df_for_display = df_for_display.loc[df_for_display.credible == "yes"]
    
    _display_top_n_feats(
        df_for_display, 25, yvar, plot_var + "_mean", xlab, ylab, "Top Features", output_dir, flip
    )

    # Create index.html
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write(_html_template.format(
            plot_var=plot_var,
            xlab=xlab,
            plot_file=f"{plot_var}_plot.png"
        ))