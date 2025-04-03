import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from qiime2 import Visualization
from qiime2.plugin import Str, Bool

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

def plot(results_artifact: ImmutableMetadata, results_dir: str, plot_var: str, flip: bool = False) -> None:
    """Create plots for BIRDMAn analysis results.
    
    Parameters
    ----------
    output_dir : str
        Directory where visualization outputs will be written
    results_dir : str
        Directory containing BIRDMAn analysis results
    plot_var : str
        Variable to plot (e.g. "host_age_)
    flip : bool, optional
        Whether to flip the plot orientation, by default False
    """
    # Read results
    input_path = os.path.join(results_dir, "results", "beta_var.tsv")
    df = pd.read_csv(input_path, sep="\t", index_col="Feature")
    
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

_html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BIRDMAn Plot for {plot_var}</title>
    <style>
        body {{
            padding: 20px;
            font-family: Arial, sans-serif;
        }}
        .plot-container {{
            text-align: center;
            margin: 20px 0;
        }}
        .plot-container img {{
            max-width: 100%;
            height: auto;
        }}
        h1 {{
            color: #2c3e50;
            text-align: center;
        }}
    </style>
</head>
<body>
    <h1>BIRDMAn Plot for {xlab}</h1>
    <div class="plot-container">
        <img src="{plot_file}" alt="BIRDMAn plot">
    </div>
</body>
</html>
""" 