import importlib.resources
from pathlib import Path
import urllib.parse
from collections import Counter

import altair as alt
import pandas as pd
import numpy as np

import qiime2
import q2templates


def _unpack_hdi_and_filter(df, col):
    """Unpack high density intervals from string format and parse credible intervals.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing HDI values
    col : str
        Column name containing HDI values in format "(lower,upper)"
        
    Returns
    -------
    pd.DataFrame
        DataFrame with unpacked HDI values and credible interval calculations
    """
    df[["lower", "upper"]] = df[col].str.split(",", expand=True)
    
    # remove ( from lower and ) from upper and convert to float
    df.lower = df.lower.str[1:].astype("float")
    df.upper = df.upper.str[:-1].astype("float")

    df["credible"] = np.where((df.lower > 0) | (df.upper < 0), "yes", "no")

    # Keep the original HDI values for error bars
    # No need to subtract mean as Altair expects the actual range values
    return df

def _plot_differentials(
        output_dir,
        data,
        taxonomy,
        title,
        feature_id_label,
        effect_size_label,
        error_label,
        feature_ids,
        effect_size_threshold,
        taxonomy_delimiter,
        label_limit,
        chart_style):
    """Create differential abundance plots using Altair.
    
    Parameters
    ----------
    output_dir : str
        Directory to save the plot
    data : qiime2.Metadata
        Data to plot
    title : str
        Plot title
    feature_id_label : str
        Label for feature IDs
    effect_size_label : str
        Label for effect sizes
    error_label : str
        Label for error values
    feature_ids : qiime2.Metadata
        Optional feature IDs to include
    effect_size_threshold : float
        Minimum absolute effect size to include
    taxonomy_delimiter : str
        Delimiter used in taxonomy strings to split taxonomic levels
    label_limit : int
        Maximum length for axis labels
    chart_style : str
        Plot style, either "bar" or "forest"
        
    Returns
    -------
    Path
        Path to the saved plot file
    """
    if len(data.to_dataframe()) == 0:
        raise ValueError("No features present in input.")

    # Process data similar to _visualizers.py
    md_df = data.to_dataframe()
    sub_df = _unpack_hdi_and_filter(md_df, effect_size_label + "_hdi")
    sub_df.rename_axis(index="Feature", inplace=True)
    sub_df = sub_df.sort_values(by=effect_size_label + "_mean")

    if feature_ids is not None:
        feature_ids_df = feature_ids.to_dataframe()
        sub_df = sub_df[sub_df.index.isin(feature_ids_df.index)]

    sub_df = sub_df[np.abs(sub_df[effect_size_label + "_mean"]) >= effect_size_threshold]
    
    filter_credible = True
    if filter_credible:
        sub_df = sub_df[sub_df.credible == "yes"]

    if len(sub_df) == 0:
        raise ValueError("No features remaining after applying filters.")

    # Sort data by effect size (ascending for forest plot, descending for bar plot)
    if chart_style == "forest":
        sub_df = sub_df.sort_values(by=effect_size_label + "_mean", ascending=True)
    else:
        sub_df = sub_df.sort_values(by=effect_size_label + "_mean", ascending=False)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    safe_title = urllib.parse.quote(title)
    fig_fn = Path(f'{safe_title}-birdman-{chart_style}-plot.html')
    fig_fp = output_dir / fig_fn

    if taxonomy and taxonomy_delimiter is not None:
        y_labels = []
        seen = Counter()
        for i, e in enumerate(sub_df.index):
            if taxonomy_delimiter in e:
                fields = [field for field in e.split(taxonomy_delimiter)
                          if not field.endswith('__')]
            else:
                fields = [e]
            most_specific = fields[-1]
            if most_specific in seen:
                y_labels.append(f"{seen[most_specific]}: {most_specific} *")
            else:
                y_labels.append(most_specific)
            seen[most_specific] += 1
        sub_df['y_label'] = y_labels
        sub_df['feature'] = [id_.replace(taxonomy_delimiter, ' ')
                         for id_ in sub_df.index]
    else:
        sub_df['y_label'] = sub_df['feature'] = sub_df.index

    sub_df['enriched'] = ["enriched" if x else "depleted"
                      for x in sub_df[effect_size_label + "_mean"] > 0]

    # Use the HDI values directly for error bars
    sub_df['error-upper'] = sub_df['upper']
    sub_df['error-lower'] = sub_df['lower']

    shared_y = alt.Y("y_label",
                     title="Feature ID",
                     sort=None)  # No need to sort in Altair since data is pre-sorted

    if chart_style == "bar":
        # Bar chart with error bars
        bars = alt.Chart(sub_df).mark_bar().encode(
            x=alt.X(effect_size_label + "_mean", title="Log Fold Change (LFC)"),
            y=shared_y,
            tooltip=alt.Tooltip(["feature", effect_size_label + "_mean",
                                 "error-lower", "error-upper"]),
            color=alt.Color('enriched', title="Relative to reference",
                            scale=alt.Scale(
                                domain=["enriched", "depleted"],
                                range=["#4c78a8", "#f58518"]),
                            sort="descending")
        )

        error = alt.Chart(sub_df).mark_rule(color='black').encode(
            x='error-lower',
            x2='error-upper',
            y=shared_y,
        )

        chart = (bars + error).properties(title=title)
    else:
        # Forest plot style
        error_bars = alt.Chart(sub_df).mark_errorbar(extent='ci').encode(
            x=alt.X('error-lower', title="Log Fold Change (LFC)"),
            x2='error-upper',
            y=shared_y,
            color=alt.Color('enriched', title="Relative to reference",
                            scale=alt.Scale(
                                domain=["enriched", "depleted"],
                                range=["#4c78a8", "#f58518"]),
                            sort="descending")
        )

        points = alt.Chart(sub_df).mark_point(
            filled=True,
            stroke='black',
            strokeWidth=1
        ).encode(
            x=alt.X(effect_size_label + "_mean"),
            y=shared_y,
            color=alt.Color('enriched', title="Relative to reference",
                            scale=alt.Scale(
                                domain=["enriched", "depleted"],
                                range=["#4c78a8", "#f58518"]),
                            sort="descending"),
            tooltip=alt.Tooltip(["feature", effect_size_label + "_mean",
                                 "error-lower", "error-upper"])
        )

        chart = (error_bars + points).properties(title=title)

    chart = chart.configure_legend(
        strokeColor='gray',
        padding=10,
        cornerRadius=10,
    )

    chart = chart.configure_axisY(titleAlign='left',
                                  titleY=-10, titleAngle=0)
    if label_limit is not None:
        chart = chart.configure_axis(labelLimit=label_limit)
    else:
        chart = chart.configure_axis(labelLimit=0)

    chart.save(fig_fp)
    return fig_fp


def da_plot(output_dir: str,
            data: qiime2.Metadata,
            taxonomy: pd.DataFrame = None,
            effect_size_label: str = 'lfc',
            feature_id_label: str = 'id',
            error_label: str = 'se',
            effect_size_threshold: float = 0.0,
            feature_ids: qiime2.Metadata = None,
            taxonomy_delimiter: str = None,
            label_limit: int = None,
            chart_style: str = "bar") -> None:
    """Generate bar plot views of differential abundance analysis output.
    
    Parameters
    ----------
    output_dir : str
        Directory where visualization outputs will be written
    data : qiime2.Metadata
        The differential abundance analysis output to be plotted
    taxonomy : pd.DataFrame, optional
        Optional taxonomy information to annotate features
    effect_size_label : str, optional
        Label for effect sizes in data, by default 'lfc'
    feature_id_label : str, optional
        Label for feature ids in data, by default 'id'
    error_label : str, optional
        Label for effect size errors in data, by default 'se'
    effect_size_threshold : float, optional
        Exclude features with an absolute value of effect size less than this
        threshold, by default 0.0
    feature_ids : qiime2.Metadata, optional
        Exclude features if their ids are not included in this index, by default None
    taxonomy_delimiter : str, optional
        Delimiter used in taxonomy strings to split taxonomic levels, by default None
    label_limit : int, optional
        Set the maximum length that will be viewable for axis labels, by default None
    chart_style : str, optional
        Style of the plot, either "bar" or "forest", by default "bar"
    """
    
    df = data.to_dataframe()
    
    # Find all columns that end with _hdi (these are our effect size columns) and exclude the Intercept column
    effect_size_columns = [col for col in df.columns if col.endswith('_hdi') and not col.startswith('Intercept')]
    
    if not effect_size_columns:
        raise ValueError("No effect size columns found in input data. Expected columns ending with '_hdi'.")
    
    # Create a list to store figure data
    figure_data = []
    
    # Create plots for each effect size column
    for hdi_col in effect_size_columns:
        # Extract the base name (e.g., 'lfc' from 'lfc_hdi')
        base_name = hdi_col.replace('_hdi', '')
        
        try:
            figure_fp = _plot_differentials(
                output_dir, data, taxonomy,
                title=base_name,
                effect_size_label=base_name,
                feature_id_label=feature_id_label,
                error_label=error_label,
                effect_size_threshold=effect_size_threshold,
                feature_ids=feature_ids,
                taxonomy_delimiter=taxonomy_delimiter,
                label_limit=label_limit,
                chart_style="forest")
            figure_fn = figure_fp.parts[-1]
            figure_data.append((True, figure_fn, base_name, None))
        except ValueError as e:
            figure_data.append((False, None, base_name, str(e)))
    
    # Create index.html
    ASSETS = importlib.resources.files('q2_birdman') / 'assets'
    index = Path(ASSETS, 'diff_abundance_plots', 'index.html')
    
    context = {
        'figures': figure_data
    }
    q2templates.render(str(index), output_dir, context=context) 