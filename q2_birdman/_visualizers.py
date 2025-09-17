import importlib.resources
from pathlib import Path
import re
import os
from collections import Counter

import altair as alt
import pandas as pd
import numpy as np
from scipy import stats

import qiime2
import q2templates
import biom

def _sanitize_filename(filename):
    """Convert a string to a safe filename by replacing invalid characters.
    
    Parameters
    ----------
    filename : str
        The string to convert to a safe filename
        
    Returns
    -------
    str
        A safe filename with invalid characters replaced
    """
    # Replace invalid filesystem characters with underscores
    # Windows: < > : " | ? * \ /
    # Unix: / (and null bytes, but we can't have those in strings)
    invalid_chars = r'[<>:"|?*\\/\[\]]'
    safe_name = re.sub(invalid_chars, '_', filename)
    
    # Remove leading/trailing spaces and dots
    safe_name = safe_name.strip(' .')
    
    # Ensure the filename isn't empty
    if not safe_name:
        safe_name = 'unnamed'
    
    return safe_name

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
    base_name = col.replace('_hdi', '')
    clean_mean_col = 'mean_effect'
    df[clean_mean_col] = df[base_name + '_mean']
    df[["lower", "upper"]] = df[col].str.split(",", expand=True)
    df.lower = df.lower.str[1:].astype("float")
    df.upper = df.upper.str[:-1].astype("float")
    df["credible"] = np.where((df.lower > 0) | (df.upper < 0), "yes", "no")
    return df

def _compute_sample_log_ratios(table_df, sub_df, effect_size_label):
    """Compute log ratios of enriched vs depleted features for each sample.
    
    Parameters
    ----------
    table_df : pd.DataFrame
        Feature table DataFrame
    sub_df : pd.DataFrame
        DataFrame containing differential abundance results
    effect_size_label : str
        Label for effect sizes in data
        
    Returns
    -------
    pd.DataFrame
        DataFrame with log ratios for each sample
    """
    enriched_features = sub_df[sub_df[effect_size_label + "_mean"] > 0].index
    depleted_features = sub_df[sub_df[effect_size_label + "_mean"] < 0].index
    result = pd.DataFrame(index=table_df.columns)
    enriched_sums = table_df.loc[enriched_features].sum(axis=0)
    depleted_sums = table_df.loc[depleted_features].sum(axis=0)
    result['log_ratio'] = np.log2((enriched_sums + 1) / (depleted_sums + 1))
    return result

def _perform_kruskal_wallis(data, groups):
    """Perform Kruskal-Wallis test on log ratios between groups.
    
    Parameters
    ----------
    data : pd.Series
        Log ratio values
    groups : pd.Series
        Group labels
        
    Returns
    -------
    tuple
        (test statistic, p-value)
    """
    groups = groups.dropna()
    data = data[groups.index]
    unique_groups = groups.unique()
    
    # check enough data
    if (len(unique_groups) < 2 or len(data) < 3 or data.nunique() <= 1):
        return None, None
    
    group_data = [data[groups == g] for g in unique_groups]
    
    # check if any group has insufficient variation
    if any(len(group) > 0 and group.nunique() <= 1 for group in group_data):
        return None, None
    
    try:
        stat, pval = stats.kruskal(*group_data)
        return stat, pval
    except ValueError as e:
        # hndle cases where Kruskal-Wallis fails due to insufficient variation
        if "All numbers are identical" in str(e):
            return None, None
        else:
            raise e

def _perform_pearson_correlation(x, y):
    """Perform Pearson correlation between two variables.
    
    Parameters
    ----------
    x : pd.Series
        First variable
    y : pd.Series
        Second variable
        
    Returns
    -------
    tuple
        (correlation coefficient, p-value)
    """
    mask = ~(x.isna() | y.isna())
    x_clean = x[mask]
    y_clean = y[mask]
    
    # Check all validation conditions at once
    if (mask.sum() < 2 or 
        x_clean.nunique() <= 1 or 
        y_clean.nunique() <= 1):
        return None, None
    
    try:
        r, pval = stats.pearsonr(x_clean, y_clean)
        return r, pval
    except ValueError as e:
        # handle cases where Pearson correlation fails due to insufficient variation
        if "All numbers are identical" in str(e):
            return None, None
        else:
            raise e

def _create_statistical_annotation(data, is_numeric):
    """Create statistical annotation text for the plot.
    
    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing log ratios and metadata values
    is_numeric : bool
        Whether the metadata is numeric
        
    Returns
    -------
    str
        Formatted statistical annotation
    """
    if is_numeric:
        r, pval = _perform_pearson_correlation(data['metadata_value'], data['log_ratio'])
        if r is not None and pval is not None:
            return f"Pearson\nr = {r:.3f}\np = {pval:.3f}"
    else:
        stat, pval = _perform_kruskal_wallis(data['log_ratio'], data['metadata_value'])
        if stat is not None and pval is not None:
            return f"Kruskal-Wallis\nH = {stat:.3f}\np = {pval:.3f}"
    return None

def _create_metadata_visualization(sub_df, table_df, metadata_df, metadata_cols, effect_size_label, effect_size_threshold=0.0, palette="category10"):
    """Create metadata visualization showing log ratios of enriched vs depleted features.
    
    Parameters
    ----------
    sub_df : pd.DataFrame
        DataFrame containing differential abundance results
    table_df : pd.DataFrame
        Feature table DataFrame
    metadata_df : pd.DataFrame
        Sample metadata DataFrame
    metadata_cols : list
        List of metadata column names
    effect_size_label : str
        Label for effect sizes in data
    effect_size_threshold : float, optional
        Minimum absolute effect size to include, by default 0.0
    palette : str
        Color scheme for metadata categories
    """
    # Filter features based on effect size threshold
    sub_df = sub_df[np.abs(sub_df[effect_size_label + '_mean']) >= effect_size_threshold]
    
    # Check all validation conditions at once
    if (len(sub_df) == 0 or len(sub_df) < 2):
        return None
        
    log_ratios = _compute_sample_log_ratios(table_df, sub_df, effect_size_label)
    log_ratios['log_ratio'] = log_ratios['log_ratio'].replace([np.inf, -np.inf], np.nan)
    log_ratios = log_ratios.dropna()
    
    # Check remaining validation conditions
    if (len(log_ratios) < 3 or log_ratios['log_ratio'].nunique() <= 1):
        return None
    
    log_ratio_min = log_ratios['log_ratio'].min()
    log_ratio_max = log_ratios['log_ratio'].max()
    log_ratio_limits = [log_ratio_min * 1.1, log_ratio_max * 1.1]
    
    plot_data = []
    for col in metadata_cols:
        col_data = log_ratios.merge(metadata_df[[col]], left_index=True, right_index=True)
        col_data['metadata_column'] = col
        col_data = col_data.rename(columns={col: 'metadata_value'})
        col_data['is_numeric'] = pd.api.types.is_numeric_dtype(metadata_df[col])
        col_data['stat_annotation'] = _create_statistical_annotation(col_data, col_data['is_numeric'].iloc[0])
        col_data['mean_log_ratio'] = col_data.groupby('metadata_value')['log_ratio'].transform('mean')
        plot_data.append(col_data)
    
    plot_data = pd.concat(plot_data, ignore_index=True)
    
    dropdown = alt.binding_select(
        options=metadata_cols,
        name="Select Metadata Column: "
    )
    
    selection = alt.param(
        name='Column',
        value=metadata_cols[0],
        bind=dropdown
    )
    
    base = alt.Chart(plot_data).encode(
        y=alt.Y('log_ratio:Q', 
                title='Log2(Enriched/Depleted)',
                axis=alt.Axis(grid=True),
                scale=alt.Scale(domain=log_ratio_limits)),
        tooltip=['metadata_column', 'metadata_value', 'log_ratio']
    ).transform_filter(
        alt.datum.metadata_column == selection
    )
    
    scatter = base.mark_circle(size=60).encode(
        x=alt.X('metadata_value:Q', 
                axis=alt.Axis(orient='bottom', title=None),
                sort=alt.SortField('mean_log_ratio', order='descending')),
        color=alt.Color('metadata_value:Q', 
                       scale=alt.Scale(scheme=palette),
                       legend=alt.Legend(title="Value"))
    ).transform_filter(
        'datum.is_numeric == true'
    ) + base.mark_circle(size=60).encode(
        x=alt.X('metadata_value:N', 
                axis=alt.Axis(orient='bottom', labelAngle=30, title=None),
                sort=alt.SortField('mean_log_ratio', order='descending')),
        color=alt.Color('metadata_value:N',
                       scale=alt.Scale(scheme=palette),
                       legend=alt.Legend(title="Category"))
    ).transform_filter(
        'datum.is_numeric == false'
    )
    
    regression = alt.Chart(plot_data).transform_filter(
        alt.datum.metadata_column == selection
    ).transform_filter(
        'datum.is_numeric == true'
    ).mark_line(
        color='gray',
        opacity=0.8,
        strokeWidth=2,
        strokeDash=[5, 5]
    ).encode(
        x=alt.X('metadata_value:Q'),
        y=alt.Y('log_ratio:Q')
    ).transform_regression(
        'metadata_value', 'log_ratio',
        method='quad'
    )
    
    annotation = alt.Chart(plot_data).mark_text(
        align='left',
        baseline='top',
        fontSize=12,
        fontWeight='lighter',
        lineHeight=20
    ).encode(
        x=alt.value(5),
        y=alt.value(5),
        text=alt.Text('stat_annotation:N')
    ).transform_filter(
        alt.datum.metadata_column == selection
    )
    
    plot = (scatter + regression + annotation).resolve_scale(
        x='independent',
        color='independent'
    )
    
    return plot.add_params(
        selection
    ).properties(
        title='Log-Ratios of Enriched/Depleted Features',
        width=500,
        height=400
    ).configure_axis(
        labelLimit=100
    ).configure_axisBottom(
        labelAngle=45,
        minExtent=50
    )

def _plot_differentials(
        output_dir,
        data,
        taxonomy,
        title,
        effect_size_label,
        effect_size_threshold,
        taxonomy_delimiter,
        label_limit,
        chart_style="bar",
        palette="category10"):
    """Create differential abundance plots using Altair.
    
    Parameters
    ----------
    output_dir : str
        Directory to save the plot
    data : qiime2.Metadata
        Data to plot
    taxonomy : pd.DataFrame
        Taxonomy information to annotate features
    title : str
        Plot title
    effect_size_label : str
        Label for effect sizes
    effect_size_threshold : float
        Minimum absolute effect size to include
    taxonomy_delimiter : str
        Delimiter used in taxonomy strings to split taxonomic levels
    label_limit : int
        Maximum length for axis labels
    chart_style : str
        Plot style, either "bar" or "forest"
    palette : str
        Color scheme for enriched/depleted features. Can be a discrete Altair scheme
        (e.g., "category10", "accent", "dark2", "paired", "set1", "set2", "set3", "tableau10", "tableau20")
        or a comma-separated pair of hex colors (e.g., "#4c78a8,#f58518")
    """
    if len(data.to_dataframe()) == 0:
        raise ValueError("No features present in input.")

    md_df = data.to_dataframe()
    sub_df = _unpack_hdi_and_filter(md_df, effect_size_label + "_hdi")
    sub_df.rename_axis(index="Feature", inplace=True)
    sub_df = sub_df.sort_values(by='mean_effect')

    sub_df = sub_df[np.abs(sub_df['mean_effect']) >= effect_size_threshold]
    
    filter_credible = True
    if filter_credible:
        sub_df = sub_df[sub_df.credible == "yes"]

    # Set up output directory and file path
    output_dir = Path(output_dir)
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        raise ValueError(f"Failed to create output directory {output_dir}: {e}")
    
    safe_title = _sanitize_filename(title)
    fig_fn = Path(f'{safe_title}-birdman-{chart_style}-plot.html')
    fig_fp = output_dir / fig_fn
    
    # Ensure the output directory is writable
    if not os.access(output_dir, os.W_OK):
        raise ValueError(f"Output directory {output_dir} is not writable")

    if len(sub_df) == 0:
        # Create an empty chart with a message
        empty_chart = alt.Chart(pd.DataFrame({'x': [0], 'y': [0]})).mark_text(
            text=f"No features remaining after applying filters (effect size threshold: {effect_size_threshold})",
            fontSize=14,
            color='gray'
        ).encode(
            x=alt.value(250),  # Center of the chart
            y=alt.value(200)   # Center of the chart
        ).properties(
            width=500,
            height=400,
            title=title
        )
        try:
            empty_chart.save(fig_fp)
            return fig_fp
        except Exception as e:
            raise ValueError(f"Failed to save empty chart to {fig_fp}: {e}")

    if taxonomy is not None and taxonomy_delimiter is not None:
        y_labels = []
        seen = Counter()
        taxonomy_df = taxonomy
        for i, e in enumerate(sub_df.index):
            if e in taxonomy_df.index:
                tax_string = taxonomy_df.loc[e].iloc[0]
                if taxonomy_delimiter in tax_string:
                    fields = [field for field in tax_string.split(taxonomy_delimiter)
                              if not field.endswith('__')]
                else:
                    fields = [tax_string]
                most_specific = fields[-1]
                if most_specific in seen:
                    y_labels.append(f"{seen[most_specific]}: {most_specific} *")
                else:
                    y_labels.append(most_specific)
                seen[most_specific] += 1
            else:
                y_labels.append(e)
        sub_df['y_label'] = y_labels
        sub_df['feature'] = [id_.replace(taxonomy_delimiter, ' ')
                         for id_ in sub_df.index]
    else:
        sub_df['y_label'] = sub_df['feature'] = sub_df.index

    sub_df['enriched'] = ["enriched" if x else "depleted"
                      for x in sub_df['mean_effect'] > 0]

    sub_df['error-upper'] = sub_df['upper']
    sub_df['error-lower'] = sub_df['lower']

    shared_y = alt.Y("y_label",
                     title="Feature ID",
                     sort=alt.SortField('mean_effect', order='descending' if chart_style == "bar" else 'ascending'))

    if chart_style == "bar":
        bars = alt.Chart(sub_df).mark_bar().encode(
            x=alt.X('mean_effect', title="Log Fold Change (LFC)"),
            y=shared_y,
            tooltip=alt.Tooltip(["feature", 'mean_effect',
                                 "error-lower", "error-upper"]),
            color=alt.Color('enriched', title="Relative to reference",
                            scale=alt.Scale(scheme=palette))
        )

        error = alt.Chart(sub_df).mark_rule(color='black').encode(
            x='error-lower',
            x2='error-upper',
            y=shared_y,
        )

        chart = (bars + error).properties(title=title)
    else:
        error_bars = alt.Chart(sub_df).mark_errorbar(extent='ci').encode(
            x=alt.X('error-lower', title="Log Fold Change (LFC)"),
            x2='error-upper',
            y=shared_y,
            color=alt.Color('enriched', title="Relative to reference",
                            scale=alt.Scale(scheme=palette, reverse=True),
                            sort="descending")
        )

        points = alt.Chart(sub_df).mark_point(
            filled=True,
            stroke='black',
            strokeWidth=1
        ).encode(
            x=alt.X('mean_effect'),
            y=shared_y,
            color=alt.Color('enriched', title="Relative to reference",
                            scale=alt.Scale(scheme=palette, reverse=True),
                            sort="descending"),
            tooltip=alt.Tooltip(["feature", 'mean_effect',
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

    if chart_style == "forest":
        sub_df = sub_df.sort_values(by='mean_effect', ascending=True)
    else:
        sub_df = sub_df.sort_values(by='mean_effect', ascending=False)

    try:
        chart.save(fig_fp)
        return fig_fp
    except Exception as e:
        raise ValueError(f"Failed to save chart to {fig_fp}: {e}")


def plot(output_dir: str,
            data: qiime2.Metadata,
            taxonomy: pd.DataFrame = None,
            table: biom.Table = None,
            metadata: pd.DataFrame = None,
            effect_size_threshold: float = 0.0,
            taxonomy_delimiter: str = ';',
            label_limit: int = None,
            chart_style: str = "bar",
            palette: str = "category10") -> None:
    """Generate bar plot views of differential abundance analysis output.
    
    Parameters
    ----------
    output_dir : str
        Directory where visualization outputs will be written
    data : qiime2.Metadata
        The differential abundance analysis output to be plotted
    taxonomy : pd.DataFrame, optional
        Taxonomy information to annotate features
    table : biom.Table, optional
        The feature table containing the samples over which feature-based differential abundance was computed
    metadata : pd.DataFrame, optional
        The sample metadata that includes the columns used in the analysis
    effect_size_label : str, optional
        Label for effect sizes in data, by default 'lfc'
    effect_size_threshold : float, optional
        Exclude features with an absolute value of effect size less than this
        threshold, by default 0.0
    taxonomy_delimiter : str, optional
        Delimiter used in taxonomy strings to split taxonomic levels, by default ';'
    label_limit : int, optional
        Set the maximum length that will be viewable for axis labels, by default None
    chart_style : str, optional
        Style of the plot, either "bar" or "forest", by default "bar"
    palette : str, optional
        Color scheme for enriched/depleted features, by default "category10"
    """
    
    df = data.to_dataframe()
    
    # Find all columns that end with _hdi (these are our effect size columns) and exclude the Intercept column
    effect_size_columns = [col for col in df.columns if col.endswith('_hdi') and not col.startswith('Intercept')]
    
    if not effect_size_columns:
        raise ValueError("No effect size columns found in input data. Expected columns ending with '_hdi'.")
    
    # Check if effect size threshold might be too restrictive
    if effect_size_threshold > 0:
        total_features = len(df)
        for hdi_col in effect_size_columns:
            base_name = hdi_col.replace('_hdi', '')
            filtered_count = len(df[np.abs(df[base_name + '_mean']) >= effect_size_threshold])
            if filtered_count == 0:
                print(f"WARNING: Effect size threshold ({effect_size_threshold}) for {base_name} filters out all {total_features} features. Consider using a lower threshold.")
            elif filtered_count < 5:
                print(f"WARNING: Effect size threshold ({effect_size_threshold}) for {base_name} leaves only {filtered_count}/{total_features} features. This may be too restrictive for meaningful analysis.")
    
    # Create lists to store figure data
    da_figure_data = []
    metadata_figure_data = []
    
    # Create plots for each effect size column
    for hdi_col in effect_size_columns:
        # Extract the base name (e.g., 'lfc' from 'lfc_hdi')
        base_name = hdi_col.replace('_hdi', '')
        
        try:
            # Create differential abundance plot
            da_figure_fp = _plot_differentials(
                output_dir, data, taxonomy,
                title=base_name,
                effect_size_label=base_name,
                effect_size_threshold=effect_size_threshold,
                taxonomy_delimiter=taxonomy_delimiter,
                label_limit=label_limit,
                chart_style=chart_style,
                palette=palette)
            da_figure_fn = da_figure_fp.parts[-1]
            da_figure_data.append((True, da_figure_fn, base_name, None))
            
            # Create metadata visualization if table and metadata are provided
            if table is not None and metadata is not None:
                print("Creating metadata visualization...")
                table_df = table.to_dataframe()
                    
                # Get metadata columns for dropdown
                metadata_df = metadata.to_dataframe()
                metadata_cols = metadata_df.columns.tolist()
                print(f"Available metadata columns: {metadata_cols}")
                
                # Create metadata visualization
                print("Creating metadata chart...")
                metadata_chart = _create_metadata_visualization(
                    df, table_df, metadata_df, metadata_cols, base_name, 
                    effect_size_threshold=effect_size_threshold, palette=palette)
                if metadata_chart is not None:
                    print("Metadata chart created successfully")
                    safe_base_name = _sanitize_filename(base_name)
                    metadata_figure_fp = Path(output_dir) / f"{safe_base_name}_metadata.html"
                    print(f"Saving metadata chart to {metadata_figure_fp}")
                    try:
                        metadata_chart.save(metadata_figure_fp)
                        metadata_figure_fn = metadata_figure_fp.parts[-1]
                        metadata_figure_data.append((True, metadata_figure_fn, base_name, None))
                    except Exception as e:
                        print(f"Failed to save metadata chart: {e}")
                        metadata_figure_data.append((False, None, base_name, f"Failed to save: {e}"))
                else:
                    # Provide more specific error message
                    filtered_df = df[np.abs(df[base_name + '_mean']) >= effect_size_threshold]
                    if len(filtered_df) == 0:
                        error_msg = f"Effect size threshold ({effect_size_threshold}) too restrictive - no features remain"
                    elif len(filtered_df) < 2:
                        error_msg = f"Effect size threshold ({effect_size_threshold}) too restrictive - only {len(filtered_df)} feature(s) remain"
                    else:
                        error_msg = "Insufficient variation in log ratios for statistical analysis"
                    print(f"Failed to create metadata chart: {error_msg}")
                    metadata_figure_data.append((False, None, base_name, error_msg))
            else:
                print("No table or metadata provided for visualization")
                metadata_figure_data.append((False, None, base_name, "No metadata or table provided"))
                
        except ValueError as e:
            da_figure_data.append((False, None, base_name, str(e)))
            metadata_figure_data.append((False, None, base_name, str(e)))
            raise ValueError(f"Error creating differential abundance plot for {base_name}: {str(e)}")

    
    # Create index.html
    ASSETS = importlib.resources.files('q2_birdman') / 'assets'
    index = Path(ASSETS, 'diff_abundance_plots', 'index.html')
    
    context = {
        'da_figures': da_figure_data,
        'metadata_figures': metadata_figure_data
    }
    q2templates.render(str(index), output_dir, context=context) 