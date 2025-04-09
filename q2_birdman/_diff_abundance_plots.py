import importlib.resources
from pathlib import Path
import urllib.parse
from collections import Counter

import altair as alt
import pandas as pd
import numpy as np
from scipy import stats

import qiime2
import q2templates
import biom

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
    # Get enriched and depleted features
    enriched_features = sub_df[sub_df[effect_size_label + "_mean"] > 0].index
    depleted_features = sub_df[sub_df[effect_size_label + "_mean"] < 0].index
    
    # Initialize result DataFrame
    result = pd.DataFrame(index=table_df.columns)
    
    # Compute sum of enriched features for each sample
    enriched_sums = table_df.loc[enriched_features].sum(axis=0)
    
    # Compute sum of depleted features for each sample
    depleted_sums = table_df.loc[depleted_features].sum(axis=0)
    
    # Compute log ratio (avoiding division by zero)
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
    if len(unique_groups) < 2:
        return None, None
    
    group_data = [data[groups == g] for g in unique_groups]
    stat, pval = stats.kruskal(*group_data)
    return stat, pval

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
    if mask.sum() < 2:
        return None, None
    r, pval = stats.pearsonr(x[mask], y[mask])
    return r, pval

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

def _create_metadata_visualization(sub_df, table_df, metadata_df, metadata_cols, effect_size_label, palette="category10"):
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
    palette : str
        Color scheme for metadata categories. Can be a discrete Altair scheme
        (e.g., "category10", "accent", "dark2", "paired", "set1", "set2", "set3", "tableau10", "tableau20")
        or a comma-separated pair of hex colors (e.g., "#4c78a8,#f58518")
    """
    # Handle color scheme
    if ',' in palette:
        # Custom hex colors
        try:
            colors = palette.split(',')
            if len(colors) != 2:
                raise ValueError("Custom palette must contain exactly two colors")
            # Validate hex colors
            for color in colors:
                if not color.startswith('#') or len(color) != 7:
                    raise ValueError(f"Invalid hex color: {color}")
        except Exception as e:
            raise ValueError(f"Invalid custom palette format: {e}")
    else:
        # Use Altair's built-in scheme directly
        colors = None  # Will use scheme directly in encoding

    # Compute log ratios for each sample
    log_ratios = _compute_sample_log_ratios(table_df, sub_df, effect_size_label)
    
    # Replace infinite values with NaN and drop them
    log_ratios['log_ratio'] = log_ratios['log_ratio'].replace([np.inf, -np.inf], np.nan)
    log_ratios = log_ratios.dropna()
    
    # Set reasonable limits for log ratios
    log_ratio_min = log_ratios['log_ratio'].min()
    log_ratio_max = log_ratios['log_ratio'].max()
    log_ratio_range = max(abs(log_ratio_min), abs(log_ratio_max))
    log_ratio_limits = [-log_ratio_range * 1.1, log_ratio_range * 1.1]
    
    # Prepare data for visualization
    plot_data = []
    for col in metadata_cols:
        # Create a DataFrame for each metadata column
        col_data = log_ratios.merge(metadata_df[[col]], left_index=True, right_index=True)
        col_data['metadata_column'] = col
        col_data = col_data.rename(columns={col: 'metadata_value'})
        # Store the column type
        col_data['is_numeric'] = pd.api.types.is_numeric_dtype(metadata_df[col])
        # Add statistical annotation
        col_data['stat_annotation'] = _create_statistical_annotation(col_data, col_data['is_numeric'].iloc[0])
        plot_data.append(col_data)
    
    # Combine all data
    plot_data = pd.concat(plot_data, ignore_index=True)
    
    # Create dropdown selection
    dropdown = alt.binding_select(
        options=metadata_cols,
        name="Select Metadata Column: "
    )
    
    # Create selection
    selection = alt.param(
        name='Column',
        value=metadata_cols[0],
        bind=dropdown
    )
    
    # Create base plot with conditional encoding
    base = alt.Chart(plot_data).encode(
        y=alt.Y('log_ratio:Q', 
                title='Log2(Enriched/Depleted)',
                axis=alt.Axis(grid=True),
                scale=alt.Scale(domain=log_ratio_limits)),
        tooltip=['metadata_column', 'metadata_value', 'log_ratio']
    ).transform_filter(
        alt.datum.metadata_column == selection
    )
    
    # Create scatter plot with conditional encoding
    scatter = base.mark_circle(size=60).encode(
        x=alt.X('metadata_value:Q', 
                axis=alt.Axis(orient='bottom', title=None)),
        color=alt.Color('metadata_value:Q', 
                       scale=alt.Scale(scheme=palette),
                       legend=alt.Legend(title="Value"))
    ).transform_filter(
        'datum.is_numeric == true'
    ) + base.mark_circle(size=60).encode(
        x=alt.X('metadata_value:N', 
                axis=alt.Axis(orient='bottom', labelAngle=45, title=None)),
        color=alt.Color('metadata_value:N',
                       scale=alt.Scale(scheme=palette),
                       legend=alt.Legend(title="Category"))
    ).transform_filter(
        'datum.is_numeric == false'
    )
    
    # Add boxplot for categorical data
    boxplot = alt.Chart(plot_data).transform_filter(
        alt.datum.metadata_column == selection
    ).transform_filter(
        'datum.is_numeric == false'
    ).mark_boxplot(
        extent='min-max',  # Show whiskers from min to max
        size=30,  # Width of the box
        opacity=0.3  # Make it semi-transparent
    ).encode(
        x=alt.X('metadata_value:N', 
                axis=alt.Axis(orient='bottom', labelAngle=45, title=None)),
        y=alt.Y('log_ratio:Q'),
        color=alt.Color('metadata_value:N',
                       scale=alt.Scale(scheme=palette),
                       legend=None)  # Hide the boxplot legend since we already have one
    )
    
    # Add regression line for numeric data
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
    
    # Create statistical annotation
    annotation = alt.Chart(plot_data).mark_text(
        align='left',
        baseline='top',
        fontSize=12,
        fontWeight='lighter',
        lineHeight=20
    ).encode(
        x=alt.value(5),  # pixels from left
        y=alt.value(5),  # pixels from top
        text=alt.Text('stat_annotation:N')
    ).transform_filter(
        alt.datum.metadata_column == selection
    )
    
    # Combine all elements
    plot = (boxplot + scatter + regression + annotation).resolve_scale(
        x='independent',
        color='independent'
    )
    
    # Add selection and title
    return plot.add_params(
        selection
    ).properties(
        title='Log-Ratios of Enriched/Depleted Features',
        width=500,
        height=400
    )

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
        chart_style,
        palette="category10"):
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
    palette : str
        Color scheme for enriched/depleted features. Can be a discrete Altair scheme
        (e.g., "category10", "accent", "dark2", "paired", "set1", "set2", "set3", "tableau10", "tableau20")
        or a comma-separated pair of hex colors (e.g., "#4c78a8,#f58518")
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
                            scale=alt.Scale(scheme=palette),
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
                            scale=alt.Scale(scheme=palette),
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
                            scale=alt.Scale(scheme=palette),
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
            table: biom.Table = None,
            metadata: pd.DataFrame = None,
            effect_size_label: str = 'lfc',
            feature_id_label: str = 'id',
            error_label: str = 'se',
            effect_size_threshold: float = 0.0,
            feature_ids: qiime2.Metadata = None,
            taxonomy_delimiter: str = None,
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
        Optional taxonomy information to annotate features
    table : biom.Table, optional
        The feature table containing the samples over which feature-based differential abundance was computed
    metadata : pd.DataFrame, optional
        The sample metadata that includes the columns used in the analysis
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
    palette : str, optional
        Color scheme for enriched/depleted features, by default "category10"
    """
    
    df = data.to_dataframe()
    
    # Find all columns that end with _hdi (these are our effect size columns) and exclude the Intercept column
    effect_size_columns = [col for col in df.columns if col.endswith('_hdi') and not col.startswith('Intercept')]
    
    if not effect_size_columns:
        raise ValueError("No effect size columns found in input data. Expected columns ending with '_hdi'.")
    
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
                feature_id_label=feature_id_label,
                error_label=error_label,
                effect_size_threshold=effect_size_threshold,
                feature_ids=feature_ids,
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
                metadata_chart = _create_metadata_visualization(df, table_df, metadata_df, metadata_cols, base_name, palette=palette)
                if metadata_chart is not None:
                    print("Metadata chart created successfully")
                    metadata_figure_fp = Path(output_dir) / f"{base_name}_metadata.html"
                    print(f"Saving metadata chart to {metadata_figure_fp}")
                    metadata_chart.save(metadata_figure_fp)
                    metadata_figure_fn = metadata_figure_fp.parts[-1]
                    metadata_figure_data.append((True, metadata_figure_fn, base_name, None))
                else:
                    print("Failed to create metadata chart")
                    metadata_figure_data.append((False, None, base_name, "Could not create metadata visualization"))
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