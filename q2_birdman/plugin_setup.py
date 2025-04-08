# ----------------------------------------------------------------------------
# Copyright (c) 2024, Lucas Patel.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib
from qiime2.plugin import Citations, Plugin, Str, Int, Metadata, Bool, Float, Range
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.metadata import ImmutableMetadata
from q2_types.feature_data import FeatureData, Taxonomy
from q2_birdman import __version__
from q2_birdman._methods import run
from q2_birdman._visualizers import plot
from q2_birdman._diff_abundance_plots import da_plot

citations = Citations.load("citations.bib", package="q2_birdman")

plugin = Plugin(
    name="birdman",
    version=__version__,
    website="https://github.com/biocore/BIRDMAn",
    package="q2_birdman",
    description="Bayesian Inferential Regression for Differential Microbiome Analysis (BIRDMAn) is a framework for performing differential abundance analysis through Bayesian inference.",
    short_description="Bayesian Inferential Regression for Differential Microbiome Analysis",
    # Please retain the plugin-level citation of 'Caporaso-Bolyen-2024'
    # as attribution of the use of this template, in addition to any citations
    # you add.
    citations=[citations['Caporaso-Bolyen-2024']]
)

importlib.import_module('q2_birdman._transformers')

plugin.methods.register_function(
    function=run,
    inputs={
        'table': FeatureTable[Frequency],
    },
    parameters={
        'metadata': Metadata,
        'formula': Str,
        'threads': Int,
        'longitudinal': Bool,
        'subject_column': Str
    },
    parameter_descriptions={
        'metadata': 'The sample metadata that includes the columns specified in the formula.',
        'formula': 'The formula used to define the model. This should be a valid Patsy formula that references columns in the metadata.',
        'threads': 'Number of threads to use for parallel processing. Increasing the number of threads can reduce computation time.',
        'longitudinal': ('Whether to use the longitudinal model with random effects for subjects. '
                        'If True, subject_column must also be specified. [default: False]'),
        'subject_column': ('Column name in metadata containing subject IDs for longitudinal analysis. '
                          'Required only if longitudinal=True. [default: None]')
    },
    outputs=[('output_dir', ImmutableMetadata)],
    input_descriptions={
        'table': 'The feature table containing the samples over which feature-based differential abundance should be computed.',
    },
    output_descriptions={
        'output_dir': 'The resulting inference results, including parameter estimates, derived from the BIRDMAn model.',
    },
    name='Run BIRDMAn',
    description=('Run BIRDMAn on a feature table with a given model formula. '
                'Supports both standard and longitudinal analyses using Negative Binomial models.'),
    citations=[]
)

plugin.visualizers.register_function(
    function=plot,
    inputs={
        'results_artifact': ImmutableMetadata,
        'taxonomy': FeatureData[Taxonomy]
    },
    parameters={
        'plot_var': Str,
        'flip': Bool
    },
    input_descriptions={
        'results_artifact': 'QIIME2 artifact containing BIRDMAn results',
        'taxonomy': 'Optional taxonomy information to annotate features'
    },
    parameter_descriptions={
        'plot_var': 'Variable to plot (e.g. "host_age")',
        'flip': 'Whether to flip the plot orientation [default: False]'
    },
    name='Plot BIRDMAn Results',
    description='Create plots for BIRDMAn analysis results, showing the top features associated with a given variable.',
    citations=[]
)

plugin.visualizers.register_function(
    function=da_plot,
    inputs={
        'data': ImmutableMetadata,
        'taxonomy': FeatureData[Taxonomy]
    },
    parameters={
        'effect_size_label': Str,
        'feature_id_label': Str,
        'error_label': Str,
        'effect_size_threshold': Float % Range(0.0, None, inclusive_start=True),
        'feature_ids': Metadata,
        'taxonomy_delimiter': Str,
        'label_limit': Int,
        'chart_style': Str
    },
    input_descriptions={
        'data': 'The differential abundance analysis output to be plotted',
        'taxonomy': 'Optional taxonomy information to annotate features'
    },
    parameter_descriptions={
        'effect_size_label': 'Label for effect sizes in data [default: lfc]',
        'feature_id_label': 'Label for feature ids in data [default: id]',
        'error_label': 'Label for effect size errors in data [default: se]',
        'effect_size_threshold': 'Exclude features with an absolute value of effect size less than this threshold [default: 0.0]',
        'feature_ids': 'Exclude features if their ids are not included in this index [default: None]',
        'taxonomy_delimiter': 'Delimiter used in taxonomy strings to split taxonomic levels [default: None]',
        'label_limit': 'Set the maximum length that will be viewable for axis labels [default: None]',
        'chart_style': 'Style of the plot, either "bar" or "forest" [default: bar]'
    },
    name='Differential Abundance Plot',
    description='Generate bar plot views of differential abundance analysis output, showing enriched and depleted features with error bars.',
    citations=[]
)

