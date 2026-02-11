# q2-birdman

BIRDMAn is a framework for performing differential abundance on microbiome data through a Bayesian lens. `q2-birdman` implements the default NegativeBinomial model for both cross-sectional and longitudinal analyses. For more complex experimental designs, users are encouraged to use [BIRDMAn](https://github.com/biocore/BIRDMAn) directly, which enables utilization of custom models and more detailed inference results.

## Installation

`q2-birdman` requires an existing [QIIME 2](https://docs.qiime2.org) environment. Follow the [QIIME 2 installation instructions](https://docs.qiime2.org) to set one up before proceeding.

Once your QIIME 2 environment is active, install the BIRDMAn dependencies:
```bash
mamba install -c conda-forge biom-format patsy xarray arviz cmdstanpy
pip install birdman==0.1.0
```

Then clone and install the plugin:
```bash
git clone https://github.com/lucaspatel/q2-birdman
cd q2-birdman
pip install -e .
```

## Using q2-BIRDMAn

`q2-birdman` provides two main commands: `run` and `plot`. The `run` command performs Bayesian differential abundance analysis on your feature table, while `plot` visualizes the results.

### Running BIRDMAn

The `run` command requires a feature table (`.qza`), metadata file (`.tsv`), and a formula specifying the model. For example, to analyze the effect of age on microbial abundances:

```bash
qiime birdman run \
  --i-table feature-table.qza \
  --m-metadata-file metadata.tsv \
  --p-formula "age" \
  --o-output-dir results.qza \
  --p-threads 32
```

The formula can include multiple variables and interactions. For longitudinal studies, you can use the `--p-longitudinal` flag and specify a subject column:
```bash
qiime birdman run \
  --i-table feature-table.qza \
  --m-metadata-file metadata.tsv \
  --p-formula "age+sex+bmi_score" \
  --o-output-dir results.qza \
  --p-longitudinal True \
  --p-subject-column "host_subject_id"
```

### Visualizing Results

After running the analysis, you can visualize the results using the `plot` command. This creates an interactive visualization showing the differential abundance of features:
```bash
qiime birdman plot \
  --i-data results.qza \
  --i-table feature-table.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization plot.qzv \
  --p-palette rainbow \
  --p-chart-style forest \
  --p-effect-size-threshold 2
```

The visualization supports various options:
- `--p-palette`: Color scheme for enriched/depleted features
- `--p-chart-style`: Plot style (e.g., "forest" or "bar")
- `--p-effect-size-threshold`: Minimum absolute effect size to include
- `--i-taxonomy`: Optional taxonomy file for feature annotation
- `--p-taxonomy-delimiter`: Delimiter used in taxonomy strings

## Testing `q2-birdman`

After completing the install steps above, confirm that everything is working as expected by running:
```shell
make test
```

You should get a report that tests were run, and you should see that all tests passed and none failed. It's usually ok if some warnings are reported.

If all of the tests pass, you're ready to use the plugin.
Start by making QIIME 2's command line interface aware of `q2-birdman` by running:
```shell
qiime dev refresh-cache
```

You should then see the plugin in the list of available plugins if you run:
```shell
qiime info
```

You should be able to review the help text by running:
```shell
qiime birdman --help
```

Have fun! 😎

## Issues

If you encounter issues with `cmdstanpy`, ensure it was installed from conda-forge (as shown above) rather than from pip. If you installed `cmdstanpy` from pip, you must compile the default Negative Binomial model directly (via Python):
```python
import cmdstanpy
cmdstanpy.CmdStanModel(stan_file="q2_birdman/src/stan/negative_binomial_single.stan")
```

## About

The `q2-birdman` Python package was [created from template](https://develop.qiime2.org/en/latest/plugins/tutorials/create-from-template.html).
To learn more about `q2-birdman`, refer to the [project website](https://github.com/biocore/BIRDMAn).
To learn how to use QIIME 2, refer to the [QIIME 2 User Documentation](https://docs.qiime2.org).
To learn QIIME 2 plugin development, refer to [*Developing with QIIME 2*](https://develop.qiime2.org).

`q2-birdman` is a QIIME 2 community plugin, meaning that it is not necessarily developed and maintained by the developers of QIIME 2.
Please be aware that because community plugins are developed by the QIIME 2 developer community, and not necessarily the QIIME 2 developers themselves, some may not be actively maintained or compatible with current release versions of the QIIME 2 distributions.
More information on development and support for community plugins can be found [here](https://library.qiime2.org).
If you need help with a community plugin, first refer to the [project website](https://github.com/biocore/BIRDMAn).
If that page doesn't provide information on how to get help, or you need additional help, head to the [Community Plugins category](https://forum.qiime2.org/c/community-contributions/community-plugins/14) on the QIIME 2 Forum where the QIIME 2 developers will do their best to help you.
