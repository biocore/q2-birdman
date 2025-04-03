# QIIME2 Environment Files for q2-birdman

This directory contains environment files for installing q2-birdman with different QIIME2 distributions.

## Available Environment Files

- `q2-birdman-qiime2-amplicon-2024.5.yml`: For QIIME2 amplicon distribution 2024.5
- `q2-birdman-qiime2-amplicon-2025.4.yml`: For QIIME2 amplicon distribution 2025.4
- `q2-birdman-qiime2-moshpit-2025.4.yml`: For QIIME2 moshpit distribution 2025.4
- `q2-birdman-qiime2-tiny-2024.5.yml`: For QIIME2 tiny distribution 2024.5

## Installation Instructions

To install q2-birdman with a specific QIIME2 distribution, use the following command:

```bash
conda env create \
 -n q2-birdman \
 -f https://raw.githubusercontent.com/biocore/q2-birdman/main/environment-files/q2-birdman-qiime2-<target-distribution>-<target-epoch>.yml
```

Replace `<target-distribution>` with either `amplicon`, `moshpit`, or `tiny` and `<target-epoch>` with either `2024.5` or `2025.4`.

For example, to install q2-birdman with QIIME2 amplicon distribution 2025.4:

```bash
conda env create \
 -n q2-birdman \
 -f https://raw.githubusercontent.com/biocore/q2-birdman/main/environment-files/q2-birdman-qiime2-amplicon-2025.4.yml
```

After installation, activate the environment:

```bash
conda activate q2-birdman
```

## Development Installation

For development purposes, you can install q2-birdman directly from the repository:

```bash
conda env create -n q2-birdman-dev -f environments/q2-birdman-qiime2-tiny-2024.5.yml
conda activate q2-birdman-dev
pip install -e .
``` 