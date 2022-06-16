
# dReads: a tool to extract damaged reads from BAM files


[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/dmg-reads?include_prereleases&label=version)](https://github.com/genomewalker/dmg-reads/releases) [![dmg-reads](https://github.com/genomewalker/dmg-reads/workflows/dReads_ci/badge.svg)](https://github.com/genomewalker/dmg-reads/actions) [![PyPI](https://img.shields.io/pypi/v/dmg-reads)](https://pypi.org/project/dmg-reads/) [![Conda](https://img.shields.io/conda/v/genomewalker/dmg-reads)](https://anaconda.org/genomewalker/dmg-reads)

A simple tool to extract damaged reads from BAM files

# Installation

We recommend having [**conda**](https://docs.conda.io/en/latest/) installed to manage the virtual environments

### Using pip

First, we create a conda virtual environment with:

```bash
wget https://raw.githubusercontent.com/genomewalker/dmg-reads/master/environment.yml
conda env create -f environment.yml
```

Then we proceed to install using pip:

```bash
pip install dmg-reads
```

### Using conda

```bash
conda install -c conda-forge -c bioconda -c genomewalker dmg-reads
```

### Install from source to use the development version

Using pip

```bash
pip install git+ssh://git@github.com/genomewalker/dmg-reads.git
```

By cloning in a dedicated conda environment

```bash
git clone git@github.com:genomewalker/dmg-reads.git
cd dmg-reads
conda env create -f environment.yml
conda activate dmg-reads
pip install -e .
```


# Usage

dReads will take a TSV file produced from [metaDMG](https://metadmg-dev.github.io/metaDMG-core/) and extract the reads from a BAM file. 
For a complete list of options:

```
$ dReads --help
usage: dReads [-h] -r METADMG_RESULTS -f METADMG_FILTER [-b BAM] [-p PREFIX] [--debug] [--version]

A simple tool to extract damaged reads from BAM files

optional arguments:
  -h, --help            show this help message and exit
  -r METADMG_RESULTS, --metaDMG-results METADMG_RESULTS
                        A file from metaDMG ran in local mode (default: None)
  -f METADMG_FILTER, --metaDMG-filter METADMG_FILTER
                        Which filter to use for metaDMG results (default: None)
  -b BAM, --bam BAM     The BAM file used to generate the metaDMG results (default: None)
  -p PREFIX, --prefix PREFIX
                        Prefix used for the output files (default: None)
  --debug               Print debug messages (default: False)
  --version             Print program version
```

One would run `dReads` as:

```bash
dReads -r 8e5a029071.tp-mdmg.local.weight-1.csv.gz -b 8e5a029071.dedup.filtered.bam -f '{ "Bayesian_D_max": 0.1, "Bayesian_z": 2.5 }' -b 8e5a029071.dedup.filtered.bam  --only-damaged
```

The filter is a JSON object where the keys are one of the metaDMG results headers.

This will produce the following files:

```bash
├── 8e5a029071.damaged.fastq.gz 
└── 8e5a029071.nondamaged.fastq.gz
```



