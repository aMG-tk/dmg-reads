
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

### Using mamba

```bash
mamba install -c conda-forge -c bioconda -c genomewalker dmg-reads
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

dReads will take a TSV file produced from [metaDMG](https://metadmg-dev.github.io/metaDMG-core/) and extract the reads from a BAM file. One can select a list of taxa and ranks to extract the reads from.

For a complete list of options:

```bash
$ dReads --help
usage: dReads [-h] -m METADMG_RESULTS -f METADMG_FILTER [-b BAM] [-p PREFIX] [--combine]
              [-T TAXONOMY_FILE] [-r RANK] [-M SORT_MEMORY] [-t THREADS] [--debug] [--version]

A simple tool to extract damaged reads from BAM files

optional arguments:
  -h, --help            show this help message and exit
  -m METADMG_RESULTS, --metaDMG-results METADMG_RESULTS
                        A file from metaDMG ran in local mode (default: None)
  -f METADMG_FILTER, --metaDMG-filter METADMG_FILTER
                        Which filter to use for metaDMG results (default: None)
  -b BAM, --bam BAM     The BAM file used to generate the metaDMG results (default: None)
  -p PREFIX, --prefix PREFIX
                        Prefix used for the output files (default: None)
  --combine             If set, the reads damaged and non-damaged will be combined in one fastq file
                        (default: False)
  -T TAXONOMY_FILE, --taxonomy-file TAXONOMY_FILE
                        A file containing the taxonomy of the BAM references in the format
                        d__;p__;c__;o__;f__;g__;s__. (default: None)
  -r RANK, --rank RANK  Which taxonomic group and rank we want to get the reads extracted. (default:
                        None)
  -M SORT_MEMORY, --sort-memory SORT_MEMORY
                        Set maximum memory per thread for sorting; suffix K/M/G recognized (default:
                        1G)
  -t THREADS, --threads THREADS
                        Number of threads (default: 1)
  --debug               Print debug messages (default: False)
  --version             Print program version
```

One would run `dReads` as:

```bash
dReads -m RISE505_MA873_L1.tp-mdmg.local.weight-1.csv.gz -b RISE505_MA873_L1.dedup.filtered.sorted.bam -f '{ "Bayesian_D_max": 0.1, "Bayesian_z": 2 }' --prefix RISE505_MA873_L1 --taxonomy-file gtdb-r202-organelles-viruses.tax.tsv --rank '{"genus": "Yersinia", "class":"Bacilli"}
```

The filter is a JSON object where the keys are one of the metaDMG results headers. If `--taxonomy-file` and `--rank` are set, the reads will be extracted from the selected taxonomic group and rank. If `--only-damaged` is set, only the damaged reads will be extracted. If `--combine` is set, the damaged, non-damaged and multi-mapped reads will be combined in one fastq file.

The previous command will produce the following files:

```bash
├── RISE505_MA873_L1.c__Bacilli.non-damaged.fastq.gz
└── RISE505_MA873_L1.g__Yersinia.damaged.fastq.gz
```

If the `--combine` and the `--only-damaged` flag are not set, `dReads` will produce three files per taxa/rank or BAM file:

- `*.damaged.fastq.gz`: The reads mapped to a reference that shows damage
- `*.non-damaged.fastq.gz`: The reads mapped to a reference that does not show damage
- `*.multi.fastq.gz`: The reads mapped to multiple references which are damaged and non-damaged


# Using taxonomies
To be able to extract reads from specific taxa and/or ranks, one needs to provide a taxonomy file. This file should be a TSV file with the following format:

```
ACCESSION\td__Bacteria;l__Bacteria;k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Yersinia;s__Yersinia pestis
```

`ACCESSION` is the reference accession in the BAM file. The taxonomy is separated by `;` and the taxonomic groups are separated by `__`. The taxonomic groups recognized by `dReads` in `--taxonomy-file` and `--rank` are:
  - **domain**: `d__`
  - **lineage**: `l__`
  - **kingdom**: `k__`
  - **phylum**: `p__`
  - **class**: `c__`
  - **order**: `o__`
  - **family**: `f__`
  - **genus**: `g__`
  - **species**: `s__`

> **Note**: The taxonomic groups are case sensitive and one can include as many as desired. For example, if one wants to extract the reads from the genus *Yersinia* and the class *Bacilli*, one would use `--rank '{"genus": "Yersinia", "class":"Bacilli"}`.


