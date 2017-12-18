[![license](https://img.shields.io/badge/license-GNU%20GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.txt)

# Generate enhancer probes for MPRA experiments

This folder contains a pipeline to select positive and negative enhancer regions using DNase seq from the epigenetic roadmap.

An "enhancer" is defined as an open DNase accessibility which is active in several samples and is not located in front of a transcription start site.


## Pipeline description

The following data is generated:

### Positive enhancers

Here we will select regions that are always active over all samples. On a MPRA experiment this region should be "activate".

### Negative enhancers

Here we will select regions of open DNase accessibility that are only active in one tissue (e.g. blood) but not in any  other tissue. On a MPRA experiment with a different tissue this region should be deactivated

## Requirements

`Git` and `conda` (Miniconda3) are needed to run the pipeline. Nothing else. `conda` will manage all programs that will be needed

## Installation from Github sources

1. Clone this repository.
2. Create a conda environment.
3. Edit the config file `config.yml`.
4. Run snakemake pipeline.

### Clone this repository

You need git installed for this! Use your terminal and move to a folder where you want to checkout the pipeline. Then run:

```
git clone https://github.com/visze/mpra_enhancer_selection.git
```

### Create a conda environment

Miniconda3 must be installed before. Then use a terminal and run the following commands:

```
# Install dependencies and programs into isolated environment
conda env create --file environment.yaml

# activate environment
source activate dnase
```

## Edit the config file `config.yml`

All settings are configured via the configuration file `config.yml`. All variables are described in there. E.g. it is important to set the path to a reference fasta file of the genome (Build hg19 or GRCh37 AND indexed) via the parameter `reference_genome`.

## run the pipeline
```
# Execute the workflow (dryrun to see what it does)
snakemake -n

# Really start it (takes around one day). You can speedup by using multiple cores: --cores 4
snakemake

# to deactivate your workflow:
source deactivate
```
