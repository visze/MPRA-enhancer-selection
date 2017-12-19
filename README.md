
[![License](https://img.shields.io/badge/license-GNU%20GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.txt)

# Generate enhancer probes for MPRA experiments

This folder contains a pipeline to select positive and negative enhancer regions using DNase seq from the epigenetic roadmap.

An "enhancer" is defined as an open DNase accessibility which is active in several samples and is not located in front of a transcription start site.

Some stats about the samples in the peigenetic rodmap:

- Number of samples: 127
- Number of Tissue Groups: 19
- Largest Tissue Groups:
  - ENCODE2012: 16 (samples from encode. This is a mixture of tissues)
  - Blood & T-cell: 14
  - Digestive: 12
  - Other: 11  
  - Brain: 10


## Pipeline description

The first part of the probe generation is split up in "positive" and "negative" groups. But in general there are steps:

1. Find common narrow peaks and filter them
2. Remove regions that are close to a  transcription start site (promoter removal)
3. Get the center of the peaks by maximum signals
4. Create probes using the new center and the probe length

The first part is slightly different between the positive enhancers (active in lot's of samples) and negative enhancers (active only in a tissue group)

### 1.a Positive enhancers

Here we will select regions that are always active over all samples. On a MPRA experiment this region should be "activate".

1. Running `bedtools intersect` with parameters `-c` `-r` and `-f`. One random sample will be selected as `-a` input. All samples (including `-a`) used as `-b` input. `bedtools` will count the overlaps of the all samples(`b`) with the regions in `a`. `-f` is the minimum fraction of the overlap (e.g. 1.0, 0.9,..). because of `-r` the overlaps must hold for both direction. Region in `a` in regions in `b` and vice versa.
2. Filtering the data using a threshold for the minimum number of samples that have the regions. In total there are 127 samples. So with a threshold of 110 the DNase narrow peak is present in 85% of all samples.

### 1.b Negative enhancers

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
# clone from github
git clone https://github.com/visze/MPRA-enhancer-selection.git
# enter the project
cd MPRA-enhancer-selection
```

### Create a conda environment

Miniconda3 must be installed before. Then use a terminal and run the following commands:

```
# Install dependencies and programs into isolated environment
conda env create --file environment.yml

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
