
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

### 1.A Positive enhancers

Here we will select regions that are always active over all samples. On any MPRA experiment with any tiisue or cell-line this region should be "activate".

1. Running `bedtools intersect` with parameters `-c` `-r` and `-f`. One random sample will be selected as `-a` input. All samples (including `-a`) used as `-b` input. `bedtools` will count the overlaps of the all samples(`b`) with the regions in `a`. `-f` is the minimum fraction of the overlap (e.g. 1.0, 0.9,..). because of `-r` the overlaps must hold for both direction. Region in `a` in regions in `b` and vice versa.
2. Filtering the data using a threshold for the minimum number of samples that have the regions. In total there are 127 samples. So with a threshold of 110 the DNase narrow peak is present in 86% of all samples.

The next steps are in common with the negative enhancers.

### 1.B Negative enhancers

Here we will select regions of open DNase accessibility that are only active in one tissue (e.g. blood) but not in any  other tissue. On a MPRA experiment with a different tissue this region should be deactivated.

1. First we select a tissue group line `Brain`.
2. Then we run the "counting" on that group. This is again done with `bedtools intersect` (see 1.A Positive enhancers). One random sample of the group will be selected as `-a` input. All samples of the same group (including `-a`) used as `-b` input.
3. Now we do the same analysis but we will change the `-b` samples to all samples that are not in the selected group (e.g. are not brain).
4. Filtering the data.
    1. First we have a threshold for the minimum number of active regions within the group. E.g. if we select as group `Brain` and set the threshold to `8`the 80% of all Brain samples should have this region actives.
    2. Because we want to have that region exclusive active in the tissue group we filter on teh second analysis with a maximum threshold. E.g. a threshold of 10 for the Brain group means that teh active region Brain region is only active in 10 or less other samples. This is only around 8.5% (10 brain samples, 117 other samples).

### Select Enhancers

Here DNase data is used. In theory this can be Promoter and Enhancer regions. because it is known that promoter are in front of a TSS we can simple exclude regions in front of genes. E.g. Here I use ensembl or Gencode transcripts and remove all regions that are in the range of 2kb in from of sthe TSS. The direction of the gene is considered here!

### Center of Peaks

We are still left with regions of one sample and of different sizes. The site can range from 200 to several thousand bp. But we need only a small defined length (e.g. 171).

So we have to define a new center. Therefore I use the imputed signals that overlaps the region and for every sample I selected the maximum as center. Then (for simplicity) the average position is chosen.

### Creating Probes

Here we need a length of the probe (e.g. 171) which will be placed around the new computed center. We will get two output files. One a multi-fatsta file with the region and the DNA of that regions. Second a BED file with just the region of the probe in the genome. Final files are in the folder `results/RoadMap/positives|negatives/design/final/`.

### Naming Convention of Files:

1. Positive enhancers: `minSamples_120.overlap_0.95.gencode_27_minDistanceToTSS_2.0.probeLength.171.bed` means:
    1. `minSamples_120`: 120 of 127 samples have this region called
    2. `overlap_0.95`:`bedtools intersect` requires an overlap of 95%
    3. `gencode_27_minDistanceToTSS_2.0`: Regions lies not within 2.0kb of a TSS and the transcription db gencode in version 27 is used.
    4. `probeLength.171`: the length of the region and the extracted sequence is 171bp
2. Negative enhancers: `Brain.minsamplesWithinGroup_8.maxOther_5.groupOverlap_0.95.ensembl_75_minDistanceToTSS_2.0.probeLength.171.bed` means:
    1. `Brain`: Tissue group is brain
    2. `minsamplesWithinGroup_8`: Minimumm of 8 from all 10 Brain samples should be active (80%)
    3. `maxOther_5`: Maximum of 5 other active samples are allowed. So around 4% (5 of 117).
    4. `overlap_0.95`:`bedtools intersect` requires an overlap of 95% (within and to other groups)
    5. `ensembl_75_minDistanceToTSS_2.0`: Regions lies not within 2.0kb of a TSS and the transcription db Ensembl in version 75 is used.
    6. `probeLength.171`: the length of the region and the extracted sequence is 171bp

### Files:

#### Bed files (0 based, start inclusive, end exclusive)

1. Positive Enhancers columns: Contig, start, stop, number of active samples.
```
chr1	121484654	121484825	114
chr1	153959295	153959466	119
chr1	156186345	156186516	121
chr1	161582272	161582443	110
chr1	223254706	223254877	119
chr2	87623843	87624014	114
chr2	92306027	92306198	120
chr2	92319060	92319231	123
chr2	92321512	92321683	118
chr2	133016842	133017013	120
chr3	125634806	125634977	115
  ```
2. Negative Enhancers columns: Contig, start, stop, number of active samples within the group, Number of active samples outsite of the group.
```
chr1	32158975	32159146	8	3
chr1	76747146	76747317	8	4
chr1	76856259	76856430	8	1
chr1	111299189	111299360	8	5
chr1	117598077	117598248	8	4
chr1	190547264	190547435	8	4
chr1	203128185	203128356	8	2
chr1	221017113	221017284	8	5
chr2	50201277	50201448	8	3
chr2	110457372	110457543	8	2
chr2	188392116	188392287	8	5
chr2	213159651	213159822	8	3
chr3	32352787	32352958	8	4
chr3	52060168	52060339	8	0
chr3	120206219	120206390	8	3
chr3	123241312	123241483	8	5
  ```

#### Fasta files:
1. Positive Enhancers
    1. Header:
      - contig:start-stop (1 based, both inclusive)
      - active_count:114 (number of active samples)
    2. Sequence
```
>chr2:87623844-87624014 active_count:114
tagccagccaggcccgccagccagccagccagcgagccaagccagccaagccagccagcctgccaagccagccggccagccaagctagccaatccactcagccactcaagccagccaagtcacccggccatccaagccagccaagccagtcagccagcccagacagccaag
>chr2:92306028-92306198 active_count:120
ctctttttgtggaatctgcaagtgcatatttagctagatttgacgatttcgttggaaacgggattacatataaaaagcagacagcagcattctcagaaactcctttgtgatgtttgcattcaagtcacagagttgaacattccctttcatagagcaggattgaaaaactct
```
1. Negative Enhancers
    1. Header:
      - contig:start-stop (1 based, both inclusive)
      - active_count:9_2 (number of active samples within group; number of active samples outsite group)
    2. Sequence
```
>chr1:4768993-4769163 active_count:9_2
ACTTGAAGAGGAAAAACAAATCGACCTCTCCCTGCCACTGTTGCAATTGGTTGGTTTTTCTGCATAACAGCTGGGTGTCTTAGAAATGAGGGGGTTTCTATAGTAACCAATTACAGCCATGATTGGTGAAAAATCACAGAAATATCCTGTGTGTGAAGTTATGCCAGCGAG
>chr1:32158976-32159146 active_count:8_3
ATTTGGCAAATCTACAGTCTCAAGCAGGGGGGTATCCCATCCCCCACACCCGTACCCAAGCACATCAGCTCACACACAGCACAGCCGGAGGCATATGGACACACACATgcacggtggcagaaacccatgctggctgagtgtcaggctgcctgacttcaaatcctgtacttg
```

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
