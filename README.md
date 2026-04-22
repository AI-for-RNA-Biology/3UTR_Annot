# 3'UTR Annotation

## Summary

The pipeline reannotates 3'UTR of expressed transcripts based on read coverage. In brief, it is divided into 2 main parts:
1. The longest expressed 3' UTR in each sample is extracted, starting from the provided annotation. 

Continuously transcribed regions are identified using a sliding window across the genome requiring a minimum coverage of 7 reads in more than 80 positions per window of 100 bp. Neighbouring regions separated by lowmappable regions are merged as described in (Miura et al., 2013). Expressed fragments are associated with matching strand overlapping 3’UTR using Ensembl Rn5 version 78 (v78) (Flicek et al., 2013). Isolated expressed fragments that do not overlap with any feature are associated with the closest 3’UTR if (1) the 3’UTR is <10kb and (2) there are no intervening annotations. The assigned expressed fragments are filtered to exclude potential intragenic transcription, overlapping transcripts, and retained introns as described in (Miura et al., 2013). If the expressed fragment continued beyond the end of the annotated 3’UTR, the fragment is taken as the new 3’ end.

2. Expressed alternative 3' UTR isoforms are identified based on changes in read coverage.

A marked change in the level of coverage in the 3’ to 5’ end direction is expected to occur at the boundaries of alternative polyadenylation sites due to the accumulation of the reads in 3’ termini of the transcripts. To identify alternative 3’UTR isoforms the Segmentor3IsBack R package (Cleynen et al., 2014) is used. Segmentor3IsBack models nucleotide level read coverage along 3’UTRs with a negative binomial distribution and identifies segments of spatially coherent coverage along individual 3’UTRs. Prior to coverage segmentation, nucleotide level read coverage is smoothed using a running median of width 150 nt. As previously reported, each 3’ end can have multiple cleavage positions in a small window (Tian et al., 2005); thus change-points located within 25 nt from one another are clustered together and the most promoter-distal change-point is selected. The algorithm is applied to every sample, both using raw coverage and log2-scaled coverage. All annotations are then merged and 3’ end located within 50 nt distance are clustered by selecting the most promoter-distal annotation. 

Ensembl annotated and newly annotated 3’UTR isoforms with PAS motifs in the -150 nt to +100 nt region surrounding the 3’ end termini are flagged. 4 flags for the presence of PAS motifs are included: 1) Presence of the 12 canonical and non-canonical PAS motifs (AATACA, ATTAAA, TATAAA, AATATA, AATAGA, AGTAAA, AATGAA, ACTAAA, CATAAA, GATAAA, AAGAAA, and AATAAA) 2) PolyA_DB (Lee et al., 2007) 3) PolyASite (Herrmann et al, 2020), and 4) PolyaID (Stroup et al, 2023). 

## Dependencies
### Conda and Snakemake 
Conda is required to run Snakemake. 
For Linux systems, this installer can be used: [Miniconda3-latest-Linux-x86_64.sh](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

Download the installer and run:

```
bash Miniconda3-latest-Linux-x86_64.sh
```

After installing conda, install mamba which is required for snakemake

```
conda install -n base -c conda-forge mamba
```

Then install snakemake and singularity

```
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

### Other Software Dependencies
All other dependencies are listed in yml files found in [./envs/](../envs/). These are automatically installed when the pipeline is invoked using Snakemake

Alternatively, conda enviroments could be installed using the yaml files before running the pipeline. The conda enviroments within each rule should then be modified accordingly. 
## Input

Path to input files are specified in the [config.yml](../config/config.yml) file.

### BAM Files (bam_dir)
Path to directory that contains bam files. The defauly suffix is ".Aligned.sortedByCoord.out.bam", which is produced by STAR. This should be modified if needed. 

### Annotation Files
The required enseml annotation files for human (GRCh38) are stored on Zenodo ([https://doi.org/10.5281/zenodo.19695412](https://doi.org/10.5281/zenodo.19695412)) and can be placed in the [./workflow/data/)](./workflow/data/) directory. 
For other species, please create the relevant annotation files. 

### Output Directory
Path to output directory. 

## Running the Pipeline

By default, the entire pipeline is run and the final annotation file will be stored in the specified output directory. 
To run only certain rules, change the "Input" rule in the [Snakefile](./workflow/Snakefile)

First, set the correct directory and activate the environment:

```
cd ./workflow
conda activate snakemake
```

To check the rules that'll be invoked based on the output:

```
snakemake --dry-run
```
To run the pipeline: 
```
snakemake --use-conda --rerun-trigger mtime
```

### HPC(SLURM)
The pipeline can be invoked to run in HPC enviroments (SLURM) using the following scripts
Option 1: 
```
sbatch run_Snakemake.sh
```
The number of concurrent jobs (ie samples processed) should be specified under --jobs within the script, along with the number of CPU cores. 

Option 2: 
```
./run_Snakemake_SLURM.sh
```
The number of concurrent jobs (ie samples processed) should be specified under --jobs within the script. 
Other specifications for SLURM (e.g. number of cores, memory, runtime) can be specified by modifying [../config/config.yaml](../config/config.yaml)

## Output
[hg38_Denovo_Ensembl_RefSeq_Filt_Reformat.gtf](hg38_Denovo_Ensembl_RefSeq_Filt_Reformat.gtf): Input for quantification pipeline

[hg38_Denovo_Ensembl_RefSeq_Filt.gtf](hg38_Denovo_Ensembl_RefSeq_Filt.gtf): Ensembl, RefSeq, and De-novo annotations that are filtered for overlap. Flags include sample of origin (if de-novo) and whether a motif is present or if it is annotated in an existing database

[hg38_Denovo.gtf](hg38_Denovo.gtf): All De-novo annotations


