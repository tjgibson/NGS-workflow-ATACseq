# Introduction

This repository contains a Snakemake workflow for performing a basic ATAC-seq analysis pipeline.

This workflow performs the following steps:

1.  Optional: Retrieval of publically available sequencing data from NCBI GEO/SRA.
2.  Optional: merge reads from the same sample that were sequenced on separate lanes or sequencing runs.
3.  Trim adapters from raw reads. I do this using [NGMerge](https://doi.org/10.1186/s12859-018-2579-2) which overlaps the two paired reads and removes overhanging sequence. This conveniently eliminates the need to provide specific adapter sequences.
4.  Align raw reads to a reference genome. The workflow will automatically retrieve any publically available reference genome and build the bowtie2 index.
5.  Filter aligned reads based on alignment quality and discard multi-mapping reads. The filtering parameters can be customized.
6.  Optional: filter aligned reads to discard reads aligning to contigs or mitochondrial genome. This filtering can also be customized.
7.  Separate aligned reads based on fragment size. Small fragments (\< 100bp) are considered to be derived from accessible regions and used for downstream peak calling.
8.  Optional: Perform peak calling using MACS2 on fragments \< 100bp.
9.  Generate bigWig files containing z-score normalized read depth to allow data visualization in the genome browser.
10. Optional: Use DESeq2 to perform differential accessibility analysis.

# Running the workflow

Follow the steps below to run this workflow:

## Step 1: Software installation

Ensure that you have a conda-based Python3 distribution installed.
I recommend [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge).

Snakemake and Snakedeploy are used to deploy and run this workflow.
These can be installed using Mamba:

```         
mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy
```

OR Conda:

```         
conda create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy
```

## Step 2: Deploy the workflow

Activate the conda environment:

```
conda activate snakemake
```

Create a new directory for your analysis and enter it:

```         
mkdir -p path/to/project-workdir
cd path/to/project-workdir
```

Deploy the workflow:

```         
snakedeploy deploy-workflow  https://github.com/tjgibson/NGS-workflow-ATACseq . --branch main 
```

This command will create all the files necessary for running this workflow.

## Step 3: Entering your sample information

*Note:* This workflow requires paired-end data.
It is not recommended to perform ATAC-seq with single-end sequencing.
The inability to separate fragments based on size will reduce sensitivity and may cause artifacts.

### The unit table

The previous step should have created the file `config/units.tsv` in your local project directory.
This is a template that can be modified to contain the information about the samples you wish to analyze.
Before modifying it, take a look at the file to get a sense of how it is formatted.
The fields of this file are described below:

**sample_name:** The name you want to use for each individual sample.
I like to use some sort of naming convention like developmental-stage_genotype_treatment_rep1, where each piece of important information about the sample is separated by an underscore.
If a single piece of information, such as developmental stage, contains multiple words, separate them with a dash: L3-larva_rep1.

**unit_name:** If you sequenced the same sample across multiple lanes or sequencing runs, these will be listed as separate units.
For example, control_rep1_A and control_rep1_B.
If you only have a single fastq file for a sample, then the `unit_name` should be identical to the `sample_name`.

**fq1:** The file path to the location on your computer where the raw fastq file is stored.
These can be located anywhere on your computer, but I suggest placing them in a directory called `data`, `raw`, or `raw_data` within your analysis directory.
Fastq files should be compressed via gzip to save space.

**fq2:** The file path to the second fastq file.

**sra:** This workflow can perform automatic retrieval of publically available data from the Short Read Archive (SRA).
Enter the SRA accession number (e.g. SRR0000001) here.
If you provide both an SRA number and a local fastq file, the workflow will use the local file.

**read_format:** This should be set to `PE` for paired-end.
See note above regarding single-end data.

**sample_group:** Indicates the name to use for each group of replicates.
This name will be used to name bigWig and peak files for merged replicates.
For example, for samples control_rep1 and control_rep2, `sample_group` could be set to control.

### The sample table

Step 2 should also have created the file `config/samples.tsv` in your local project directory.
This is a template that can be modified to include information about the your experimental design and the different conditions tested.
This information will be used to determine the comparisons made by DESeq2 for differential accessibility analysis.
Take a look at the file to get a sense of how it is formatted.
A description of the fields is listed below:

**sample_name:** The name of each individual sample.
This must match the sample names used in the unit table above.

**condition:** The experimental condition each sample is assigned to, for example wild-type vs. mutant.
This will determine which groups are compared by DESeq2 for differential accessibility analysis.

**experiment:** A name for the experiment, for example embryo_ATACseq.
This is useful in the case that multiple experiments are analyzed in the same workflow.
For samples from different experiments, DESeq2 analysis will be performed separately.

## Step 4: Configuring the workflow

The workflow can be customized to suit your needs by editing the file `config/config.yaml`.
Below is a description of the different options that can be customized within the config file.

### Merging of technical replicates

If you have samples that were sequenced on multiple runs or sequencing lanes, change the following lines of the config.yaml file:

```         
mergeReads:
  activate: True
```

This will result in merging of the fastq files prior to alignment.

### Filtering of reads based on alignment to specific chromosomes

When working with Drosophila data (or other organisms with really good genome assemblies), I typically discard reads aligning to unplaced contigs or scaffolds.
For Drosophila, this entails only retaining reads aligning to the major chromosome arms (chr2L, chr2R, Chr3L, chr3R, chrX, and chrY).
This can be done with the following lines of the config file:

```         
filter_chroms: True
keep_chroms: config/keep_chroms.txt
```

You can modify the keep_chroms.txt file to reflect a different filtering strategy.
Make sure that the chromosome naming convention used in the file matches that of your chosen reference genome (see below).
For example, chr2L for UCSC vs. 2L for Ensembl.

### Setting the reference genome

This is where you set the reference genome to use for read alignment.
Give the genome a name and provide a link to the fasta file for the genome.
Reference genomes can be found through the [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/downloads.html), [Ensembl](https://ensemblgenomes.org/), or organism-specific genome databases (e.g. [Flybase](https://ftp.flybase.net/genomes/)).
The pipeline will automatically retrieve the reference genome from the provided link and build the bowtie2 index required for alignment.
Modify the following lines as necessary for your desired reference genome:

```         
ref_genome:
  name: dm6
  link: https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
```

If you want to use a custom genome that is not publically available or a bowtie2 index that you have already built, create a folder called `resources/` in your working directory and copy your fasta file or the bowtie2 index files to this directory.
If using a custom fasta file, rename the file to be `ref_genome.fasta.gz`.
If using a custom bowtie2 index, use the prefix `genome` for the index files (e.g. `genome.1.bt2`).

### DESeq2 options

This section of the config file allows for control over how differential accessibility analysis is performed by DESeq2.
If you want to perform differential accessibility analysis, set the following:

```         
run_diff_accessibility: TRUE
```

If performing DESeq2 analysis, the following lines provide information on how DESeq2 analysis will be done:

```         
diff_accessibility:
  experiments: 
    Experiment_1:
      contrasts:
        wt-vs-mutant:
          - wt
          - mutant
    Experiment_2:
      contrasts:
        control-vs-drugA:
          - control
          - drugA
        control-vs-drugB:
          - control
          - drugB
  model: ~condition
  count_threshold: 50
  padj_cutoff: 0.05
  log2FC_cutoff: 1
```

It is important to note that the information in this section needs to match the information in the sample table.
Each experiment listed here must match and experiment listed in the `experiment` column of the sample table.
The contrasts listed (e.g. control, drugA, drugB) must match the values in the `condition` table of the sample table.

### modifying parameters for specific steps

The final section of the config file allows customization of the parameters for individual steps in the workflow:

```         
params:
  bowtie2_index: ""
  bowtie2_align: "-k 2 --very-sensitive --no-mixed --no-discordant -X 5000"
  filter_multireads: "-f bam -F '(mapping_quality >= 30) and ([XS] == null)'" 
  bigwigs_ind: "--binSize 10"
  bigwigs_merged: "--binSize 10"
  macs2_call_atac_peaks_ind: "-f BAMPE --keep-dup all -g dm --call-summits"
  macs2_call_atac_peaks_merged: "-f BAMPE --keep-dup all -g dm --call-summits"
```

These can be modified as desired to alter these steps.
For example, by default the workflow will run bowtie2 with parameters `-k 2 --very-sensitive --no-mixed --no-discordant -X 5000`.
The parameters listed here in the config file will be passed directly to the corresponding step when it is run.

*Note:* For ATAC-seq analysis, I recommend using the above parameters for peak calling to maximize sensitivity.

## Step 5: Running the workflow

You are now ready to run the workflow.
Before running the workflow, it is recommended to do a "dry run", in which Snakemake will check that all the required input and output files exist, and provide a summary of what the workflow will do.
To do a dry run, invoke the following command:

```         
snakemake -n
```

Examining the output of this dry run can help make sure that samples names and locations of raw files were entered correctly.

If everything looks good, run the workflow with the following command:

```         
snakemake --use-conda -c 1
```

The `-c 1` parameter tells snakemake to use a single core for running the workflow.
Depending on the number of processors/cores available on your computer, this number can be increased, which will speed up execution by allowing multiple steps to run in parallel and allowing steps such as read alignment to use multiple cores.

This workflow can also be run on a computing cluster or using cloud computing services.
See the relevant sections of the Snakemake documentation for more information: [Cluster execution](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) [Cloud execution](https://snakemake.readthedocs.io/en/stable/executing/cloud.html)

# Navigating the output

Once finished, the workflow will have produced various output files that will be located in a folder named `results/`.
The results will be organized into the following subdirectories:

**aligned_reads:** BAM files containing aligned, sorted, filtered reads and BAI bam indices, separated by fragment size.
sample_small.bam will contain fragments \< 100bp, sample_large.bam will contain fragments \> 150bp, and sample_total.bam will contain all fragments.
See the [original ATAC-seq paper](https://doi.org/10.1038/nmeth.2688) for more details on fragment size filtering.
A number of different intermediate BAM files will be generated as the workflow runs.
Upon completion, intermediate files will be deleted and only only the filtered, sorted BAM files will be retained.

**bigwigs:** z-score normalized bigWig files for individual replicates and merged replicates.
Separate bigWigs are generated for small, large, and total fragments.

**peaks:** narrowPeak and bed files containing peaks called by MACS2.
This includes peaks called on individual and merged replicates.
Peak calling is performed on small fragments.

**count_tables:** The count tables produced by featureCounts.
the rows correspond to individual ATAC-seq peaks, and the columns correspond to samples.
These tables will be used by DESeq for performing differential accessibility analysis.

**DEseq2:** The results of differential accessibility analysis using DESeq2.
This will contain results tables summarizing differential accessibility.
There will be one table per experiment and the file names will end in results.tsv.
This folder will also contain a .dds file for each experiment.
This file contains the DESeq R object and can be imported into R to extract additional data and information for each experiment.

# Combining workflows for multiple data types

It may be desirable to process multiple data types (e.g. ChIP-seq, RNA-seq, ATAC-seq) and then integrate those results into downstream analysis.
This can be done in a single Snakemake workflow that contains multiple "modules" for the different data types.
For an example of how to build a complex workflow including multiple datasets as well as custom downstream analysis, see [this repository](https://github.com/tjgibson/S2_pioneers_manuscript) for an example.

# To-do

In the future, I hope to add some of the following features:

-   Add small, test fastq files to repo to allow for automated testing of rules

-   Create an HTML summary file containing read alignment and filtering statistics, number of peaks called, and other useful information about the workflow.

-   The optional step to take unaligned reads, BLAST them, and provide a summary of which organisms the unaligned reads likely originated from.

-   Generate fragment size distributions for each sample

-   Generate PCA and volcano plots
