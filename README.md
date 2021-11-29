# Hackaton

Analysis of mutations at codon 625 of SF3B1 gene in uveal melanoma.

## Dependencies

The pipeline runs on [nextflow](https://www.nextflow.io/) a domain-specific language created to automate data-analysis pipelines whilst maximising reproducibility.
Nextflow enables scientists to focus on their analyses, isolating different parts of the pipeline into processes whose dependencies can be dealt with using containers and virtual environments with technologies such as [Docker](https://www.docker.com/), [Singularity](https://singularity.hpcng.org/), and [Anaconda](https://www.anaconda.com/products/individual).

The recommended way to install `nextflow` is via `conda`, using [the environment file](https://github.com/bio-TAGI/Hackathon/blob/main/nextflow_conda_env.yml).
```bash
conda env create -f nextflow_conda_env.yml # will create an env called "nextflow"
conda activate nextflow
# You can edit the file at your choice, specially if the environment name conflicts
# with a preexisting conda env on your system
```

## Workflow DAG

Nextflow workflows should form a _DAG_ (i.e. directed acyclic graph), which represents the flow of data through the different steps 
required to produce the final result. 

This pipeline will generate a set of figures, representing differential gene expression analysis of RNA-Seq data.

![dag](https://user-images.githubusercontent.com/28574085/143959065-28154038-71db-487c-a816-6777934db3d3.png)

## Hardware requirements

The workflow

```
nextflow run nextflow/main.nf
```

## Caveats

* A good internet connection is required. Retrieving `fastq` can be really slow and is thus a bottleneck.
* [`fasterq-dump`](https://github.com/ncbi/sra-tools/wiki) will randomly segfault. 
At first we thought this was caused by connection problems, but running `ping` ruled this out. Apparently, the segfault is [a known issue](https://github.com/ncbi/sra-tools/issues/518).

