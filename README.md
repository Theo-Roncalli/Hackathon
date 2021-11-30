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

Docker should be installed as well:
```bash
sudo apt install docker
```

Once nexflow is installed, it will automatically retrieve the docker images used within the pipeline.

## Workflow DAG

Nextflow workflows should form a _DAG_ (i.e. directed acyclic graph), which represents the flow of data through the different steps 
required to produce the final result. 

This pipeline will generate a set of figures, representing differential gene expression analysis of RNA-Seq data.

![dag](https://user-images.githubusercontent.com/28574085/143959065-28154038-71db-487c-a816-6777934db3d3.png)

## Hardware requirements

A machine with at least 32 GB of **FREE** RAM (to create the index on the reference genome).
Recommended configuration is 64 GB, by default the index creation process is configured to use 50 GB.

## Executing The Workflow

1. Clone the repo to your machine
```bash
git clone https://github.com/bio-TAGI/Hackathon.git
cd Hackathon
```
2. Create and activate the virtual environment
```bash
conda env create -f nextflow_conda_env.yml
conda activate nextflow
```
3. Run the wokflow with default parameters.
```bash
cd Nextflow
nextflow run main.nf
```
4. If you had to stop the workflow run, or if some error occurred, you can always resume the execution as follows:
```bash
nextflow run main.nf -resume
```
5. Specifying parameters from the command line
```bash
nextflow run main.nf --param1 value1\
--param2 value2\
--paramn valuen # these are generic names, not actual parameters for the pipeline
```

### Optional parameters

* `index_cpus` (number of cpus reserved for the genome indexation process,   `default=14`)
* `mapping_cpus` (idem. for the mapping process, used to create BAM files, `default=14`)
* `counting_cpus` (idem. for the counting process. `default=7`)
* `mapping_memory` (RAM reserved for mapping . `default=50GB`)


If you already possess some of the files needed to execute the pipeline, you can specify them as follows:

* `reads` (path pointing to a directory containing the `fasterq` files)
* `genome` (path pointing to a directory containing the genome FASTA file)
* `index` (Répertoire contenant les fichiers d’index)
* `mapping` (Répertoire contenant les fichiers BAM)
* `counting` (Chemin d’accès entier au fichier de comptage – comprend le fichier lui-même)
* `metadata` (Chemin d’accès entier au fichier de métadonnées – comprend le fichier lui-même)

If unspecified, they will these files will be downloaded following the default values from the config file : [nextflow.config](https://github.com/bio-TAGI/Hackathon/blob/main/Nextflow/nextflow.config). These too, can be tweaked and overriden:

* `ids` List of SRR accession number to fetch paired-end fastq files. 
  * default=`['SRR628582', 'SRR628583', 'SRR628584', 'SRR628585', 'SRR628586', 'SRR628587', 'SRR628588', 'SRR628589']`
* `genome_url` URL to download the reference genome. 
  *  default `ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`
* `annotation_url` URL to donwload the reference genome's annotation. 
  * default `ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz`

## Caveats

* A good internet connection is required. Retrieving `fastq` can be really slow and is thus a bottleneck.
* [`fasterq-dump`](https://github.com/ncbi/sra-tools/wiki) will randomly segfault. 
At first we thought this was caused by connection problems, but running `ping` ruled this out. Apparently, the segfault is [a known issue](https://github.com/ncbi/sra-tools/issues/518).
* The workflow will inevitably fail if you try building the genome's index on a machine with less than ~30 GB of RAM available.
  * As a general rule, tweak all parameters to reasonable values that fit your setup and needs. We don't know your hardware, you do ;)
