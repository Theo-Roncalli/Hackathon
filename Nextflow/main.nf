/* nextflow run main.nf --reads ../Data/Reads --genome ../Data/Genome --index ../Data/Index \
--index_cpus 7 \
--mapping_cpus 7 \
--mapping_memory '12GB'
*/


/* nextflow run main.nf --reads ../Data/ReadsTmp \
--genome_url http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz \
--annotation_url ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz \
--index_cpus 7 \
--mapping_cpus 7 \
--mapping_memory '12GB' \
--index /home/ubuntu/Documents/AMI2B/Hackathon/Projet/Data/IndexTmp \
--genome /home/ubuntu/Documents/AMI2B/Hackathon/Projet/Data/GenomeTmp

*/

/*
    Nexflow pipeline to perform a full RNA-seq analysis (differential expression)
    from a series of SRA accession numbers and a reference genome.

    All parameters are defined in the `params` scope
    within the "nextflow.config" configuration file.
    However the default values can be overriden by 
    specifing them as command line arguments.

    Usage:
    -----
    * Launch the pipeline using the default params
        $ nextflow run main.nf

    * Override params.ids (SRA accession numbers)
        TODO



    ######################################################################
                    WORKFLOW DIAGRAM
    ######################################################################

    ---------------     --------------------
    | SRA entries |     | Reference genome |
    ---------------     --------------------
        ||                     ||
        ||                     ||
        ||              ------------------
        ||              | Index building |
        ||              ------------------
        \\                     ||
         \\                    ||
          \\           ------------------------
           \\==========| Mapping RNA-seq data |
                       | to reference genome  |
                       ------------------------
                               ||
                               ||
                        -------------------------------
                        | Building a count matrix     |
                        | (for genes accross samples) |
                        -------------------------------
                               ||
                               ||
                        -----------------------
                        |Perform differential |
                        |expression analysis  |
                        -----------------------
                               ||
                               ||
                        -----------------
                        | Build reports |
                        -----------------

    ######################################################################
*/

nextflow.enable.dsl=2

log.info """\
D I F F R N A - N F  v0.1.0 
===========================
genome_url       :  $params.genome_url
annotations_url  :  $params.annotation_url
SRA ids          :  $params.ids
readlength-1     :  $params.sjdbOverhang
reads_dir        :  $params.reads
genome_dir       :  $params.genome
index_dir        :  $params.index
"""

process Fasterq {
    /*
    Use ncbi sra-tools' fasterq-dump to rapdily retrieve
    and extract fastq files from SRA-accession-numbers.

    arguments:
    ---------
        ids: a SRA accession number (a channel containing many of them)

    output:
    ------
        A chanel of tuples, containing the id (SRA accession)
        and the path to the fastq files (id, [id_1.fastq, id_2.fastq])

    For further information refer too fasterq-dump's documentation:
    https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
    */

    tag "Downloading ${ids}..."

    input:
        val ids

    output:
        tuple val("${ids}"), path("*_{1,2}.fastq")

    script:
    """
    fasterq-dump ${ids}
    """
}

process Genome {
    /*
    Use `wget` to retrieve a genome, then expand it using `gunzip`.

    arguments:
    ---------
        url: url pointing to the desired reference genome

    output:
    ------
        A path (glob), of all the uncompressed parts of the genome.
    */
    tag "Retrieving genome: ${genome_url}, annotation: ${annotation_url}"

    input:
        val genome_url
        val annotation_url

    output:
        path "*.f*a"
        path "*.gtf"
        // DISCUSSION :
        // gunzip expands the patterns `*.fna.gz` and `*fa.gz`
        // shouldn't this wildcard be the same ?

    script:
    """
    #!/usr/bin/env bash
    echo "Downloading genome..."
    wget ${genome_url}
    [[ ${genome_url} == *.gz ]] && gunzip *.gz || echo "File already unzip."
    echo "Downloading annotations..."
    wget ${annotation_url}
    [[ ${annotation_url} == *.gz ]] && gunzip *.gz || echo "File already unzip."
    """
}

process Index {
    /*
	Create an index for the desired reference genome.

    arguments:
    ---------
        genome_file: a path, pointing to the genomeFastaFiles

    output:
    ------
        path: A directory containing the genome index generated by STAR

    params:
    ------
        params.index_cpus: an integer, specifying the number of threads to be used
                           whilst creating the index.
    */

    cpus params.index_cpus
    tag "Creation of the index"

    input:
        path genome_path
        path annotation_path

    output:
        path "GenomeDir"
    
    script:
    """
    #!/usr/bin/env bash
    STAR --runThreadN ${task.cpus} \
         --runMode genomeGenerate \
         --genomeFastaFiles ${genome_path} \
         --sjdbGTFfile ${annotation_path} \
         --sjdbOverhang ${params.sjdbOverhang}
    """

}

process Mapping {
    /*
	Create the mapping for the RNA-seq data.
    */

    memory params.mapping_memory
    cpus params.mapping_cpus
    tag "Mapping ${fastq_files[0]} to reference genome index."

    input:
        each fastq_files
        path index_path

    output:
        path "${fastq_files[0]}.bam"
    
    script:
    """
    #!/usr/bin/env bash
    echo "Mapping computation for ${fastq_files[0]}..."
    STAR  --runThreadN ${task.cpus} \
    	  --outFilterMultimapNmax 10 \
    	  --genomeDir ${index_path} \
    	  --readFilesIn ${fastq_files[1][0]} ${fastq_files[1][1]} \
    	  --outSAMtype BAM SortedByCoordinate
    mv Aligned.sortedByCoord.out.bam ${fastq_files[0]}.bam
    echo "Done for ${fastq_files[0]}"
    """
}

process Counting {
    /*
	Create the counting matrix for the RNA-seq data.
    */

    cpus params.counting_cpus
    tag "Counting the number of reads per gene."

    input:
    path annotation_file
    path bam_files

    output:
    path "counts.txt"

    script:
    """
    #!/usr/bin/env bash
    echo ${bam_files}
    featureCounts -p -T ${params.counting_cpus} -t gene -g gene_id -s 0 -a ${annotation_file} -o counts.txt ${bam_files}
    """

}


//     echo "here the bam files :${bam_files[0]}..."
//    echo "here the bam files :${bam_files}..."


// samtools index ${fastq_files[0]}.bam

workflow {

    // SEE EXPLANATION OF THE NEW PROGRAM STRUCTURE AFTER THE WORKFLOW

    // Retrieve RNA-seq data (fastq files / SRA accession numbers)
    ids = Channel.fromList(params.ids)
    fastq_files = (
        params.reads == null ?
        Fasterq(ids) :
        Channel.fromFilePairs("${params.reads}/SRR*_{1,2}.fastq*", checkIfExists:true)
    )
//    fastq_files.view()

    // Retrieve genome and annotations
    if (params.genome == null) {
    	Genome(params.genome_url, params.annotation_url)
        path_genome = Genome.out[0]
        path_annotation = Genome.out[1]
    } else {
    	Channel
            .fromPath("${params.genome}/*.fa", checkIfExists: true)
            .set{ path_genome }
    	Channel
            .fromPath("${params.genome}/*.gtf", checkIfExists: true)
            .set{ path_annotation }
    }
//    path_genome.view()
//    path_annotation.view()

    // Create genome index
    path_index = (
        params.index == null ?
        Index(path_genome, path_annotation) :
        Channel.fromPath("${params.index}", checkIfExists:true)
    )
//    path_index.view()

    //fastq_files.view()
//    mapping_path = Mapping(fastq_files, path_index)
    mapping_path = (
        params.mapping == null ?
        Mapping(fastq_files, path_index) :
        Channel.fromPath("${params.mapping}/*.bam", checkIfExists:true)
    )
    mapping_path.view()
    path_annotation.view()

    // Create counting matrix
    counting_path = Counting(path_annotation,mapping_path.toSortedList())
    counting_path.view()

}

// featureCounts -T CPUS -t gene -g gene_id -s 0 -a annotations.gtf -o counts.txt -M bam1 bam2 bam3 bam4
