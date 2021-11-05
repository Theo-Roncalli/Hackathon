// nextflow run main.nf --reads ../Data/Reads --genome ../Data/Genome/GRCh38.primary_assembly.genome.fa

nextflow.enable.dsl=2

params.reads = null
params.genome = null
params.ids = ['SRR15678351','SRR15289297']
params.ftp = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.14_GRCh37.p13/GO_TO_CURRENT_VERSION/GCA_000001405.28_GRCh38.p13_genomic.fna.gz"
// params.ids = ['SRR628582', 'SRR628583', 'SRR628584', 'SRR628585', 'SRR628586', 'SRR628587', 'SRR628588', 'SRR628589']

process Fasterq {

    tag "Importation of ${ids}"

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

    tag "Importation of ${ftp}"

    input:
    val ftp

    output:
    path "*.fna"

    script:
    """
    wget ${ftp}
    gunzip *.f*a.gz
    """
}

process Index {

    tag "Creation of the index"
    //conda "STAR"

    input:
        path genome_file

    output:
        path "Index/"
        

    script:
    """
    STAR --runThreadN ${params.index_cpus} --runMode genomeGenerate –genomeDir Index --genomeFastaFiles ${genome_file}
    """
}

workflow {

    // SEE EXPLANATION OF THE NEW PROGRAM STRUCTURE AFTER THE WORKFLOW

// START getting entry files
    ids = Channel.fromList(params.ids).ifEmpty(null)
    fasterq_files = (
        ids == null ?
        Channel.fromFilePairs("${params.reads}/SRR*_{1,2}.fastq", checkIfExists:true)
        : Fasterq(ids)
    )
    //fasterq_files.view()
// END getting entry files

// START getting genome
    ftp = Channel.value(params.ftp).ifEmpty(null)
    genome_file = (
        ftp == null ?
        Channel.fromPath("${params.genome}", checkIfExists:true)
        : Genome(ftp)
    )
    //genome_file.view()
// END getting genome

// START creating genome index
    path_index = Index(genome_file)
    //path_index.view()
// END creating genome index

}

// NEW PROGRAM STRUCTURE
/*
    Instead of a series of if statements (imperative programming)
    to build the parameters, channels are built with the ternary 
    operation which allows yielding two different values based on a condition.

    The reformulation is as follows

    Original :
        if some_condition:
            x = 5
        else:
            x = 6

    Explicit ternary operator (valid Python, try it for yourself):
        x = 5 if some_condition else 6

    C-style (valid in Groovy) ternary operator
        x = some_condition ? 5 : 6
    
    It's basically asking a question, with the convention that the first 
    value after the question mark is the "yes" and the second (after the colon :)
    is "no".

*/


// MANUAL TESTING COMMANDS HISTORY :

// docker pull combinelab/salmon
// wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz
// wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
// cat gencode.v29.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz
// salmon index -t gentrome.fa.gz --decoys GRCh38.primary_assembly.genome.fa.gz -p 12 -i salmon_index --gencode

// gunzip *.f*a.gz
// salmon index -t gencode.v29.transcripts.fa.gz -i index --gencode

// STAR --runThreadN 6 --runMode genomeGenerate –genomeDir Index --genomeFastaFiles Genome/GRCh38.primary_assembly.genome.fa
