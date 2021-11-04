// nextflow run main.nf --fastq ../Data

nextflow.enable.dsl=2

params.SRAPath = null
params.genomePath = null
params.idSRA = ['SRR15678351','SRR15289297']
params.ftpGenome = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.14_GRCh37.p13/GO_TO_CURRENT_VERSION/GCA_000001405.28_GRCh38.p13_genomic.fna.gz"
// params.idSRA = ['SRR628582', 'SRR628583', 'SRR628584', 'SRR628585', 'SRR628586', 'SRR628587', 'SRR628588', 'SRR628589']

process Fasterq {

    tag "Importation of ${idSRA}"

    input:
    val idSRA

    output:
    tuple val("${idSRA}"), path("*_{1,2}.fastq")

    script:
    """
    fasterq-dump ${idSRA}
    """
}

process Genome {

    tag "Importation of ${ftp}"

    input:
    val ftp

    output:
    path("*.fna")

    script:
    """
    wget ${ftp}
    gunzip *.f*a.gz
    """
}

workflow {

    if (params.SRAPath == null){
        idSRA = Channel.fromList(params.idSRA)
        fasterq_files = Fasterq(idSRA)
//        fasterq_files.view()
    }
    else if (params.SRAPath instanceof String){
        fasterq_files = Channel.fromFilePairs("${params.SRAPath}/SRR*_{1,2}.fastq",checkIfExists:true)
        fasterq_files.view()
    }
    else {
        throw new Exception("The path for fastq files is not a string.")
    }

    if (params.genomePath == null){
        ftp = Channel.value(params.ftpGenome)
        genome_file = Genome(ftp)
        genome_file.view()
    }
    else if (params.genomePath instanceof String){
        genome_file = Channel.fromPath("${params.genomePath}/*.f*a",checkIfExists:true)
        genome_file.view()
    }
    else {
        throw new Exception("The path for genome file is not a string.")
    }

}

