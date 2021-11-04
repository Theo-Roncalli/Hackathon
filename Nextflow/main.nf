nextflow.enable.dsl=2

params.idSRA = ['SRR15678351','SRR15289297']
params.fastq = null
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

workflow {
    if (params.fastq == null){
        idSRA = Channel.fromList(params.idSRA)
        fasterq_files = Fasterq(idSRA)
//        fasterq_files.view()
    }
    else if (params.fastq instanceof String){
        fasterq_files = Channel.fromFilePairs("${params.fastq}/SRR*_{1,2}.fastq",checkIfExists:true)
//        fasterq_files.view()
    }
    else {
        throw new Exception("fastq parameter is not a string.")
    }
}

