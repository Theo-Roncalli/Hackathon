nextflow.enable.dsl=2

//params.idSRA = ["SRR15678351","SRR15678351"]
params.fastq = null
// params.idSRA = ["SRR628582", "SRR628583"]
params.idSRA = ['SRR628582', 'SRR628583', 'SRR628584', 'SRR628585', 'SRR628586', 'SRR628587', 'SRR628588', 'SRR628589']

process Fasterq {

    tag "Importation of ${idSRA}"

    input:
        val idSRA

    output:
        path "*_{1,2}.fastq"

    script:
    """
    fasterq-dump ${idSRA}
    """
}

process Import {

    tag "Scan of fasterq files"

    input:
        path params.fastq

    output:
        stdout

    script:
    """

    """
}

workflow {
    if (params.fastq == null){
//      idSRA = Channel.of("SRR628582", "SRR628583")
        idSRA = Channel.fromList(params.idSRA).flatten()
        fasterq_files = Fasterq(idSRA)
        fasterq_files.view()
        println fasterq_files.getClass()
        println fasterq_files
    }
    else if (params.fastq instanceof String){
        fasterq_files = Import(params.fastq)
/*        x = Channel.fromList([["/home/ubuntu/Documents/AMI2B/Hackathon/Projet/Nextflow/work/6c/5028f23ac1cb97bb5e679417ef1057/SRR15678351_1.fastq",
        "/home/ubuntu/Documents/AMI2B/Hackathon/Projet/Nextflow/work/6c/5028f23ac1cb97bb5e679417ef1057/SRR15678351_2.fastq"],
        ["/home/ubuntu/Documents/AMI2B/Hackathon/Projet/Nextflow/work/bd/064b3faf470aa8cdc4db2b0a2700c9/SRR15678351_1.fastq",
        "/home/ubuntu/Documents/AMI2B/Hackathon/Projet/Nextflow/work/bd/064b3faf470aa8cdc4db2b0a2700c9/SRR15678351_2.fastq"]])
        x.view()
        println x.getClass()*/
    }

}