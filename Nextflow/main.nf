nextflow.enable.dsl=2

params.input = "data/*.fastq"
params.ids = "SRR15678351"
//params.ids = "SRR628582"
// params.ids = ['SRR628582', 'SRR628583', 'SRR628584', 'SRR628585', 'SRR628586', 'SRR628587', 'SRR628588', 'SRR628589']

process Fasterq {

    tag "File ${ids} importation"

    input:
    val ids

//    output:
//    path "${projectDir}"

    script:
    """
    echo ${workDir}
    fasterq-dump -q ${ids}
    echo ${workDir}
    """


}

workflow {
    println "workflow launched."
    ids = Channel.of(params.ids)
    Fasterq(ids)

}