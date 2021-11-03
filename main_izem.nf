nextflow.enable.dsl=2

params.ids = "SRR15678351"
// params.ids = ['SRR628582', 'SRR628583', 'SRR628584', 'SRR628585', 'SRR628586', 'SRR628587', 'SRR628588', 'SRR628589']

process Fasterq {

    tag "Importation of the file ${ids}"

    input:
    val ids

    script:
    """
    fasterq-dump ${ids}
    """

}




process MapThemAll {
  
   tag "Creation of index ..."

    '''
    
    mkdir genomeDir

    STAR --runThreadN 3 \
    --runMode genomeGenerate \
    --genomeDir genomeDir \
    --genomeFastaFiles hg19.fa
    
    '''      

}


workflow {
    println "workflow launched."
    ids = Channel.of(params.ids)
    Fasterq(ids)
    println "Second Step , Creating index"
    MapThemAll()
    println "Let's map it ! "

}
