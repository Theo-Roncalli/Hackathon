#!/usr/bin/env nextflow


str = Channel.from('SRR15678351','SRR15289297')

//Contrairement à mon camarade théo , je n'utilise pas le DSL 2 .


process Download {
    
    input :
    env SRR from str

    
    output : 
    file "*_{1,2}.fastq" into result
    
    '''
    fasterq-dump $SRR     
    '''

}
    
process IndexCreation {

    input:
    path x from '/home/izem/Desktop/Master/Master1/Nextflow/hg19.fa'
    
    """
    mkdir genomeDir
      
    STAR --runThreadN 16 --runMode genomeGenerate --genomeDir genomeDir --genomeFastaFiles $x

    """
}




