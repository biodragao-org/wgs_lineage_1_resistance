#!/usr/bin/env nextflow


/*
################
NEXTFLOW Global Config
################
*/

ch_refGbk = Channel.value("$baseDir/NC000962_3.gbk")
ch_refFasta = Channel.value("$baseDir/NC000962_3.fasta")

/*
################
download genomes
################
*/

Channel.fromFilePairs("./*_{1,2}.fastq.gz", flat: true)
        .into { ch_fastqGz; ch_snippy }


/*
#==============================================
# gzip
#==============================================
*/



process gzip {
    container 'abhi18av/biodragao_base'
    publishDir 'results/gzip'
//    echo true

    input:
    set genomeName, file(genomeReads) from ch_gzip

    output:
    tuple genomeName, path(genome_1_fq), path(genome_2_fq) into ch_trimmomatic
    tuple genomeName, path(genome_1_fq), path(genome_2_fq) into ch_in_rdAnalyzer

    script:
    genome_1_fq = genomeReads[0].name.split("\\.")[0] + '.fastq'
    genome_2_fq = genomeReads[1].name.split("\\.")[0] + '.fastq'

    """
    gzip -dc ${genomeReads[0]} > ${genome_1_fq}
    gzip -dc ${genomeReads[1]} > ${genome_2_fq}

    """
}





/*
#==============================================
# mtbseq
quay.io/biocontainers/mtbseq:1.0.4--pl526_0
#==============================================
*/




/*
#==============================================
# TODO mykrobe
quay.io/biocontainers/mykrobe:0.8.1--py37ha80c686_0
#==============================================
*/




/*
#==============================================
# TODO quast
# quay.io/biocontainers/quast:5.0.2--1
#==============================================
*/

