#!/usr/bin/env nextflow


/*
################
NEXTFLOW Global Config
################
*/

/*

docker pull quay.io/biocontainers/mtbseq:1.0.4--pl526_0
docker pull quay.io/biocontainers/mykrobe:0.8.1--py37ha80c686_0
*/


genomeIDs = [
'ERR036201',
'ERR036213',
]



ch_refGbk = Channel.value("$baseDir/NC000962_3.gbk")
ch_refFasta = Channel.value("$baseDir/NC000962_3.fasta")

/*
################
download genomes
################
*/


/*
--------------------------------------------
 if the genomes are available on a cloud drive
--------------------------------------------
*/

//process rclone {
//
//    output:
//    stdout downloadGenomes_result
//
//    shell:
//
//    """
//rclone copy onedrive-emab:G04868-done/G04868_analysis/G04868_1.fastq.gz . -vv
//rclone copy onedrive-emab:G04868-done/G04868_analysis/G04868_2.fastq.gz . -vv
//    """
//}
//

/*
--------------------------------------------
 if the genomes are available locally
--------------------------------------------
*/

//Channel.fromFilePairs("./*_{1,2}.fastq.gz", flat: true)
//        .into { ch_fastqGz; ch_snippy }


/*
--------------------------------------------
 if the genomes are available on NCBI
--------------------------------------------
*/

ids = [
'ERR036201',
'ERR036213',
'ERR036217',
]

//Channel
//    .fromSRA(ids)
//    .set { ch_gzip }


Channel
    .fromSRA(ids, cache: true)
    .into{ch_gzip; ch_tbProfiler_in}


/*
################
gzip these files
################
*/

process gzip {
//    echo true

    publishDir 'results/gzip'
    container 'abhi18av/biodragao_base'

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

//###############
//snippy_command
//###############
//*/
//
//process runSnippy {
//
//    container 'quay.io/biocontainers/snippy:4.6.0--0'
//    publishDir 'results/snippy', pattern: './*vcf.gz'
//
//    input:
//    path refGbk from ch_refGbk
//    tuple (genomeInfo, path(read_1_gz), path(read_2_gz)) from ch_snippy
//
//    script:
//    genomeName = read_1_gz.name.split("\\.")[0].split("\\_")[0]
//
//    """
//    snippy --cpus 4 --outdir $genomeName --ref $refGbk --R1 $read_1_gz --R2 $read_2_gz
//    """
//}





/*
###############
trimmomatic
###############
*/


process runTrimmomatic {
//    echo true
// TODO add regexp to only publish the *paired* files
    publishDir 'results/trimmomatic'
    container 'quay.io/biocontainers/trimmomatic:0.35--6'


    input:
    tuple genomeName, path(fq_1), path(fq_2) from ch_trimmomatic

    output:
    tuple genomeName, path(fq_1), path(fq_1_paired), path(fq_1_unpaired), path(fq_2), path(fq_2_paired), path(fq_2_unpaired) into pairedCh
    tuple genomeName, path(fq_1_paired), path(fq_1_unpaired), path(fq_2_paired), path(fq_2_unpaired) into(ch_velvet_41, ch_velvet_49, ch_velvet_55)
    val genomeName into ch_genomeName
    tuple genomeName, path(fq_1_paired), path(fq_2_paired) into ch_fastqc
    tuple genomeName, path(fq_1_paired) into ch_in_spotyping

    script:

    fq_1_paired = genomeName + '_1_paired.fastq'
    fq_1_unpaired = genomeName + '_1_unpaired.fastq'
    fq_2_paired = genomeName + '_2_paired.fastq'
    fq_2_unpaired = genomeName + '_2_unpaired.fastq'

    """
    trimmomatic \
    PE -phred33 \
    $fq_1 \
    $fq_2 \
    $fq_1_paired \
    $fq_1_unpaired \
    $fq_2_paired \
    $fq_2_unpaired \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}


/*
###############
fastqc
###############
*/





process fastqcAfterTrimmomatic {
    publishDir 'results/fastqcAfterTrimmomatic'
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    input:
    tuple  genomeName, path(fq_1_paired), path(fq_2_paired) from ch_fastqc

    output:
    file("fastqc_${genomeName}_logs") into ch_multiqc

    script:
    """
    mkdir fastqc_${genomeName}_logs
    fastqc -o fastqc_${genomeName}_logs ${fq_1_paired}
    fastqc -o fastqc_${genomeName}_logs ${fq_2_paired}
    """
}


/*
#==============================================
# TB_Profiler
#==============================================
*/

process tbProfiler {
//    publishDir path: output_tbProfiler, pattern: "./results/*json", mode: "copy"
//    conda conda_emilyn_py3

//    publishDir 'results/tbProfiler'
    container 'quay.io/biocontainers/tb-profiler:2.8.6--pypy_0'

    input:
    set genomeName, file(genomeReads) from ch_tbProfiler_in

    script:
    read_1_gz = genomeReads[0]
    read_2_gz = genomeReads[1]

    """
    tb-profiler profile -1 $read_1_gz -2 $read_2_gz  -t 4 -p $genomeName
    """
}

/*
#==============================================
# RD_Analyzer
#==============================================
*/

process rdAnalyzer {
//    publishDir path: output_rdAnalyzer, pattern: "*result", mode: "copy"
//    conda conda_emilyn_py2

//    publishDir 'results/rdAnalyzer'
    container 'abhi18av/rdanalyzer'


    input:
    tuple genomeName, path(fq_1), path(fq_2) from ch_in_rdAnalyzer

    script:

    """
    python  /RD-Analyzer/RD-Analyzer.py  -o ${genomeName} ${fq_1} ${fq_2}
    """
}


/*
#==============================================
# Spotyping
#==============================================
*/

process spotyping {
//    publishDir path: output_spotyping, pattern: "*txt", mode: "copy"
//    publishDir path: output_spotyping, pattern: "*xls", mode: "copy"
//    conda conda_emilyn_py2


//    publishDir 'results/spotyping'
    container 'abhi18av/spotyping'


    input:
    tuple genomeName, path(fq_1_paired) from ch_in_spotyping

    script:

    """
    python /SpoTyping-v2.0/SpoTyping-v2.0-commandLine/SpoTyping.py ${fq_1_paired} -o ${genomeName}.txt
    """
}
