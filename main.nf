#!/usr/bin/env nextflow


/*
#==============================================
# config
#==============================================
*/

params.outdir = "results"
ch_refGbk = Channel.value("$baseDir/NC000962_3.gbk")
ch_refFasta = Channel.value("$baseDir/NC000962_3.fasta")


/*
#==============================================
# read genomes
#==============================================
*/

Channel.fromFilePairs("./*_{R1,R2}.fastq.gz")
        .into { ch_gzip; ch_in_tbProfiler; ch_in_snippy }





/*
#==============================================
# snippy
#==============================================
*/

process snippy {

    container 'quay.io/biocontainers/snippy:4.6.0--0'
    publishDir 'results/snippy_results'

    input:
    path refGbk from ch_refGbk
    input:
    set genomeName, file(genomeReads) from ch_in_snippy

    script:
    read_1_gz = genomeReads[0]
    read_2_gz = genomeReads[1]

    """
    snippy --cpus 4 --outdir $genomeName --ref $refGbk --R1 $read_1_gz --R2 $read_2_gz
    """
}




/*
#==============================================
# TB_Profiler
#==============================================
*/


process tbProfiler {

    container 'quay.io/biocontainers/tb-profiler:2.8.6--pypy_0'
    publishDir 'results/tbProfiler_results'

    input:
    set genomeName, file(genomeReads) from ch_in_tbProfiler

    script:
    read_1_gz = genomeReads[0]
    read_2_gz = genomeReads[1]

    """
    tb-profiler profile -1 $read_1_gz -2 $read_2_gz  -t 4 -p $genomeName
    """
}



/*
#==============================================
# gzip
#==============================================
*/

process gzip {
    container 'abhi18av/biodragao_base'
    publishDir 'results/gzip_results'

    input:
    set genomeName, file(genomeReads) from ch_gzip

    output:
    tuple genomeName, path(genome_1_fq), path(genome_2_fq) into ch_in_trimmomatic
    tuple genomeName, path(genome_1_fq), path(genome_2_fq) into ch_in_spades

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
trimmomatic
#==============================================
*/


process trimmomatic {
    publishDir 'results/trimmomatic_results' // TODO add regexp to only publish the *paired* files
    container 'quay.io/biocontainers/trimmomatic:0.35--6'


    input:
    tuple genomeName, path(fq_1), path(fq_2) from ch_in_trimmomatic

    output:
    tuple  genomeName, path(fq_1_paired), path(fq_2_paired) into ch_in_fastqc

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
#==============================================
# fastqc
#==============================================
*/

process fastqc {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    publishDir 'results/fastqc_results'

    input:
    tuple  genomeName, path(fq_1_paired), path(fq_2_paired) from ch_in_fastqc

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
# spades
#==============================================
*/

process spades {
    container 'quay.io/biocontainers/spades:3.14.0--h2d02072_0'
    publishDir 'results/spades_results'

    input:
    tuple genomeName, path(fq_1), path(fq_2) from ch_in_spades

    output:
    path """${genomeName}_spades/scaffolds.fasta""" into ch_in_prokka


    script:

    """
    spades.py -k 21,33,55,77 --careful --only-assembler --pe1-1 ${fq_1} --pe1-2 ${fq_2} -o ${genomeName}_spades -t 2
    """
}


/*
#==============================================
# prokka
#==============================================
*/

process prokka {
    container 'quay.io/biocontainers/prokka:1.14.6--pl526_0'
    publishDir 'results/prokka_results'

    input:
    path refFasta from ch_refFasta
    path bestContig from ch_in_prokka

    script:
    genomeName = bestContig.getName().split("\\_")[0]
    contigName = bestContig + "_NC000962_3.fasta.fasta"

//prokka --outdir 23_prokka --prefix 23 23_scaffolds.fasta
    """
    prokka --outdir ./${genomeName}_prokka --prefix $genomeName $contigName
    """
}



/*
#==============================================
# mtbseq
quay.io/biocontainers/mtbseq:1.0.4--pl526_0
#==============================================
*/

/*
// in a folder with all fqgz files
MTBseq --step TBfull 
*/


/*
#==============================================
# TODO mykrobe
quay.io/biocontainers/mykrobe:0.8.1--py37ha80c686_0
#==============================================
*/

/*
mykrobe predict G04868 tb --output G04868_mykrobe.csv --seq G04868_MTBSeq_nextseq_151bp_R1.fastq.gz G04868_MTBSeq_nextseq_151bp_R2.fastq.gz

*/


/*
#==============================================
# TODO quast
# quay.io/biocontainers/quast:5.0.2--1
#==============================================
*/

/*
   cp ../spades_results/work/ae/4d2677d866d9030d613d882de8b639/210_spades/scaffolds.fasta 210_scaffolds.fasta
 
   cp ../spades_results/work/b0/a4240ada1affe5b4217da78f21bc8f/23_spades/scaffolds.fasta 23_scaffolds.fasta

quast.py 23_scaffolds.fasta 210_scaffolds.fasta

*/
