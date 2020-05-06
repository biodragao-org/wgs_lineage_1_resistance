#!/usr/bin/env nextflow


/*
################
NEXTFLOW Global Config
################
*/

/*
created two conda envs for this analysis

# this version doesn't exist
conda create --name py254 python=2.5.4
conda_emilyn_py2 = "/data1/users/emilyn.costa/.conda/envs/emilyn-py2"

# https://askubuntu.com/questions/204510/how-to-install-python-2-5-4
# http://ubuntuhandbook.org/index.php/2013/08/install-python-2-6-2-5-3-3-in-ubuntu-12-04-12-10-13-04/

sudo apt-add-repository ppa:fkrull/deadsnakes
sudo apt-get update
sudo apt-get install python2.5

*/


malawiLineage1_genomeIDs = [
'ERR036201',
'ERR036213',
'ERR036217',
'ERR036222',
'ERR036223',
'ERR036230',
'ERR036236',
'ERR037473',
'ERR037476',
'ERR037489',
'ERR037491',
'ERR037496',
'ERR037510',
'ERR037511',
'ERR037519',
'ERR037520',
'ERR037521',
'ERR037523',
'ERR037536',
'ERR037555',
'ERR124640',
'ERR126602',
'ERR126605',
'ERR126607',
'ERR126621',
'ERR126631',
'ERR126636',
'ERR126639',
'ERR161031',
'ERR161059',
'ERR161063',
'ERR161085',
'ERR161089',
'ERR161100',
'ERR161106',
'ERR161121',
'ERR161124',
'ERR161129',
'ERR161131',
'ERR161133',
'ERR161138',
'ERR161144',
'ERR161161',
'ERR161180',
'ERR161181',
'ERR161190',
'ERR161201',
'ERR163952',
'ERR163963',
'ERR163972',
'ERR163980',
'ERR163982',
'ERR163994',
'ERR164000',
'ERR164002',
'ERR164005',
'ERR164019',
'ERR176453',
'ERR176456',
'ERR176476',
'ERR176483',
'ERR176485',
'ERR176490',
'ERR176491',
'ERR176517',
'ERR176524',
'ERR176564',
'ERR176570',
'ERR176588',
'ERR176597',
'ERR176602',
'ERR176618',
'ERR176622',
'ERR176630',
'ERR176637',
'ERR176654',
'ERR176661',
'ERR176677',
'ERR176689',
'ERR176712',
'ERR176726',
'ERR176735',
'ERR176747',
'ERR176757',
'ERR176770',
'ERR176775',
'ERR176788',
'ERR176795',
'ERR176800',
'ERR176819',
'ERR176823',
'ERR176827',
'ERR181680',
'ERR181682',
'ERR181687',
'ERR181723',
'ERR181724',
'ERR181725',
'ERR181728',
'ERR181739',
'ERR181760',
'ERR181775',
'ERR181778',
'ERR181779',
'ERR181793',
'ERR181794',
'ERR181795',
'ERR181798',
'ERR181799',
'ERR181807',
'ERR181819',
'ERR181820',
'ERR181831',
'ERR181844',
'ERR181860',
'ERR181862',
'ERR181907',
'ERR181910',
'ERR181912',
'ERR181920',
'ERR181929',
'ERR181930',
'ERR181936',
'ERR181949',
'ERR181951',
'ERR181961',
'ERR181963',
'ERR181969',
'ERR181971',
'ERR181972',
'ERR182008',
'ERR182009',
'ERR182010',
'ERR182017',
'ERR182024',
'ERR182042',
'ERR182054',
'ERR190336',
'ERR190337',
'ERR190344',
'ERR190345',
'ERR190349',
'ERR190353',
'ERR190361',
'ERR190365',
'ERR190377',
'ERR190382',
'ERR190389',
'ERR190403',
'ERR211991',
'ERR212022',
'ERR212038',
'ERR212048',
'ERR212054',
'ERR212066',
'ERR212067',
'ERR212073',
'ERR212079',
'ERR212084',
'ERR212088',
'ERR212089',
'ERR212090',
'ERR212093',
'ERR212094',
'ERR212096',
'ERR212097',
'ERR212102',
'ERR212103',
'ERR212105',
'ERR212106',
'ERR212108',
'ERR212109',
'ERR212110',
'ERR212111',
'ERR212113',
'ERR212116',
'ERR212131',
'ERR212136',
'ERR212137',
'ERR212150',
'ERR212155',
'ERR212157',
'ERR212160',
'ERR212163',
'ERR212164',
'ERR212166',
'ERR212178',
'ERR212179',
'ERR212180',
'ERR216902',
'ERR216903',
'ERR216904',
'ERR216927',
'ERR216940',
'ERR216976',
'ERR216988',
'ERR221557',
'ERR221562',
'ERR221563',
'ERR221569',
'ERR221570',
'ERR221578',
'ERR221582',
'ERR221593',
'ERR221596',
'ERR221606',
'ERR221607',
'ERR245647',
'ERR245648',
'ERR245651',
'ERR245652',
'ERR245654',
'ERR245656',
'ERR245657',
'ERR245662',
'ERR245670',
'ERR245672',
'ERR245674',
'ERR245677',
'ERR245679',
'ERR245691',
'ERR245695',
'ERR245698',
'ERR245699',
'ERR245703',
'ERR245705',
'ERR245717',
'ERR245720',
'ERR245721',
'ERR245727',
'ERR245734',
'ERR245739',
'ERR245741',
'ERR245748',
'ERR245751',
'ERR245760',
'ERR245761',
'ERR245763',
'ERR245775',
'ERR245781',
'ERR245782',
'ERR245784',
'ERR245811',
'ERR245813',
'ERR245817',
'ERR245820',
'ERR245822',
'ERR245825',
'ERR245829',
'ERR245836',
'ERR245843',
'ERR245844',
'ERR245846',
'ERR323051',
'ERR323079',
'ERR323092',
'ERR323093',
'ERR323106',
'ERR323111',
'ERR473282',
'ERR473297',
'ERR473299',
'ERR473305',
'ERR473320',
'ERR473345',
'ERR473349',
'ERR473353',
'ERR736807',
'ERR736809',
'ERR736820',
'ERR736821',
'ERR736827',
'ERR736828',
'ERR736844',
'ERR736860',
'ERR736866',
'ERR736868',
'ERR736869',
'ERR736879',
'ERR773789',
'ERR773793',
'ERR773799',
'ERR773800',
'ERR773802',
'ERR773810',
'ERR779656',
'ERR779657',
'ERR779660',
'ERR779664',
'ERR779679'
]



params.outdir = "results"

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
