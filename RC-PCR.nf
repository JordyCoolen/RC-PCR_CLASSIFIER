#!/usr/bin/env nextflow
params.threads = 4
params.outDir = "./output"
params.reads = "$baseDir/test/test_OUT01_R{1,2}.fastq.gz"
params.UMILEN = 18
params.minreadlength = 100
params.database = "SILVA"
params.abricate = true

// Parsing the input parameters
outDir           = "$params.outDir"
threads          = "$params.threads"
UMILEN           = "$params.UMILEN"
minreadlength    = "$params.minreadlength"
database         = "$params.database"
primerfile       = "${baseDir}/db/${database}/primers/${database}_primers.fasta"
bedfile          = "${baseDir}/db/${database}/primers/${database}_primers.bed"
KMAdb            = "${baseDir}/db/${database}/KMA/${database}"
blastdbpath      = "${baseDir}/db/${database}/"
STARfasta        = "${baseDir}/db/${database}/STAR/CYP51A.fa"
def samplename   = file("$params.reads").simpleName[0].split('_')[0]
abricate         = "$params.abricate"

// Tools paths and command prefixes
reporter         = "$baseDir/scripts/final_report.py"
vcf2table        = "$baseDir/scripts/vcf2table.py"

// Parsing the input parameters
outDir           = "$params.outDir"+"/${samplename}"
threads          = "$params.threads"

// special channel for the fastQ reads
Channel
      .fromFilePairs( params.reads )
      .ifEmpty { "cannot find read pairs in path"}
      .set  { reads_ch1 }

log.info """

NEXTFLOW SILVA Classification RC-PCR V0.1
================================
samplename : $samplename
reads      : $params.reads
outDir     : $params.outDir
codeBase   : $baseDir
threads    : $params.threads
abricate   : $params.abricate
UMILEN     : $params.UMILEN
minreadlenth : $params.minreadlength
database   : $params.database

~~~~~~~~~~~Databases~~~~~~~~~~~
primerfile : $primerfile
bedfile    : $bedfile
blastdb    : $blastdbpath+"blast"
KMAdb      : $KMAdb

~~~~~~~~~~~Authors~~~~~~~~~~~~~~
        J.P.M. Coolen
        B. van Wessel
================================
"""

// Clean reads (adapter and read length filter)
process '1A_clean_reads' {
    tag '1A'
    conda 'bioconda::fastp=0.20.1 bioconda::pyfastx=0.6.12 conda-forge::simplejson=3.17.0'
    publishDir outDir + '/fastp', mode: 'copy'
    input:
        set pairID, file(reads) from reads_ch1
    output:
        set file("${samplename}_R1_fastp.fastq.gz"), file("${samplename}_R2_fastp.fastq.gz") into fastp_2A, fastp_3A, fastp_3C, fastp_6
        file "${samplename}.fastp.json"
        file "${samplename}.fastp.html"
        file ".command.*"
    script:
        """
        fastp -i ${reads[0]} -o ${samplename}_R1_fastp.fastq.gz \
                                -I ${reads[1]} -O ${samplename}_R2_fastp.fastq.gz \
                                --umi_len=${UMILEN} --umi --umi_loc=per_read --umi_prefix=UMI \
                                --html ${samplename}.fastp.html --json ${samplename}.fastp.json \
                                --length_required ${minreadlength} --trim_poly_g --trim_poly_x
        """
}

// Clean reads (adapter and read length filter)
process '2A_measure_amplicons' {
    tag '1A'
    conda 'conda-forge::pandas=1.2.4 bioconda::pysam=0.15.3 anaconda::openpyxl'
    publishDir outDir, mode: 'copy'
    input:
        set file(forward_read), file(reverse_read) from fastp_2A
    output:
        file "*.csv"
        file ".command.*"
    script:
        """
        python $baseDir/scripts/main.py --input $forward_read --primers ${primerfile}
        """
}

// create KMA tool to detect 16S
// Process 3A: KMA
process '3A_KMA' {
  tag '3A'
  time "30m"
  conda 'bioconda::kma=1.3.28'
  publishDir outDir + '/kma', mode: 'copy'
  input:
  file reads from fastp_3A
  output:
    file "${samplename}*"
    file "${samplename}.fsa" into consensus_4A
    file "${samplename}.res" into kma_5A, kma_7B
    file ".command.*"
  script:
    """
    kma -t_db ${KMAdb} -ipe ${reads[0]} ${reads[1]} -t ${threads} \
    -ef -ex_mode -1t1 -and -apm f -o ${samplename} 2>/dev/null || exit 0
    #kma -t_db ${KMAdb} -ipe ${reads[0]} ${reads[1]} -t ${threads} \
    #-a -ex_mode -ef -1t1 -and -apm f -o ${samplename} -sam 4 > ${samplename}.sam 2>/dev/null || exit 0
    """
}

// create abricate to detect 16S/18S SILVA
// Process 4A: abricate
process '4A_abricate' {
  tag '4A'
  conda 'bioconda::abricate=1.0.1'
  publishDir outDir + '/abricate', mode: 'copy'
  input:
    file consensus from consensus_4A
  output:
    file "${samplename}_blast.txt" into blast_7B
    file ".command.*"
  script:
    if(abricate==true)
    """
    abricate --datadir ${blastdbpath} --db blast ${consensus} --mincov 30 --minid 60 --threads ${threads}  > ${samplename}_blast.txt
    """
    else if(abricate==false)
    """
    echo 'none' > ${samplename}_blast.txt
    """
}

// SILVA database circle packing visualization using nodejs
process '5A_circle_packing_viz' {
    tag '5A'
    conda "${baseDir}/conda/env-nodejs"
    publishDir outDir + '/report/viz', mode: 'copy'
    input:
        file kma from kma_5A
    output:
        file ".command.*"
    script:
        """
        python ${baseDir}/conda/env-nodejs/Circle-packing-visualization/data-processing/circle-packing-parsing.py -i "${outDir}/kma" \
        -o ${baseDir}/conda/env-nodejs/Circle-packing-visualization/Data.csv

        cd ${baseDir}/conda/env-nodejs/Circle-packing-visualization/

        #excute nodejs code
        npm run build

        # copy files manual to report folder
        mkdir -p ${outDir}/report/viz
        cp index.html ${outDir}/report/viz
        cp -r dist ${outDir}/report/viz/
        cp style.css ${outDir}/report/viz
        cp ${baseDir}/conda/env-nodejs/Circle-packing-visualization/Data.csv ${outDir}/report/viz
        """
}

// Process 6: multiQC
process '6_multiQC' {
  tag '6'
  conda 'bioconda::multiqc'
  publishDir outDir + '/QC', mode: 'copy'
  input:
  file reads from fastp_6
  output:
    file "*.html"
    file ".command.*"
  script:
    """
    multiqc ${outDir}/fastp/
    """
}

// Process 7A: obtain run parameters
process '7A_parameters' {
    tag '7A'
    publishDir outDir + '/report', mode: 'copy'
    input:
    output:
        file "parameters.txt" into params_7B
    script:
        """
        touch parameters.txt
        echo "Parameter\tValue" >> parameters.txt
        echo "Reads\t$params.reads" >> parameters.txt
        echo "Database\t$params.database" >> parameters.txt
        echo "abricate:\t$params.abricate" >> parameters.txt
        echo "minreadlength:\t$params.minreadlength" >> parameters.txt
        echo "UMI length:\t$params.UMILEN" >> parameters.txt
        """
}

// Process 7B: generate a report for interpretation by the clinician (or for research purposes)
process '7B_report' {
    tag '7B'
    conda "${baseDir}/conda/env-025066a104bf8ce5621e328d8009733a"
    publishDir outDir + '/report', mode: 'copy'
    input:
        file params from params_7B
        file blast from blast_7B
        file kma from kma_7B
    output:
        file "${samplename}.html"
        file "${samplename}.pdf"
        file ".command.*"
    script:
        """
        $reporter --sampleName ${samplename} --params ${params} --blast ${blast} \
        --kma ${kma}
        """
}