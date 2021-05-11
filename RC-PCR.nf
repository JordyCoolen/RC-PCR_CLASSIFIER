#!/usr/bin/env nextflow
params.threads = 4
params.outDir = "./output"
params.reads = "$baseDir/test/test_OUT01_R{1,2}.fastq.gz"
params.UMILEN = 18
params.minreadlength = 100
params.database = "18S"
params.abricate = true

// Parsing the input parameters
outDir           = "$params.outDir"
threads          = "$params.threads"
UMILEN           = "$params.UMILEN"
minreadlength    = "$params.minreadlength"
database         = "$params.database"
primerfile       = "${baseDir}/db/primers/${database}_primers.fasta"
KMAdb            = "${baseDir}/db/KMA/${database}"
blastdbpath      = "${baseDir}/db/blast_db/"
def samplename   = file("$params.reads").simpleName[0].split('_')[0]
abricate         = "$params.abricate"

// Parsing the input parameters
outDir           = "$params.outDir"+"/${samplename}"
threads          = "$params.threads"

// special channel for the fastQ reads
Channel
      .fromFilePairs( params.reads )
      .ifEmpty { "cannot find read pairs in path"}
      .set  { reads_ch1 }

log.info """

NEXTFLOW RC-PCR V0.1
================================
samplename : $samplename
reads      : $params.reads
outDir     : $params.outDir
codeBase   : $baseDir
threads    : $params.threads
abricate   : $params.abricate

~~~~~~~~~~~Databases~~~~~~~~~~~
primerfile : $primerfile
blastdb    : $blastdbpath+"/"$database
KMAdb      : $KMAdb

~~~~~~~~~~~Authors~~~~~~~~~~~~~~
        J.P.M. Coolen
            V. Karvink
            T. Baltussen
================================
"""

// Clean reads (adapter and read length filter)
process '1A_clean_reads' {
    tag '1A'
    conda 'bioconda::fastp=0.20.1 bioconda::pyfastx=0.6.12 conda-forge::simplejson=3.17.0'
    publishDir outDir + '/fastp', mode: 'copy', pattern: "*.fastp.json"
    input:
        set pairID, file(reads) from reads_ch1
    output:
        set file("${samplename}_R1_fastp.fastq.gz"), file("${samplename}_R2_fastp.fastq.gz") into fastp_2A, fastp_3A, fastp_6
        file "${samplename}.fastp.json"
        file "${samplename}.fastp.html"
        file ".command.*"
    script:
        """
        fastp -i ${reads[0]} -o ${samplename}_R1_fastp.fastq.gz \
                                -I ${reads[1]} -O ${samplename}_R2_fastp.fastq.gz \
                                --umi_len=${UMILEN} --umi --umi_loc=per_read --umi_prefix=UMI \
                                --html ${samplename}.fastp.html --json ${samplename}.fastp.json \
                                --length_required ${minreadlength}
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
        file "*.xlsx"
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
  conda 'bioconda::kma=1.3.15'
  publishDir outDir + '/kma', mode: 'copy'
  input:
  file reads from fastp_3A
  output:
    file "${samplename}*"
    file "${samplename}.vcf.gz" into vcf_5A
    file "${samplename}.fsa" into consensus_4A
    file "${samplename}.sam" into kma_3B
    file ".command.*"
  script:
    """
    kma -t_db ${KMAdb} -ipe ${reads[0]} ${reads[1]} -t ${threads} -a -ex_mode -ef -dense -1t1 -ref_fsa -vcf 2 -and -apm f -o ${samplename} -sam 4 > ${samplename}.sam 2>/dev/null || exit 0
    #kma -t_db ${KMAdb} -ipe ${reads[0]} ${reads[1]} -t ${threads} -a -ex_mode -ef -dense -1t1 -ref_fsa -mem_mode -and -apm f -o kma 2>/dev/null || exit 0
    """
}

// create KMA tool to match primers
// Process 3B: KMA
process '3B_process_KMA' {
  tag '3B'
  conda 'bioconda::samtools=1.12'
  errorStrategy 'ignore'
  publishDir outDir + '/kma', mode: 'copy'
  input:
  file sam from kma_3B
  output:
    file "${samplename}.sorted.bam"
    file "${samplename}.sorted.bam.bai"
    file ".command.*"
  script:
    """

    # merge sam files

    # sam --> bam
    samtools view -b ${sam} > ${samplename}.bam
    # sort bam
    samtools sort ${samplename}.bam > ${samplename}.sorted.bam
    # index bam
    samtools index ${samplename}.sorted.bam

    # additionally filter reference on only hits
    # this would make it possible to quickly evaluate the results

    """
}

// create abricate to detect 16S
// Process 4A: abricate
process '4A_abricate' {
  tag '4A'
  conda 'bioconda::abricate=1.0.1'
  publishDir outDir + '/abricate', mode: 'copy'
  input:
    file consensus from consensus_4A
  output:
    file "${samplename}_blast.txt"
    file ".command.*"
  script:
    if(abricate==true)
    """
    abricate --datadir ${blastdbpath} --db ${database} ${consensus} --mincov 30 --minid 60 --threads ${threads}  > ${samplename}_blast.txt
    """
    else if(abricate==false)
    """
    echo 'none' > ${samplename}_blast.txt
    """
}

// 5A: annotation of the genome/consensus fasta
process '5A_annotation' {
    tag '5A'
    conda 'bioconda::snpeff=5.0 bioconda::bcftools=1.12'
    publishDir outDir + '/annotation', mode: 'copy'
    input:
        file vcf from vcf_5A
    output:
        file "${samplename}.final.vcf"
        file ".command.*"
  script:
        if(database=="CYP51A")
        """
        bcftools view -f PASS ${vcf} > ${samplename}.pass.vcf
        snpEff CYP51A ${samplename}.pass.vcf -hgvs1LetterAa > ${samplename}.final.vcf
        """
        else if(database!="CYP51A")
        """
        echo 'none' > ${samplename}.final.vcf
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