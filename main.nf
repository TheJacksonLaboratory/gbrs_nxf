#!/usr/bin/env nextflow

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Processes ~ ~ ~ ~ ~
/*
    align fq to "transcripts" bowtie)
    convert to bam (samtools)
    gbrs convert bam to emase (gbrs bam2emase)
    gbrs compress emase file (gbrs bam2emase)
    gbrs merge compress emase file (gbrs compress)
    gbrs quantify and reconstruct emase file (gbrs quantify)
    gbrs quantification round 2 (gbrs quantify, export-genoprob)
    Quality metrics for RNA mapped against B6 with STAR
        (tools: star, picard, starAligner.pl)
*/

import Helpers
import Logos

//FIXME: logo = new Logo()
//println logo.show()

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Parameter Defaults ~ ~ ~ ~ ~ ~
def setParamDefaults() {
    params.help    = false
    params.threads = 1
    params.tmpdir  = "$TMPDIR"
    params.email   = null

    // Configurable variable parameters specific to individual runs:
    params.fastqR1     = null // path to input fastq Read1
    params.fastqR2     = null // path to input fastq Read2
    params.outputDir   = null // Base directory where output will be stored.
    params.generation  = null //TODO: QUESTION: 'DO_GENERATION' needs to have a default value.
    params.sex         = null //TODO: QUESTION: 'SEX' needs to have a default value.

    params.libStrand   = 'NONE'
    params.gbrsSupport = "${projectDir}/suppFiles"
    params.gbrsModel   = 4

}
setParamDefaults()

def helpMessage() {
    log.info"""
    =========================================
      ${workflow.manifest.name} v${workflow.manifest.version}
    =========================================
    ${workflow.manifest.description}

    Usage:
    The typical command for running the pipeline is as follows:
    >   nextflow run ${workflow.projectDir} -params-file path/to/params.yaml -profile sumner -resume
    NOTE: The params.yaml file needs to have the following mandatory parameters
            OR they need to specified on the command line.

    Mandatory:
        --fastqR1 [path]        Paths to input FASTQ Read 1
        --fastqR2 [path]        Paths to input FASTQ Read 2
        --outputDir [path]      Base directory where output will be stored.
        --generation '[G0-G40]' DO generation
        --sex [M|F]             Specify (M)ale of (F)emale

    Optional:
        --gbrsSupport [path]    GBRS support, e.g. bowtie index filesw
                                    (Default: $params.gbrsSupport)
        --gbrsModel             GBRS MODEL (Default: $params.gbrsModel)
                                    [ choices: 1,2,3,4 ]
                                    1: reads are apportioned among genes first,
                                        then between alleles, and then among isoforms.
                                    2: reads are apportioned among genes first,
                                        then among isoforms, and then between alleles.
                                    3: reads are apportioned among genes first,
                                        then among each isoform-allele combination
                                        which are treated equally.
                                    4: assumes no hierarchy and multi-reads are
                                        apportioned equally among genes, isoforms, and
                                        alleles

        --libStrand             Strandness of RNA Library (Default: $params.libStrand)
                                    [ choices: 'NONE', 'FIRST_READ_TRANSCRIPTION_STRAND',
                                               'SECOND_READ_TRANSCRIPTION_STRAND' ]
        --email [email]         The email address to send the pipeline report.
                                    (Default: $params.email)
        --threads [int]         Threads (Must be numeric, Default: $params.threads)
        --tmpdir [int]          Path to hold temporary files for software execution.

        -name [str]             Name for the pipeline run. If not specified Nextflow will
                                    automatically generate a random mnemonic.

        -profile [str]          Environment config to use.
                                    [ choices: standard (local), sumner ]

        --help                  Shows this help message.
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Param File and Format Checking ~ ~ ~ ~ ~ ~
////// Required parameters \\\\\\
if ( ! params.fastqR1 || ! params.fastqR2 ) {
    exit 1, "Parameter ERROR: fastq inputs ($params.fastqR1 & $params.fastqR2) must be specified."
}
if ( ! params.outputDir ) {
    exit 1, "Parameter ERROR: Output directory parameter must be specified."
}
if ( ! file(params.outputDir).exists() ) {
    exit 1, "Parameter ERROR: Missing base root output directory ($params.outputDir) is not found: check if path is correct."
}

//NOTE: below is EXAMPLE:
// Check bowtieIndex, append file extensions to check exists
//if ( ! params.bowtieIndex ) {
//    exit 1, "Parameter ERROR: bowtieIndex absolute path must be specified."
//} else {
//    bowtieFile = params.bowtieIndex + '.1.bt2'
//    if ( ! file(bowtieFile).exists() ) {
//        exit 1, "Parameter ERROR: bowtie2 index file for param ($params.bowtieIndex) does not exist."
//    }
//}

////// Optional params \\\\\\
if ( ! params.gbrsSupport ) {
    exit 1, "Parameter ERROR: Support Files paramete (--gbrsSupport) must be specified."
}
if ( ! file(params.gbrsSupport).exists() ) {
    exit 1, "Parameter ERROR: Support Files directory ($params.gbrsSupport) does not exist."
}

// Check threads is numeric value
def number_vars = [ 'threads' ]
CheckParams.is_numeric(params, number_vars)


//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Summary Info ~ ~ ~ ~ ~ ~
// Header info
def summary = [:]
summary['Pipeline']         = workflow.manifest.name
summary['Description']      = workflow.manifest.description
if(workflow.revision) {
    summary['Pipeline Release'] = workflow.revision
}
summary['Run Name']         = workflow.runName
summary['User']             = workflow.userName
summary['Config Profile']   = workflow.profile
summary['Config Files']     = workflow.configFiles
summary['Command Line']     = workflow.commandLine
summary['Nextflow Info']    = "v${nextflow.version}, build: ${nextflow.build}, on ${nextflow.timestamp}"
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Workflow dir']     = workflow.projectDir
if(workflow.containerEngine) {
    summary['Container Engine'] = "$workflow.containerEngine"
    //summary['Containers'] = "$workflow.container" // only works with Docker. :(
}

// Pipeline Params:
summary['Parameters......']  = ''
if(params.email) {
    summary['.  Email']      = params.email
}
summary['.  Output dir']     = params.outputDir
summary['.  Threads']        = params.threads
summary['.  TMP dir']        = params.tmpdir
summary['.  GBRS Support']   = params.gbrsSupport

summary['.  FASTQ inputs']   = "${params.fastqR1}, ${params.fastqR2}"
summary['.  GBRS Model']     = params.gbrsModel
summary['.  DO Generation']  = params.generation
summary['.  Sex']            = params.sex
summary['.  Library Strand'] = params.libStrand

summary['Run Start Time']    = workflow.start // place last of items

// Print Summary Header:
//FIXME: println Summary.show(summary)
println Summary.show(summary)

//~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Opening Variables and Channels ~ ~ ~ ~ ~
def timestamp = new Date().format("yyyyMMdd'T'hhmmSSS") // `date +"%Y%M%dT%H%M%N"`
def tmpdir    = params.tmpdir

//~~~~~ Determine sampleID and outputDir ~~~~
//NOTE: this techinque b/c only explicit fastq filenames in param, not 'globs' for fromFilePairs

//FIXME: Here we presume '_R1' is used in fastq name! (e.g. not '-R1', '.R1', etc)
def fqPairRE = ~/_R1.*\..*f[ast]*q.*$/
fqR1 = file(params.fastqR1)
fqR2 = file(params.fastqR2)
fqR1path = fqR1.toAbsolutePath().toString() // bc relative paths cause symlink breakage
fqR2path = fqR2.toAbsolutePath().toString() // bc relative paths cause symlink breakage

def sampleID = ( fqR1.name - fqPairRE )


//~~~~~~~~~~~~~~~~ Primary publishDir == sample_outdir ~~~~~
def sample_outdir = "${params.outputDir}/${timestamp}_${sampleID}/"

//FIXME: abspath iff needed....
//~~~~~ Get Abspath of 'outdir' ~~~~~\\
//abs_outdir = file(outdir).toAbsolutePath().toString()


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Extra Reference Files ~~~~~
//FIXME: use ref_dir? or input from gbrsSupport ?
def supp_dir = params.gbrsSupport
ref_tranprob    = "${supp_dir}/tranprob" // directory of tranprob.DO.${generation}.${sex}.npz files
ref_transcripts = "${supp_dir}/transcripts"
ref_trans_info  = "${supp_dir}/ref.transcripts.info"
ref_avecs       = "${supp_dir}/avecs.npz"

ref_GenomeGrid    = "${supp_dir}/ref.genome_grid.69k.noYnoMT_KBEdit.txt"
ref_StarGenomeDir = "${supp_dir}/star_idx_mm10_gencodeM23"
ref_StarIdx       = "${supp_dir}/star_idx_mm10_gencodeM23/gencode.vM23.annotation.gtf"
ref_StarIdxGenome = "${supp_dir}/star_idx_mm10_gencodeM23/mm10p6/GRCm38.p6.genome.fa"

ref_PMetricsInterval = "${supp_dir}/picard_metric_files/gencode.mm10.v23.rRNA.interval_list"
ref_PMetricsFlat     = "${supp_dir}/picard_metric_files/refFlat.txt"

ref_Gene2Transcripts     = "${supp_dir}/ref.gene2transcripts.tsv"
ref_GbrsHybridTargets    = "${supp_dir}/gbrs.hybridized.targets.info"
ref_GenePosOrdered       = "${supp_dir}/ref.gene_pos.ordered.npz"
ref_GenePosOrdered_ver   = "${supp_dir}/ref.gene_pos.ordered_0.1.6.npz"

//~~~~~~~~~~ Initial Channel of SampleID and Fastq Data ~~~~
Channel.of( sampleID, fqR1path, fqR2path )
    .toList()
    .into { ch_sample_fastqs; ch_fastq_star }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Align fastq Bowtie ~~~~~
process align_fq {
    label 'bowtie'
    label 'high_mem'
    publishDir path:sample_outdir, mode:'move', pattern:"*.log"
    publishDir path:sample_outdir, mode:'copy', pattern:"*.fastq"
    
    input:
    tuple val(sampleID), path(fqR1), path(fqR2) from ch_sample_fastqs

    path(ref_transcripts)
    //path(ref_transcripts) from ref_transcripts

    output:
    tuple sampleID, path(align_1), path(align_2) into ch_aligned_sam
    path "${sampleID}*.log"

    script:
    align_1 = "align_1.sam"
    align_2 = "align_2.sam"
    """
    zcat ${fqR1} | bowtie \
      -p ${task.cpus} \
      -q -a --best --strata \
      --sam ${align_1} \
      -v 3 \
      ${ref_transcripts} \
      - \
      2>${sampleID}.bowtie_R1.log
    zcat ${fqR2} | bowtie \
      -p ${task.cpus} \
      -q -a --best --strata \
      --sam ${align_2} \
      -v 3 \
      ${ref_transcripts} \
      - \
      2>${sampleID}.bowtie_R2.log
    """
}
process aligns_to_bams {
    label 'samtools'

    input:
    tuple sampleID, path(sam_1), path(sam_2) from ch_aligned_sam

    output:
    tuple sampleID, path(bam_1), path(bam_2) into ch_bams

    script:
    """
    samtools view -bS ${sam_1} > ${bam_1}
    samtools view -bS ${sam_2} > ${bam_2}
    """
}


//~~~~~~~~~~~~~~~~~~~~~~ bam to emase, compress, merge ~~~~~
process bams_to_emase_to_merge {
    label 'gbrs'
    label 'high_mem'
    publishDir path:sample_outdir, mode:'move', pattern:"*.log"
    publishDir path:sample_outdir, mode:'copy', pattern:"*.merged_compressed.h5"
    //TODO: QUESTION: keep which files??

    input:
    tuple sampleID, path(align1), path(align2) from ch_bams
    path(ref) from params.gbrsSupport
    //path(ref_transcripts_info) //TODO: need ref file path as input?

    output:
    path "*.log"
    tuple sampleID, path(merged) into ch_aln_merged

    script:
    emase1 = ${sampleID}.emase1.h5
    emase2 = ${sampleID}.emase2.h5
    ecomp1 = ${sampleID}.compressed.emase1.h5
    ecomp1 = ${sampleID}.compressed.emase2.h5
    merged = ${sampleID}.merged_compressed.h5
    """
    gbrs bam2emase -i ${align1} \
                   -m ${ref}/ref.transcripts_info \
                   -s A,B,C,D,E,F,G,H \
                   -o ${emase1}
    gbrs compress -i ${emase1} -o ${ecomp1}

    gbrs bam2emase -i $align2 \
                   -m ${ref}/ref.transcripts_info \
                   -s A,B,C,D,E,F,G,H \
                   -o ${emase2}
    gbrs compress -i ${emase2} -o ${ecomp2}

    gbrs compress -i ${ecomp1},${ecomp2} -o ${merged}

    ls -l > gbrs_merge_dir.fof.log
    """
}


//~~~~~~~~~~~~~~~~  Quantification and Reconstruction  ~~~~~
process quantify {
    label 'gbrs'
    label 'high_mem'
    publishDir path:sample_outdir, mode:'move', pattern:"*.log"
    publishDir path:sample_outdir, mode:'move', pattern:"*.pdf"
    publishDir path:sample_outdir, mode:'copy' //, pattern:"*"

    input:
    tuple sampleID, path(merged) from ch_aln_merged
    path(ref) from params.gbrsSupport

    output:
    path "*.log"
    path "*.pdf"
    file "${sampleID}.gbrs.interpolated.genoprobs.npz" into ch_genoprobs
    tuple sampleID, file("${sampleID}.multiway.isoforms.tpm") into ch_genes_tpm

    script:
    model      = params.model
    sex        = params.sex
    generation = params.generation
    """
    gbrs quantify \
        -i ${merged} \
        -g ${ref_Gene2Transcripts} \
        -L ${ref_GbrsHybridTargets} \
        -M ${model} \
        --report-alignment-counts \
        -o ${sampleID}

    gbrs reconstruct \
        -e ${sampleID}.multiway.genes.tpm \
        -t ${ref_tranprob}/tranprob.DO.${generation}.${sex}.npz \
        -x ${ref_avecs} \
        -g ${ref_GenePosOrdered} \
        -o ${sampleID}

    gbrs quantify \
        -i ${merged} \
        -G ${sampleID}.genotypes.tsv \
        -g ${ref_Gene2Transcripts} \
        -L ${ref_GbrsHybridTargets} \
        -M ${model} \
        --report-alignment-counts \
        -o ${sampleID}

    gbrs interpolate \
        -i ${sampleID}.genoprobs.npz \
        -g ${ref_GenomeGrid} \
        -p ${ref_GenePosOrdered_ver} \
        -o ${sampleID}.gbrs.interpolated.genoprobs.npz

    gbrs plot \
        -i ${sampleID}.gbrs.interpolated.genoprobs.npz \
        -o ${sampleID}.gbrs.plotted.genome.pdf \
        -n ${sampleID}

    ls -l > quantify_dir.fof.log
    """

}

process exportgeno {
    label 'export'
    publishDir path:sample_outdir, mode:'move', pattern:"*.log"
    publishDir path:sample_outdir, mode:'copy', pattern:"*"

    input:
    path(ref_GenomeGrid)
    file(genoprob) from ch_genoprobs

    output:
    path "*.log"
    file "*" into export_out

    script:
    """
    export-genoprob-file
        -i ${genoprob} \
        -s A,B,C,D,E,F,G,H \
        -g ${ref_GenomeGrid}
    ls -l > export_dir.fof.log
    """
}


/* NOTE: if making star index is needed (remotely) adapt this section:
Channel
  .fromPath(params.gbrsSupport + "/C57BL6J/C57BL6J.fa")
  .into{ ch_star_genome; ch_rsem_genome }
Channel
  .fromPath(params.gbrsSupport + "/C57BL6J/C57BL6J.gtf")
  .into{ ch_star_gtf; ch_rsem_gtf }

// Generate STAR index
process makeSTARindex{
    label 'star'
    label 'high_mem'
    publishDir path:sample_outdir, mode:'copy', pattern:"*"

    input:
    file fasta from ch_star_genome
    file gtf from ch_star_gtf

    output:
    path "*.log"
    file "star/*" into ch_star_index

    script:
    def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    mkdir star
    STAR \\
      --runMode genomeGenerate \\
      --runThreadN ${task.cpus} \\
      --sjdbGTFfile $gtf \\
      --sjdbOverhang 149 \\
      --genomeDir star \\
      --genomeFastaFiles $fasta \\
      $avail_mem
    """
 }
 */

//~~ Quality metrics for RNA mapped against B6 with STAR ~~~
process star {
    label 'star'
    publishDir path:sample_outdir, mode:'move', pattern:"*.log"
    //publishDir path:sample_outdir, mode:'copy', pattern:"*"

    input:
    path(genomeDir) from ref_StarGenomeDir
    path(starIdx)   from ref_StarIdx
    //path(staridx) from star_index.collect()
    //tuple sampleID, path(fqR1), path(fqR2) from ch_fastq_star   
    tuple sampleID, fqR1, fqR2 from ch_fastq_star

    output:
    path "*.log"
    tuple sampleID, path("${sampleID}.Aligned.sortedByCoord.out.bam") into ch_star_bam


    script:
    """
    STAR \
        --runThreadN ${task.cpus} \
        --genomeDir ${genomeDir} \
        --sjdbGTFfile ${starIdx} \
        --outFileNamePrefix ${sampleID}. \
        --readFilesIn ${fqR1} ${fqR2} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --quantMode GeneCounts \
        --outFilterMultimapNmax 5
    ls -l > star_dir.fof.log
    """
}


process picard {
    label 'picard'
    label 'high_mem'
    publishDir path:sample_outdir, mode:'move', pattern:"*.log"

    input:
    tuple sampleID, path("${sampleID}.Aligned.sortedByCoord.out.bam") from ch_star_bam
    path(ref_StarIdxGenome)
    path(ref_PicardMetricsFlat)

    output:
    path "*.log"
    tuple sampleID, path("${sampleID}.alignment.read.group.reorder.sorted.bam") into ch_picard_bam
    tuple sampleID, path("${sampleID}.picard_aln_metrics.txt") into ch_picard_log

    script:
    """
    picard AddOrReplaceReadGroups \
        INPUT=${sampleID}.Aligned.sortedByCoord.out.bam \
        OUTPUT=${sampleID}.alignment.read.group.bam \
        SORT_ORDER=coordinate CREATE_INDEX=true \
        RGID=ID RGLB=LIB RGPL=ILLUMINA RGSM=SAMPLE RGPU=Lane

    picard ReorderSam \
        INPUT=${sampleID}.alignment.read.group.bam \
        OUTPUT=${sampleID}.alignment.read.group.reorder.bam \
        REFERENCE=${ref_StarIdxGenome} \
        CREATE_INDEX=true

    picard SortSam \
        SO=coordinate  \
        INPUT=${sampleID}.alignment.read.group.reorder.bam \
        OUTPUT=${sampleID}.alignment.read.group.reorder.sorted.bam \
        VALIDATION_STRINGENCY=SILENT \
        CREATE_INDEX=true

    picard CollectRnaSeqMetrics \
        I=${sampleID}.alignment.read.group.reorder.sorted.bam \
        O=${sampleID}.picard_aln_metrics.txt \
        REF_FLAT=${ref_PicardMetricsFlat} \
        RIBOSOMAL_INTERVALS=${ref_PicardMetrics}
        STRAND=${LIBRARY}
    ls -l > picard_dir.fof.log
    """
}


process qc_metrics {
    label 'perl'
    publishDir path:sample_outdir, mode:'move', pattern:"*.log"
    publishDir path:sample_outdir, mode:'move', pattern:"*_stats.txt"
    publishDir path:sample_outdir, mode:'move', pattern:"*_metrics.txt"
    //FIXME: which process does "${sampleID}.Log.final.out" come from??

    input:
    tuple sampleID, path(picard_aln_metrics) from ch_picard_log

    output:
    path "*.log"
    path{"${sampleID}.summary_stats.txt"}

    script:
    log_final = "${sampleID}.Log.final.out"
    """
    perl ${projectDir}/bin/summary_QC_metrics_starAligner.pl \
        ${log_final} \
        ${picard_aln_metrics} \
        > ${sampleID}.summary_stats.txt
    """
}



// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Closing Summary ~ ~ ~ ~ ~
workflow.onComplete {
    wfEnd = [:]
    wfEnd['Completed at'] = workflow.complete
    wfEnd['Duration']     = workflow.duration
    wfEnd['Exit status']  = workflow.exitStatus
    wfEnd['Success']      = workflow.success
    if(!workflow.success){
        wfEnd['!!Execution failed'] = ''
        wfEnd['.    Error']   = workflow.errorMessage
        wfEnd['.    Report']  = workflow.errorReport
    }
    Summary.show(wfEnd)
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// vim: set ft=groovy.nextflow ts=4 sw=0 tw=100 et fdm=syntax:
