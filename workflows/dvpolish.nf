/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowDvpolish.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DVPOLISH_CHUNKFA             } from "$projectDir/modules/local/dvpolish/chunkfa"
include { DVPOLISH_PBMM2_INDEX         } from "$projectDir/modules/local/dvpolish/pbmm2_index"
include { DVPOLISH_PBMM2_ALIGN         } from "$projectDir/modules/local/dvpolish/pbmm2_align"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { SAMTOOLS_FAIDX              } from "$projectDir/modules/nf-core/samtools/faidx/main"
include { SAMTOOLS_VIEW               } from "$projectDir/modules/nf-core/samtools/view/main"
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTER } from "$projectDir/modules/nf-core/samtools/index/main"
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MERGE  } from "$projectDir/modules/nf-core/samtools/index/main"
include { SAMTOOLS_MERGE              } from "$projectDir/modules/nf-core/samtools/merge/main"
include { DEEPVARIANT                 } from "$projectDir/modules/nf-core/deepvariant/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

def meta_map = [
    id: "Species_name",
    single_end: true
]

workflow DVPOLISH {

    //ch_versions = Channel.empty()

    Channel.fromPath(params.fasta_file, checkIfExists: true)
        .map{asm -> [ meta_map, asm ]}
        .collect()          // to make a value channel, otherwise it will only be used once in the alignmnent step and than its consumed !!!
        .set{asm_file}

    Channel.fromPath(params.reads_file, checkIfExists: true)
                .map{reads -> [meta_map, reads ]}
        .set{reads_file}

    //
    // MODULE: Run SAMTOOLS_FAIDX
    //
    SAMTOOLS_FAIDX (
        asm_file,
        [meta_map, []]
    )

    //
    // MODULE: Run SPLIT_FA
    //
    DVPOLISH_CHUNKFA (
        SAMTOOLS_FAIDX.out.fai
    )

    //
    // MODULE: PBMM2_INDEX
    //
    DVPOLISH_PBMM2_INDEX (
        asm_file
    )


    reads_file.view{it: println("reads_files: (plain)" + it) }
    reads_file.groupTuple().view{it: println("reads_files: (tuple)" + it) }

    //
    // MODULE: PBMM2_ALIGN
    //
    DVPOLISH_PBMM2_ALIGN (
        reads_file,
        asm_file
    )

    //DVPOLISH_PBMM2_ALIGN.out.bam.view {it: println("DVPOLISH_PBMM2_ALIGN (bam): " + it)}



def path_closure = {meta, files -> files.collect(){[meta, it ]}}

DVPOLISH_PBMM2_ALIGN.out.bam
    .flatMap(path_closure)
    .combine(DVPOLISH_CHUNKFA.out.bed.flatMap(path_closure), by:0)
    .map { map, bam, bed -> [ map + [mergeID: bed.baseName], bam, bed]

    }
    .set { bam_bed_ch }

//    bam_bed_ch.view()


DVPOLISH_PBMM2_ALIGN.out.bai
    .flatMap(path_closure)
    .combine(DVPOLISH_CHUNKFA.out.bed.flatMap(path_closure), by:0)
    .map { meta, bai, bed -> bai }
    .set { bai_ch }

    //bai_ch.view()

    //
    // MODULE: Run SAMTOOLS_VIEW
    //
    SAMTOOLS_VIEW (bam_bed_ch, [[],[]], bai_ch)

    //
    // MODULE: Run SAMTOOLS_INDEX
    //
    SAMTOOLS_INDEX_FILTER(SAMTOOLS_VIEW.out.bam)

//    SAMTOOLS_VIEW.out.bam.view()

//    SAMTOOLS_VIEW.out.bam.groupTuple(by:1).view()

SAMTOOLS_VIEW.out.bam
    .groupTuple(by:0)
    .branch { meta, bam_list ->
        merge: bam_list.size() > 1
        link: true
    }
    .set { bam_merge_ch }

// bam_merge_ch.merge.view{ println("falk merge: " + it + " size: " + it.size() + ", " + it[0].size() + ", " + it[1].size())}
// bam_merge_ch.do_nothing.view{ println("falk do nothing: " + it)}

    SAMTOOLS_MERGE(
        bam_merge_ch.merge,
        [[],[]],
        [[],[]]
    )
    SAMTOOLS_INDEX_MERGE(SAMTOOLS_MERGE.out.bam)

    //bam_merge_ch.link.view()
    //SAMTOOLS_INDEX_FILTER.out.bai.view()

    bam_merge_ch.link
    .map { meta, bam -> [ meta, *bam ]} // the spread operator (*) flattens the bam lsit
    .join(SAMTOOLS_INDEX_FILTER.out.bai, by:0)
    .mix(SAMTOOLS_MERGE.out.bam
        .join(SAMTOOLS_INDEX_MERGE.out.bai, by:0)
    )
    .join(bam_bed_ch
    .map { meta, bam, bed -> [meta, bed]}
    .unique())
    .set {deepvariant_ch}


    deepvariant_ch.view()

//    bam_bed_ch
//    .map { meta, bam, bed -> [meta, bed]}
//    .unique().view()
//
    // run deepVariant
    DEEPVARIANT(
        deepvariant_ch,
        asm_file,
        SAMTOOLS_FAIDX.out.fai,
        [[],[]]     // tuple val(meta4), path(gzi)
    )



    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    //INPUT_CHECK (
        // file(params.input)
    // )
    // ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //
    //FASTQC (
        // INPUT_CHECK.out.reads
    // )

    //
    // MODULE: Run CAT FASTQC
    //
    // CAT_FASTQ (
        // INPUT_CHECK.out.reads
    // )
//
//
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())
//
    // CUSTOM_DUMPSOFTWAREVERSIONS (
        // ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    //
    // MODULE: MultiQC
    //
    // workflow_summary    = WorkflowDvpolish.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)
//
    // methods_description    = WorkflowDvpolish.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    // ch_methods_description = Channel.value(methods_description)
//
    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
//
    // MULTIQC (
        // ch_multiqc_files.collect(),
        // ch_multiqc_config.toList(),
        // ch_multiqc_custom_config.toList(),
        // ch_multiqc_logo.toList()
    // )
    // multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
