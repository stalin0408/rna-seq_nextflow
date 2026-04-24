nextflow.enable.dsl=2

// -----------------------------
// IMPORT MODULES
// -----------------------------
include { NANO_PLOT }        from '../modules/qc/nanoplot.nf'
include { FILTLONG }         from '../modules/filtering/filtlong.nf'
include { NANO_PLOT_POST }   from '../modules/post_qc/nanoplot_post.nf'
// Optional (when you add parsing):
 include { PARSE_NANOPLOT } from '../modules/qc/parse_nanoplot.nf'


// -----------------------------
// WORKFLOW: LONG READ PIPELINE
// -----------------------------
workflow workflow_long {

    take:
    reads_ch   // tuple: [ "long", file ]

    main:

    // Extract only FASTQ file from tuple
    reads = reads_ch.map { it[1] }

    // -------------------------
    // STEP 1: QC (NanoPlot)
    // -------------------------
    qc1 = reads | NANO_PLOT

    // -------------------------
    // (OPTIONAL) STEP 1b: Parse QC
    // -------------------------
    // Extract summary.txt from NanoPlot output
    summary = qc1.map { sample_id, files ->
    def f = files.find { it.name == 'summary.txt' }
    tuple(sample_id, f)
}
    metrics = summary | PARSE_NANOPLOT

    if metrics["reads_above_Q15"]["percent"] < 40:
    print("Low quality dataset → increase filtering stringency")

    // -------------------------
    // STEP 2: FILTERING (Filtlong)
    // -------------------------
    filtered = reads | FILTLONG

    // -------------------------
    // STEP 3: POST-QC
    // -------------------------
    qc2 = filtered | NANO_PLOT_POST

    emit:
    qc_initial = qc1
    filtered_reads = filtered
    qc_post = qc2
    metrics_initial = metrics   // enable when parser is used

    // Initial QC
    qc1 = reads | NANO_PLOT
    summary1 = qc1.map { it.find { f -> f.name == 'summary.txt' } }
    metrics_initial = summary1 | PARSE_NANOPLOT

    // Filtering (with params already added earlier)
    filtered = reads.combine(params_ch) | FILTLONG

    // Post QC
    qc2 = filtered | NANO_PLOT_POST
    summary2 = qc2.map { it.find { f -> f.name == 'summary.txt' } }
    metrics_post = summary2 | PARSE_NANOPLOT

    comparison_input = metrics_initial.join(metrics_post)

    comparison = comparison_input | COMPARE_QC

    emit:
    qc_initial = qc1
    filtered_reads = filtered
    qc_post = qc2
    metrics_initial = metrics_initial
    metrics_post = metrics_post
    comparison_report = comparison
}