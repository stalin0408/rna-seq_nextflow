process NANO_PLOT_POST {

    tag "${reads.simpleName}"

    publishDir "results/post_qc", mode: 'copy'

    input:
    path reads

    output:
    path "nanoplot_post/*"

    script:
    """
    mkdir -p nanoplot_post
    NanoPlot --fastq ${reads} -o nanoplot_post
    """
}