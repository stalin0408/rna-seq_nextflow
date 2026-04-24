process NANO_PLOT {

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("nanoplot")

    script:
    """
    mkdir -p nanoplot
    NanoPlot --fastq ${reads} -o nanoplot
    """
}