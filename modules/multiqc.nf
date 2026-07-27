process MULTIQC{
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path results_dir

    output:
    path "multiqc_report.html"
    script:
    """
    multiqc ${results_dir} -o .
    """

}