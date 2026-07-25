process DESEQ2_ANALYSIS {

    publishDir "${params.outdir}/", mode: 'copy'

    conda "rnaseq"

    input:
    path counts
    path metadata
    path script_file

    output:
    path "deseq2_results"

    script:
    """
    mkdir -p deseq2_results

    Rscript ${script_file} \
        ${counts} \
        ${metadata} \
        deseq2_results
    """
}