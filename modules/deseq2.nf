process DESEQ2_ANALYSIS {

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    path counts
    path metadata
    path script_file
    val design

    output:
    path "deseq2_results"

    script:
    """
    echo "===== DEBUG ====="
    mkdir -p deseq2_results
    which Rscript
    Rscript -e "print(.libPaths())"
    Rscript ${script_file} \
        ${counts} \
        ${metadata} \
        "${design}" \
        deseq2_results
    """
}