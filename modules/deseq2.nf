process DESEQ2_ANALYSIS {

    publishDir "${params.outdir}/", mode: 'copy'

    input:
    path counts
    path metadata
    path script_file
    val  design

    output:
    path "deseq2_results"

    script:
    """
    # Display execution status
    echo "=== STARTING DESEQ2 ==="

    # Verify R installation
    which Rscript

    # Create output directory
    mkdir -p deseq2_results

    # Display R executable and library paths
    which Rscript
    Rscript -e "print(.libPaths())"

    # Run DESeq2 analysis
    Rscript ${script_file} \
        ${counts} \
        ${metadata} \
        "${design}" \
        deseq2_results
    """
}