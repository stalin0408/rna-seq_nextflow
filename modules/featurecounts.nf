process FEATURECOUNTS {

    publishDir "results/featurecounts/", mode: 'copy'

    cpus params.featureCounts_cpus

    input:
    path bams
    path genome_index

    output:
    path "gene_counts.txt",         emit: counts
    path "gene_counts.txt.summary"

    script:
    """
    # Display the BAM files received for quantification
    echo " BAM files received:"
    printf "%s\n" ${bams}

    # Generate gene-level read counts using featureCounts
    featureCounts \
        -T ${task.cpus} \
        -a ${genome_index} \
        -F GTF \
        -t exon \
        -g gene_id \
        -o gene_counts.txt \
        ${bams}
    """
}