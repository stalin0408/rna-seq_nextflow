process FEATURECOUNTS {
    
    publishDir "results/featurecounts/", mode: 'copy'

    cpus params.featureCounts_cpus

    input:
    path bams
    path gtf

    output:
    path "gene_counts.txt", emit: counts

    script:
    """
    echo " BAM files received:"
    printf "%s\n" ${bams}
    featureCounts \
        -T ${task.cpus} \
        -a ${gtf} \
        -F GTF \
        -t exon \
        -g gene_id \
        -o gene_counts.txt \
        ${bams}
    """
}
