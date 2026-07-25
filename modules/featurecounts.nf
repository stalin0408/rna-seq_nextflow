process FEATURECOUNTS {
    publishDir "results/featurecounts/", mode: 'copy'

    input:
    path bams
    path gtf

    output:
    path "gene_counts.txt", emit: counts

    script:
script:
"""
echo " BAM files received:"
printf "%s\n" ${bams}
featureCounts \
    -T 8 \
    -a ${gtf} \
    -F GTF \
    -t exon \
    -g gene_id \
    -o gene_counts.txt \
    ${bams}
"""
}
