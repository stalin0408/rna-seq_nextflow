process STAR{

    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_0'

    publishDir "results/star", mode: 'copy'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id),
          path("${sample_id}.Aligned.sortedByCoord.out.bam"),
          emit: aligned_reads

    path("${sample_id}.Log.final.out")
    path("${sample_id}.Log.out")
    path("${sample_id}.Log.progress.out")

    script:
    """
    STAR --runThreadN ${task.cpus} \
        --genomeDir ${params.genome} \
        --readFilesIn  ${read} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sample_id}. \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --outSAMattrRGline ID:${sample_id} SM:${sample_id} \
        --quantMode GeneCounts
    """
}
