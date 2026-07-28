process STAR {

    publishDir "results/star_align", mode: 'copy'

    tag "${sample_id}"

    cpus params.star_cpus

    input:
    tuple val(sample_id), path(read), path(genome_index)

    output:
    tuple val(sample_id),
          path("${sample_id}.Aligned.sortedByCoord.out.bam"),
          emit: aligned_reads

    tuple val(sample_id),
          path("${sample_id}.ReadsPerGene.out.tab"),
          emit: gene_counts

    path("${sample_id}.Log.final.out"),    emit: final_log
    path("${sample_id}.Log.out"),          emit: star_log
    path("${sample_id}.Log.progress.out"), emit: progress_log

    script:
    """
    STAR --runThreadN ${task.cpus} \
         --genomeDir ${genome_index} \
         --readFilesIn ${read} \
         --readFilesCommand zcat \
         --outFileNamePrefix ${sample_id}. \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMstrandField intronMotif \
         --outSAMattrRGline ID:${sample_id} SM:${sample_id} \
         --quantMode GeneCounts
    """
}