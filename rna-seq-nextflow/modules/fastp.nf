process FASTP {

    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    publishDir "results/fastp", mode: 'copy'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id),
          path("${sample_id}.trimmed.fastq.gz"),
          emit: trimmed_reads

    path "${sample_id}.fastp.html"
    path "${sample_id}.fastp.json"

    script:
    """
    fastp \
        -i ${read} \
        -o ${sample_id}.trimmed.fastq.gz \
        --adapter_sequence auto\
        --cut_tail \
        --cut_mean_quality 20 \
        --length_required 30 \
        --thread 4 \
        --html ${sample_id}.fastp.html \
        --json ${sample_id}.fastp.json
    """
}
