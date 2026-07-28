process FASTP {

    // container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    publishDir "results/fastp", mode: 'copy'

    tag "${sample_id}"

    cpus params.fastp_cpus

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
    # Run FASTP for quality filtering and adapter trimming
    fastp \
        -i ${read} \
        -o ${sample_id}.trimmed.fastq.gz \
        --adapter_sequence ${params.fastp_adapter} \
        --trim_poly_g \
        --qualified_quality_phred ${params.min_quality} \
        --unqualified_percent_limit ${params.unqualified_percent} \
        --n_base_limit ${params.n_base_limit} \
        --length_required ${params.min_length} \
        --thread ${task.cpus} \
        --html ${sample_id}.fastp.html \
        --json ${sample_id}.fastp.json
    """
}