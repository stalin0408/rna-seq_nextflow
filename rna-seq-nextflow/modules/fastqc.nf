process FASTQC{
    container 'biocontainers/fastqc:v0.11.9_cv8'
    publishDir "results/fastqc", mode: 'copy'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id), path("*.html")

    script:
    """
    fastqc $read
    """
}