process FASTQC{
        
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