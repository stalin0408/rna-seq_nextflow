process FILTLONG {

    input:
    tuple path(reads), val(params)

    output:
    path "filtered.fastq.gz"

    script:
    """
    filtlong --min_length ${params.min_len} --min_mean_q ${params.min_q} ${reads} | gzip > filtered.fastq.gz
    """
}