process PARSE_NANOPLOT {

    input:
    tuple val(sample_id), path(summary_file)

    output:
    tuple val(sample_id), path("metrics.json")

    script:
    """
    python3 bin/parse_nanoplot.py ${summary_file} > metrics.json
    """
}