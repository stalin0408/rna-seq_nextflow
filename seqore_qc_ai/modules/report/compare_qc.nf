process COMPARE_QC {

    input:
    tuple val(sample_id), path(before), path(after)

    output:
    tuple val(sample_id), path("comparison.json")

    script:
    """
    python3 bin/compare_nanoplot.py ${before} ${after} > comparison.json
    """
}