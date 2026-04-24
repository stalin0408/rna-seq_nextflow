#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = "data/*.fastq.gz"

// Input channel
Channel
    .fromPath(params.input)
    .map { file -> tuple(file.baseName, file) }
    .set { reads_ch }


// Detect read type
process DETECT_READ_TYPE {

    input:
    path read

    output:
    tuple val(type), path(read)

    script:
    """
    type=\$(python3 bin/detect_read_type.py ${read})
    echo \$type
    """
}


// Main workflow
workflow {

    detected = reads_ch | DETECT_READ_TYPE

    detected
        .branch {
            long: it[0] == "long"
            short: it[0] == "short"
        }
        .set { out_ch }

    workflow_long(out_ch.long)
    workflow_short(out_ch.short)
}


// Import sub-workflows
include { workflow_long } from './workflows/long_read.nf'
include { workflow_short } from './workflows/short_read.nf'