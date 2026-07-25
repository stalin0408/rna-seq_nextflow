//import module
include { DESEQ2_ANALYSIS } from '../modules/deseq2.nf'
//include { GENERATE_REPORT } from '../modules/report.nf'


workflow RNASEQ_WORKFLOW {

    take:
    counts
    metadata
    script_file

    main:
    results = DESEQ2_ANALYSIS(
        counts,
        metadata,
        script_file
    )

//    GENERATE_REPORT(results.collect())

    emit:
    results
//    final_report = GENERATE_REPORT.out
}

