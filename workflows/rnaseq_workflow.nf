//import module
include { DESEQ2_ANALYSIS } from '../modules/deseq2.nf'
//include { GENERATE_REPORT } from '../modules/report.nf'


workflow RNASEQ_WORKFLOW {

    take:
    counts
    metadata
    script_file
    design

    main:

    println ">>> Entered RNASEQ_WORKFLOW"

    counts.view { "Counts file = $it" }

    println ">>> Calling DESEQ2_ANALYSIS"
    
    results = DESEQ2_ANALYSIS(
        counts,
        metadata,
        script_file,
        design
    )
    println ">>> Returned from DESEQ2_ANALYSIS"


//    GENERATE_REPORT(results.collect())

    emit:
    results
//    final_report = GENERATE_REPORT.out
}

