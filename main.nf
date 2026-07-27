nextflow.enable.dsl=2

include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/fastp'
include { STAR } from './modules/align'
include { FEATURECOUNTS } from './modules/featurecounts'
include { MULTIQC } from './modules/multiqc'

//include workflow
include { RNASEQ_WORKFLOW } from './workflows/rnaseq_workflow'

//setting channel
workflow {
    reads = Channel
        .fromPath("${params.input_dir}/*.{fastq,fq,fastq.gz,fq.gz}")
        .map{
            file ->
            tuple(file.baseName, file)
        }
    reads.view()


    println "INPUT_DIR = ${params.input_dir}"
    qc_raw = FASTQC(reads)

    trimmed = FASTP(reads)
    
    aligned_reads = STAR(trimmed.trimmed_reads)
    bam_files = aligned_reads.aligned_reads
        .map{ sample_id, bam -> bam }
        .collect()
    featurecounts_out = FEATURECOUNTS(bam_files, file(params.gtf))
    deseq2_result = RNASEQ_WORKFLOW(
        featurecounts_out.counts,
        file(params.meta),
        file(params.script_file),
        params.design)
    MULTIQC(file(params.outdir))
   

}
