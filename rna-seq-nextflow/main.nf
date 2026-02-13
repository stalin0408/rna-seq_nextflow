nextflow.enable.dsl=2

include { FASTQC } from './modules/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/fastqc'
include { FASTP } from './modules/fastp'

workflow {
    reads = Channel
        .fromPath("data/*.fastq.gz")
        .map{
            file ->
            tuple(file.baseName, file)
        }
    qc_raw = FASTQC(reads)
    trimmed = FASTP(reads)
    qc_trimmed = FASTQC_TRIM(trimmed.trimmed_reads)
}