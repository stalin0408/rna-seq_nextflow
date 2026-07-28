nextflow.enable.dsl = 2

// include modules
include { FASTQC }        from './modules/fastqc'
include { FASTP }         from './modules/fastp'
include { STAR_INDEX }    from './modules/star_index'
include { STAR }          from './modules/align'
include { FEATURECOUNTS } from './modules/featurecounts'
include { MULTIQC }       from './modules/multiqc'

// include workflow
include { RNASEQ_WORKFLOW } from './workflows/rnaseq_workflow'

// setting channel
workflow {

    reads = Channel
        .fromPath("${params.input_dir}/*.{fastq,fq,fastq.gz,fq.gz}")
        .map { file ->
            tuple(file.baseName, file)
        }

    FASTQC(reads)

    trimmed = FASTP(reads)

    //Build STAR index if requested
    
    if (params.star_index) {

        star_index = Channel.value(file(params.star_index))

    }
    else if (params.genome_fasta && params.annotation) {

        STAR_INDEX(
            Channel.value(file(params.genome_fasta)),
            Channel.value(file(params.annotation))
        )

        star_index = STAR_INDEX.out.data_index

    }
    else {

        error """
        Please provide either:

        --star_index <STAR index directory>

        OR

        --genome_fasta <genome.fa> \
        --annotation <genes.gtf>
        """
    }

    star_input = trimmed.trimmed_reads
        .combine(star_index)
        .map { sample_id, read, genome_index ->
            tuple(sample_id, read, genome_index)
        }

    STAR(star_input)

    bam_files = STAR.out.aligned_reads
        .map { sample_id, bam -> bam }
        .collect()

    FEATURECOUNTS(
        bam_files,
        file(params.annotation)
    )

    RNASEQ_WORKFLOW(
        FEATURECOUNTS.out.counts,
        file(params.meta),
        file(params.script_file),
        params.design
    )

    MULTIQC(file(params.outdir))

}