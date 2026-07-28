process STAR_INDEX {

    tag "STAR Genome Index"

    publishDir "${params.outdir}/star_index", mode: 'copy'

    cpus params.star_index_cpus

    input:
    path genome
    path gtf

    output:
    path "star_index", emit: data_index

    script:
    """
    # Create the output directory for the STAR genome index
    mkdir -p star_index

    # Build the STAR genome index
    STAR \
        --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir star_index \
        --genomeFastaFiles ${genome} \
        --sjdbGTFfile ${gtf} \
        --sjdbOverhang ${params.sjdbOverhang}
    """
}