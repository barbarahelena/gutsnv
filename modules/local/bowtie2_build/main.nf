nextflow.enable.types = true

/*
 * Build a Bowtie2 index of the host genome.
 */
process BOWTIE2_BUILD {
    tag fasta.name
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d4990a9e13c1aad32d60f19e0ced204f2f1e231fe856290ce75713505a1bfe44/data'
        : 'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:dbd18842ad7cb278' }"

    input:
    fasta: Path

    output:
    file('bowtie2')

    topic:
    tuple(task.process, 'bowtie2', eval("bowtie2 --version 2>&1 | sed -n 's/.*bowtie2-align-s version //p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bowtie2
    bowtie2-build ${args} --threads ${task.cpus} ${fasta} bowtie2/${fasta.baseName}
    """

    stub:
    """
    mkdir bowtie2
    touch bowtie2/${fasta.baseName}.{1,2,3,4,rev.1,rev.2}.bt2
    """
}
