nextflow.enable.types = true

/*
 * Relative abundance of the reference genomes (one FASTA per species) with CoverM.
 */
process COVERM_GENOME {
    tag id
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c4/c4402aed2ff8bd062a013369c9049e68bdd159e849b271a227a8ce2727e941e1/data'
        : 'community.wave.seqera.io/library/coverm:0.7.0--c57cbb9db4d660fb' }"

    input:
    record(id: String, single_end: Boolean, reads: List<Path>)
    genomes: Path
    extension: String

    output:
    record(id: id, abundance: file("${task.ext.prefix ?: id}.coverm.tsv"))

    topic:
    tuple(task.process, 'coverm', eval("coverm --version | sed 's/coverm //'")) >> 'versions'

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: id
    def in_reads = single_end ? "--single ${reads[0]}" : "--coupled ${reads[0]} ${reads[1]}"
    """
    TMPDIR=.

    coverm genome \\
        --threads ${task.cpus} \\
        ${in_reads} \\
        --genome-fasta-directory ${genomes} \\
        --genome-fasta-extension ${extension} \\
        ${args} \\
        --output-file ${prefix}.coverm.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: id
    """
    printf 'Genome\\t${prefix} Relative Abundance (%%)\\nunmapped\\t100\\n' > ${prefix}.coverm.tsv
    """
}
