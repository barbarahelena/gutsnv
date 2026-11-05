nextflow.enable.types = true

/*
 * Read counts of the reads that enter genotyping (sequencing-depth covariate).
 */
process SEQKIT_STATS {
    tag id
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/29/297e985dd395d87f9089b314dbe685ddce59fa0770c4d1f5e1009bc7b6ec7a8d/data'
        : 'community.wave.seqera.io/library/seqkit:2.13.0--3bacfb68370dd835' }"

    input:
    record(id: String, reads: List<Path>)

    output:
    record(id: id, stats: file("${task.ext.prefix ?: id}.seqkit_stats.tsv"))

    topic:
    tuple(task.process, 'seqkit', eval("seqkit version | sed 's/^.*v//'")) >> 'versions'

    script:
    def args   = task.ext.args ?: '--all'
    def prefix = task.ext.prefix ?: id
    """
    seqkit stats \\
        --tabular \\
        --threads ${task.cpus} \\
        ${args} \\
        ${reads.join(' ')} > ${prefix}.seqkit_stats.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: id
    """
    printf 'file\\tformat\\ttype\\tnum_seqs\\tsum_len\\n' > ${prefix}.seqkit_stats.tsv
    for f in ${reads.join(' ')}; do printf '%s\\tFASTQ\\tDNA\\t1000\\t150000\\n' \$f >> ${prefix}.seqkit_stats.tsv; done
    """
}
