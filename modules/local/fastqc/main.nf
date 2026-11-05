nextflow.enable.types = true

/*
 * Raw read quality control with FastQC.
 */
process FASTQC {
    tag id
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/816bda6dd014e810753f642611232443ddeabb2b47f9681ad43f9044816f00ff/data'
        : 'community.wave.seqera.io/library/fastqc:0.13.0--b342db7b694c9af4' }"

    input:
    record(id: String, reads: List<Path>)

    output:
    record(id: id, html: files('*.html'), zip: files('*.zip'))

    topic:
    tuple(task.process, 'fastqc', eval('fastqc --version | sed "/FastQC v/!d; s/.*v//"')) >> 'versions'

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: id
    // FastQC allocates --memory per thread; keep within the task allocation (100-10000 MB)
    def mem_mb = task.memory ? (task.memory.toUnit('MB') / task.cpus).intValue() : 0
    def memory = mem_mb > 0 ? "--memory ${Math.max(100, Math.min(10000, mem_mb))}" : ''
    """
    # Rename inputs to <sample>[_<n>].fastq.gz for consistent naming in MultiQC
    n=0
    files=()
    for f in ${reads.join(' ')}; do
        n=\$((n + 1))
        if [ ${reads.size()} -eq 1 ]; then new="${prefix}.fastq.gz"; else new="${prefix}_\${n}.fastq.gz"; fi
        [ -e "\$new" ] || ln -s "\$f" "\$new"
        files+=("\$new")
    done

    fastqc \\
        ${args} \\
        --threads ${task.cpus} \\
        ${memory} \\
        "\${files[@]}"
    """

    stub:
    def prefix = task.ext.prefix ?: id
    """
    touch ${prefix}_fastqc.html ${prefix}_fastqc.zip
    """
}
