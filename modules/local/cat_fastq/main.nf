nextflow.enable.types = true

/*
 * Concatenate FASTQ files of samples sequenced over multiple runs/lanes.
 */
process CAT_FASTQ {
    tag id
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8' }"

    input:
    record(id: String, single_end: Boolean, reads_1: List<Path>, reads_2: List<Path>)

    stage:
    stageAs reads_1, 'r1_*/*'
    stageAs reads_2, 'r2_*/*'

    output:
    record(id: id, single_end: single_end, reads: files('*.merged.fastq.gz'))

    topic:
    tuple(task.process, 'cat', eval("cat --version 2>&1 | head -n 1 | sed 's/^.*coreutils) //; s/ .*\$//'")) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: id
    def unzip  = reads_1[0].name.endsWith('.gz') ? '' : '| gzip'
    if (single_end) {
        """
        cat ${reads_1.join(' ')} ${unzip} > ${prefix}.merged.fastq.gz
        """
    } else {
        """
        cat ${reads_1.join(' ')} ${unzip} > ${prefix}_1.merged.fastq.gz
        cat ${reads_2.join(' ')} ${unzip} > ${prefix}_2.merged.fastq.gz
        """
    }

    stub:
    def prefix = task.ext.prefix ?: id
    def files = single_end ? "${prefix}.merged.fastq.gz" : "${prefix}_1.merged.fastq.gz ${prefix}_2.merged.fastq.gz"
    """
    for f in ${files}; do echo '' | gzip > \$f; done
    """
}
