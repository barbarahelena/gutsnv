nextflow.enable.types = true

/*
 * Remove host (human) reads by aligning to the host genome with Bowtie2 and keeping
 * the read pairs that do not align concordantly. Alignments themselves are discarded.
 */
process BOWTIE2_HOSTREMOVAL {
    tag id
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d4990a9e13c1aad32d60f19e0ced204f2f1e231fe856290ce75713505a1bfe44/data'
        : 'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:dbd18842ad7cb278' }"

    input:
    record(id: String, single_end: Boolean, reads: List<Path>)
    index: Path

    output:
    record(
        id: id,
        single_end: single_end,
        reads: files('*.host_removed*.fastq.gz'),
        log: file("${task.ext.prefix ?: id}.bowtie2.log")
    )

    topic:
    tuple(task.process, 'bowtie2', eval("bowtie2 --version 2>&1 | sed -n 's/.*bowtie2-align-s version //p'")) >> 'versions'

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: id
    def in_reads = single_end ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def unmapped = single_end ? "--un-gz ${prefix}.host_removed.fastq.gz" : "--un-conc-gz ${prefix}.host_removed_%.fastq.gz"
    """
    INDEX=\$(find -L ${index}/ -name "*.rev.1.bt2" | sed 's/\\.rev\\.1\\.bt2\$//')
    [ -z "\$INDEX" ] && INDEX=\$(find -L ${index}/ -name "*.rev.1.bt2l" | sed 's/\\.rev\\.1\\.bt2l\$//')
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found in ${index}" >&2 && exit 1

    bowtie2 \\
        -x \$INDEX \\
        ${in_reads} \\
        --threads ${task.cpus} \\
        ${unmapped} \\
        ${args} \\
        -S /dev/null \\
        2>| >(tee ${prefix}.bowtie2.log >&2)
    """

    stub:
    def prefix = task.ext.prefix ?: id
    def files  = single_end ? "${prefix}.host_removed.fastq.gz" : "${prefix}.host_removed_1.fastq.gz ${prefix}.host_removed_2.fastq.gz"
    """
    for f in ${files}; do echo '' | gzip > \$f; done
    touch ${prefix}.bowtie2.log
    """
}
