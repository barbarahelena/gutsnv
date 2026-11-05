nextflow.enable.types = true

/*
 * Adapter and quality trimming with fastp.
 */
process FASTP {
    tag id
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/77/779b25fffc4544e4d4e3bebf1c1a1f2fee9166515bfa15908b2fee95f513a402/data'
        : 'community.wave.seqera.io/library/fastp:1.4.0--df415eab5e33cb11' }"

    input:
    record(id: String, single_end: Boolean, reads: List<Path>)

    output:
    record(
        id: id,
        single_end: single_end,
        reads: files('*.fastp.fastq.gz'),
        json: file("${task.ext.prefix ?: id}.fastp.json"),
        html: file("${task.ext.prefix ?: id}.fastp.html"),
        log: file("${task.ext.prefix ?: id}.fastp.log")
    )

    topic:
    tuple(task.process, 'fastp', eval('fastp --version 2>&1 | sed -e "s/fastp //g"')) >> 'versions'

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: id
    // fastp --thread excludes its reader threads; keep >= 2 worker threads for gzip compression
    def threads = Math.max(2, task.cpus - (single_end ? 1 : 2))
    if (single_end) {
        """
        [ -e ${prefix}.fastq.gz ] || ln -s ${reads[0]} ${prefix}.fastq.gz

        fastp \\
            --in1 ${prefix}.fastq.gz \\
            --out1 ${prefix}.fastp.fastq.gz \\
            --thread ${threads} \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            ${args} \\
            2>| >(tee ${prefix}.fastp.log >&2)
        """
    } else {
        """
        [ -e ${prefix}_R1.fastq.gz ] || ln -s ${reads[0]} ${prefix}_R1.fastq.gz
        [ -e ${prefix}_R2.fastq.gz ] || ln -s ${reads[1]} ${prefix}_R2.fastq.gz

        fastp \\
            --in1 ${prefix}_R1.fastq.gz \\
            --in2 ${prefix}_R2.fastq.gz \\
            --out1 ${prefix}_R1.fastp.fastq.gz \\
            --out2 ${prefix}_R2.fastp.fastq.gz \\
            --detect_adapter_for_pe \\
            --thread ${threads} \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            ${args} \\
            2>| >(tee ${prefix}.fastp.log >&2)
        """
    }

    stub:
    def prefix = task.ext.prefix ?: id
    def files  = single_end ? "${prefix}.fastp.fastq.gz" : "${prefix}_R1.fastp.fastq.gz ${prefix}_R2.fastp.fastq.gz"
    """
    for f in ${files}; do echo '' | gzip > \$f; done
    touch ${prefix}.fastp.json ${prefix}.fastp.html ${prefix}.fastp.log
    """
}
