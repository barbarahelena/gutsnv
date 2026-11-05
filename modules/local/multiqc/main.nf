nextflow.enable.types = true

/*
 * Aggregate QC reports with MultiQC.
 * Software versions are passed in as YAML text.
 * NOTE: this process must not publish to the `versions` topic: its input depends on that
 * topic being closed, so publishing to it would make the pipeline hang.
 */
process MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc' }"

    input:
    reports: Set<Path>
    config: Path
    extra_config: Path?
    logo: Path?
    versions: String

    stage:
    stageAs reports, '?/*'

    output:
    record(
        report: file('multiqc_report.html'),
        data: file('multiqc_data'),
        plots: file('multiqc_plots', optional: true),
        versions: file('software_versions.yml')
    )

    script:
    def args      = task.ext.args ?: ''
    def extra     = extra_config ? "--config ${extra_config}" : ''
    def logo_conf = logo ? "--cl-config 'custom_logo: \"${logo}\"'" : ''
    """
    cat <<'END_VERSIONS' > software_versions.yml
${versions}
END_VERSIONS
    cp software_versions.yml software_versions_mqc.yml

    multiqc \\
        --force \\
        ${args} \\
        --config ${config} \\
        ${extra} \\
        ${logo_conf} \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_report.html
    cat <<'END_VERSIONS' > software_versions.yml
${versions}
END_VERSIONS
    """
}
