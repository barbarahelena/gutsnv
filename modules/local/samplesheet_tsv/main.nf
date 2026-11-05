nextflow.enable.types = true

/*
 * Write the sample metadata table (sample, cohort, subject, timepoint) used by the
 * per-species genotype QC and the summary.
 */
process SAMPLESHEET_TSV {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8' }"

    input:
    rows: List<String>

    output:
    file('samples.tsv')

    script:
    """
    cat <<'END_SAMPLES' > samples.tsv
sample\tcohort\tsubject\ttimepoint
${rows.toSorted().join('\n')}
END_SAMPLES
    """

    stub:
    """
    cat <<'END_SAMPLES' > samples.tsv
sample\tcohort\tsubject\ttimepoint
${rows.toSorted().join('\n')}
END_SAMPLES
    """
}
