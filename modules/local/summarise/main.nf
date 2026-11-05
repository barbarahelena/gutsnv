nextflow.enable.types = true

/*
 * Combine per-bin, per-sample and per-species statistics into summary tables that show
 * how many samples, quasi-phaseable samples and SNVs are available per species and cohort.
 * See templates/summarise.R.
 */
process SUMMARISE {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/r-ape_r-base_r-data.table:73189a23d396ba7a'
        : 'community.wave.seqera.io/library/r-ape_r-base_r-data.table:dc4136918d8aef77' }"

    input:
    species_cohort: Set<Path>
    sample_species: Set<Path>
    gtpro_stats: Set<Path>
    read_stats: Set<Path>
    abundance: Set<Path>
    kinship_samples: Set<Path>
    samples: Path
    genome2species: Path?

    stage:
    stageAs species_cohort, 'species_cohort/*'
    stageAs sample_species, 'sample_species/*'
    stageAs gtpro_stats, 'gtpro_stats/*'
    stageAs read_stats, 'read_stats/*'
    stageAs abundance, 'abundance/*'
    stageAs kinship_samples, 'kinship/*'

    output:
    record(
        sample_summary: file('sample_summary.tsv'),
        species_cohort_summary: file('species_cohort_summary.tsv'),
        sample_species: file('sample_species.tsv.gz'),
        species_abundance_gtpro: file('species_abundance_gtpro.tsv', optional: true),
        species_abundance: file('species_abundance.tsv', optional: true),
        abundance_comparison: file('abundance_comparison.tsv', optional: true),
        multiqc: files('*_mqc.tsv')
    )

    topic:
    file('versions.tsv') >> 'versions_files'

    script:
    template 'summarise.R'

    stub:
    """
    touch sample_summary.tsv species_cohort_summary.tsv species_abundance_gtpro.tsv gtpro_samples_mqc.tsv
    echo '' | gzip > sample_species.tsv.gz
    printf '${task.process}\\tr-base\\tstub\\n' > versions.tsv
    """
}
