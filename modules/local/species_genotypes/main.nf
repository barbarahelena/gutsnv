nextflow.enable.types = true

/*
 * Per-species allele-frequency matrices, strain quasi-phaseability and per-cohort SNV QC
 * for all species in one bin. See templates/species_genotypes.R for definitions.
 */
process SPECIES_GENOTYPES {
    tag "bin${bin}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/findutils_gzip_r-ape_r-base_r-data.table:2ecb840cdd227774'
        : 'community.wave.seqera.io/library/findutils_gzip_r-ape_r-base_r-data.table:c788f6d95371f0cf' }"

    input:
    record(bin: String, counts: Set<Path>, dict: Path)
    samples: Path
    opts: Map

    stage:
    stageAs counts, 'counts/*'

    output:
    record(
        bin: bin,
        af: files('*.af.tsv.gz', optional: true),
        snvs: files('*.snvs.tsv.gz', optional: true),
        samples: files('*.samples.tsv.gz', optional: true),
        sample_species: file("bin${bin}.sample_species.tsv.gz"),
        species_cohort: file("bin${bin}.species_cohort.tsv")
    )

    topic:
    file('versions.tsv') >> 'versions_files'

    script:
    template 'species_genotypes.R'

    stub:
    """
    printf 'snv_id\\ts1\\n100${bin}_102\\t0\\n' | gzip > 100${bin}.af.tsv.gz
    printf 'snv_id\\tspecies\\tpos\\n100${bin}_102\\t100${bin}\\t102\\n' | gzip > 100${bin}.snvs.tsv.gz
    printf 'sample\\tcohort\\tspecies\\tquasi_phaseable\\ns1\\tall\\t100${bin}\\tTRUE\\n' | gzip > 100${bin}.samples.tsv.gz
    printf 'sample\\tcohort\\tspecies\\tquasi_phaseable\\ns1\\tall\\t100${bin}\\tTRUE\\n' | gzip > bin${bin}.sample_species.tsv.gz
    printf 'species\\tcohort\\tspecies_pass\\n100${bin}\\tall\\tTRUE\\n' > bin${bin}.species_cohort.tsv
    printf '${task.process}\\tr-base\\tstub\\n' > versions.tsv
    """
}
