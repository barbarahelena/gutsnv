nextflow.enable.types = true

/*
 * Strain genetic distance, kinship matrix and PCoA for one species.
 * See templates/kinship.R for definitions.
 */
process KINSHIP {
    tag species
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/findutils_gzip_r-ape_r-base_r-data.table:2ecb840cdd227774'
        : 'community.wave.seqera.io/library/findutils_gzip_r-ape_r-base_r-data.table:c788f6d95371f0cf' }"

    input:
    record(species: String, af: Path, samples: Path)
    opts: Map

    output:
    record(
        species: species,
        kinship: file("${species}.kinship.tsv.gz", optional: true),
        distance: file("${species}.distance.tsv.gz", optional: true),
        nshared: file("${species}.nshared.tsv.gz", optional: true),
        pcs: file("${species}.pcs.tsv", optional: true),
        eigenvalues: file("${species}.pcoa_eigenvalues.tsv", optional: true),
        kinship_samples: file("${species}.kinship_samples.tsv")
    )

    topic:
    file('versions.tsv') >> 'versions_files'

    script:
    template 'kinship.R'

    stub:
    """
    for f in kinship distance nshared; do echo 'sample' | gzip > ${species}.\$f.tsv.gz; done
    touch ${species}.pcs.tsv ${species}.pcoa_eigenvalues.tsv
    printf 'sample\\tstatus\\n' > ${species}.kinship_samples.tsv
    printf '${task.process}\\tr-ape\\tstub\\n' > versions.tsv
    """
}
