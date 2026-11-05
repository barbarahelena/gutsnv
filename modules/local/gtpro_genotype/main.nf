nextflow.enable.types = true

/*
 * Metagenotyping of known gut microbial SNVs with GT-Pro (k-mer exact matching),
 * followed by decoding of the raw output into per-site reference/alternative allele
 * counts, split into species bins for downstream per-species processing.
 *
 * GT-Pro: https://github.com/zjshi/gt-pro (Shi et al., Nat Biotechnol 2022)
 */
process GTPRO_GENOTYPE {
    tag id
    label 'process_medium'

    // GT-Pro is not available from Bioconda; image built from the Dockerfile in the container repository
    // docker.io/barbarahelena/gt-pro:1.0.1 = sha256:979f24d50aaeea090b979fc686253d62b56858917a1802c58f95199109274373
    container 'docker.io/barbarahelena/gt-pro:1.0.1'

    input:
    record(id: String, reads: List<Path>)
    db: Path
    db_name: String
    n_bins: Integer

    output:
    record(
        id: id,
        bins: files("${task.ext.prefix ?: id}.bin*.tsv.gz"),
        stats: file("${task.ext.prefix ?: id}.gtpro_stats.tsv"),
        raw: files('raw/*')
    )

    topic:
    tuple(task.process, 'gt-pro', eval("GT_Pro -v | sed 's/GT-Pro version //'")) >> 'versions'

    script:
    template 'gtpro_genotype.sh'

    stub:
    def prefix = task.ext.prefix ?: id
    """
    mkdir raw
    touch raw/${prefix}_R1.tsv.lz4
    printf '100003\\t101\\t${prefix}\\t3\\t0\\n' | gzip > ${prefix}.bin003.tsv.gz
    printf 'sample\\tn_sites\\tn_species\\ttotal_count\\n${prefix}\\t1\\t1\\t3\\n' > ${prefix}.gtpro_stats.tsv
    """
}
