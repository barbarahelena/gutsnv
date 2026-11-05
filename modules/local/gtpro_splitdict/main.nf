nextflow.enable.types = true

/*
 * Split the GT-Pro SNP dictionary (species, global_pos, contig, local_pos, ref, alt)
 * into the same species bins as the allele counts, keeping species, position and alleles.
 */
process GTPRO_SPLITDICT {
    tag snp_dict.name
    label 'process_single'

    // GT-Pro is not available from Bioconda; image built from the Dockerfile in the container repository
    // docker.io/barbarahelena/gt-pro:1.0.1 = sha256:979f24d50aaeea090b979fc686253d62b56858917a1802c58f95199109274373
    container 'docker.io/barbarahelena/gt-pro:1.0.1'

    input:
    snp_dict: Path
    n_bins: Integer

    output:
    files('dict.bin*.tsv.gz')

    script:
    def reader = snp_dict.name.endsWith('.gz') ? 'gzip -dc' : 'cat'
    """
    ${reader} ${snp_dict} \\
    | awk -v nb=${n_bins} 'BEGIN { OFS = "\\t" }
        \$1 ~ /^[0-9]+\$/ {
            print \$1, \$2, \$5, \$6 | ("gzip -c > " sprintf("dict.bin%03d.tsv.gz", (\$1 + 0) % nb))
        }'
    """

    stub:
    """
    printf '100003\\t101\\tA\\tG\\n' | gzip > dict.bin003.tsv.gz
    """
}
