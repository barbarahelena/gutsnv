nextflow.enable.types = true

include { Reads             } from '../../../types/records'
include { GTPRO_GENOTYPE    } from '../../../modules/local/gtpro_genotype'
include { GTPRO_SPLITDICT   } from '../../../modules/local/gtpro_splitdict'
include { SPECIES_GENOTYPES } from '../../../modules/local/species_genotypes'
include { KINSHIP           } from '../../../modules/local/kinship'

/*
 * SNV metagenotyping with GT-Pro, per-species allele-frequency matrices with strain
 * quasi-phaseability and per-cohort QC, and per-species strain kinship matrices.
 *
 * Allele counts are split into species bins (species_id modulo n_bins), so the
 * per-species work runs in n_bins parallel jobs without one huge merge step.
 */
workflow GENOTYPING {
    take:
    reads: Channel<Reads>
    samples_tsv: Value<Path>
    gtpro_db: Path
    gtpro_db_name: String
    gtpro_snp_dict: Path
    n_bins: Integer
    qc_opts: Map
    skip_kinship: Boolean
    kinship_opts: Map

    main:
    gtpro = GTPRO_GENOTYPE(reads, gtpro_db, gtpro_db_name, n_bins)

    dict_bins = GTPRO_SPLITDICT(channel.of(gtpro_snp_dict), n_bins)
        .flatMap { files -> files.collect { f -> record(bin: binOf(f), dict: f) } }

    count_bins = gtpro
        .flatMap { r -> r.bins.collect { f -> tuple(binOf(f), f) } }
        .groupBy()
        .map { bin, files -> record(bin: bin, counts: files.toSet()) }

    species_bins = SPECIES_GENOTYPES(count_bins.join(dict_bins, by: 'bin'), samples_tsv, qc_opts)

    af      = species_bins.flatMap { r -> r.af.collect { f -> record(species: speciesOf(f), af: f) } }
    snvs    = species_bins.flatMap { r -> r.snvs.collect { f -> record(species: speciesOf(f), snvs: f) } }
    samples = species_bins.flatMap { r -> r.samples.collect { f -> record(species: speciesOf(f), samples: f) } }
    genotypes = af.join(snvs, by: 'species').join(samples, by: 'species')

    kinship = skip_kinship ? channel.empty() : KINSHIP(genotypes, kinship_opts)

    emit:
    gtpro: Channel<Record> = gtpro
    species_bins: Channel<Record> = species_bins
    genotypes: Channel<Record> = genotypes
    kinship: Channel<Record> = kinship
}

// "<sample>.bin007.tsv.gz" / "dict.bin007.tsv.gz" -> "007" (sample names may contain dots)
def binOf(f: Path) -> String {
    def parts = f.name.tokenize('.')
    return parts[parts.size() - 3].substring(3)
}

// "<species>.af.tsv.gz" -> "<species>"
def speciesOf(f: Path) -> String {
    return f.name.tokenize('.')[0]
}
