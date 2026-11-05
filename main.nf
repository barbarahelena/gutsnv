#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    gutsnv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Metagenomic SNV genotyping (GT-Pro), strain quasi-phaseability, per-species QC and
    strain kinship matrices for microbial SNV association studies.
    Github : https://github.com/barbarahelena/gutsnv
----------------------------------------------------------------------------------------
*/

nextflow.enable.types = true

include { Sample          } from './types/records'
include { PREPROCESSING   } from './subworkflows/local/preprocessing'
include { GENOTYPING      } from './subworkflows/local/genotyping'
include { COVERM_GENOME   } from './modules/local/coverm_genome'
include { SAMPLESHEET_TSV } from './modules/local/samplesheet_tsv'
include { SUMMARISE       } from './modules/local/summarise'
include { MULTIQC         } from './modules/local/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {
    // Input: samplesheet (CSV) with columns sample, fastq_1, fastq_2, cohort, subject, timepoint
    input: Channel<Sample>

    // Read preprocessing
    skip_fastqc: Boolean = false
    skip_qc: Boolean = false                // skip adapter/quality trimming with fastp
    skip_host_removal: Boolean = false      // skip host-read removal with Bowtie2
    host_fasta: Path?                       // host genome FASTA (an index is built)
    host_index: Path?                       // or: directory with a pre-built Bowtie2 index
    save_clean_reads: Boolean = false

    // Species abundance
    run_coverm: Boolean = false             // also estimate abundance with CoverM (a GT-Pro-based estimate is always computed)
    abundance_genomes: Path?                // directory with one reference genome FASTA per species (for CoverM)
    abundance_genome_extension: String = 'fna'
    genome2species: Path?                   // optional TSV without header: genome<TAB>species_id

    // GT-Pro metagenotyping
    gtpro_db: Path                          // directory with the GT-Pro database files
    gtpro_db_name: String = '20190723_881species'
    gtpro_snp_dict: Path                    // GT-Pro SNP dictionary (variants_main.covered.hq.snp_dict.tsv)
    save_gtpro_raw: Boolean = false
    species_bins: Integer = 64              // number of parallel species bins
    pos_offset: Integer = 1                 // added to GT-Pro global positions in SNV IDs

    // Allele frequencies and strain quasi-phaseability
    min_depth: Integer = 1                  // minimum reads at a site for an allele frequency
    qp_min_depth: Integer = 2               // minimum reads at a site for quasi-phaseability
    qp_site_major_freq: Float = 0.9         // a site is 'pure' if the major allele has > this fraction of reads
    qp_min_frac: Float = 0.8                // a sample is quasi-phaseable if > this fraction of sites are pure
    qp_min_sites: Integer = 100             // minimum sites to assess quasi-phaseability

    // Per-cohort SNV and species QC
    snv_max_missing: Float = 0.9            // maximum fraction of cohort samples without data
    snv_min_maf: Float = 0.05
    species_min_snvs: Integer = 1000
    min_species_samples: Integer = 10

    // Strain kinship
    skip_kinship: Boolean = false
    kinship_qp_only: Boolean = true         // compute kinship on quasi-phaseable samples only
    kinship_min_shared: Integer = 1000      // minimum shared SNVs per sample pair
    kinship_max_snvs: Integer = 0           // subsample SNVs for the distance (0 = use all)
    kinship_n_pcs: Integer = 10

    // MultiQC
    multiqc_config: Path?
    multiqc_logo: Path?
    multiqc_title: String?
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    main:
    if (!params.skip_host_removal && !params.host_fasta && !params.host_index) {
        error("Host removal requires --host_fasta or --host_index (or use --skip_host_removal)")
    }
    if (params.run_coverm && !params.abundance_genomes) {
        error("CoverM (--run_coverm) requires --abundance_genomes")
    }

    samples = params.input

    sample_rows = samples
        .map { s -> [s.sample, s.cohort ?: 'all', s.subject ?: s.sample, s.timepoint ?: ''].join('\t') }
        .unique()
        .collect()
    samples_tsv = SAMPLESHEET_TSV(sample_rows.map { rows -> rows.toList() })

    //
    // Read preprocessing
    //
    pre = PREPROCESSING(
        samples,
        params.skip_fastqc,
        params.skip_qc,
        params.skip_host_removal,
        params.host_fasta,
        params.host_index
    )

    //
    // Metagenotyping, per-species QC and kinship
    //
    geno = GENOTYPING(
        pre.reads,
        samples_tsv,
        params.gtpro_db,
        params.gtpro_db_name,
        params.gtpro_snp_dict,
        params.species_bins,
        [
            pos_offset: params.pos_offset,
            min_depth: params.min_depth,
            qp_min_depth: params.qp_min_depth,
            qp_site_major_freq: params.qp_site_major_freq,
            qp_min_frac: params.qp_min_frac,
            qp_min_sites: params.qp_min_sites,
            min_species_samples: params.min_species_samples,
            snv_max_missing: params.snv_max_missing,
            snv_min_maf: params.snv_min_maf,
            species_min_snvs: params.species_min_snvs
        ],
        params.skip_kinship,
        [
            qp_only: params.kinship_qp_only,
            min_shared: params.kinship_min_shared,
            max_snvs: params.kinship_max_snvs,
            n_pcs: params.kinship_n_pcs
        ]
    )

    //
    // Species abundance
    //
    abundance = !params.run_coverm
        ? channel.empty()
        : COVERM_GENOME(pre.reads, params.abundance_genomes, params.abundance_genome_extension)

    //
    // Summary tables
    //
    summary = SUMMARISE(
        geno.species_bins.map { r -> r.species_cohort }.collect().map { fs -> fs.toSet() },
        geno.species_bins.map { r -> r.sample_species }.collect().map { fs -> fs.toSet() },
        geno.gtpro.map { r -> r.stats }.collect().map { fs -> fs.toSet() },
        pre.read_stats.map { r -> r.stats }.collect().map { fs -> fs.toSet() },
        abundance.map { r -> r.abundance }.collect().map { fs -> fs.toSet() },
        geno.kinship.map { r -> r.kinship_samples }.collect().map { fs -> fs.toSet() },
        samples_tsv,
        params.genome2species
    )

    //
    // Software versions and MultiQC
    //
    // Bash processes publish (process, tool, version) tuples to the `versions` topic;
    // R template processes write them to a versions.tsv file in the `versions_files` topic
    version_files = channel.topic('versions_files')
        .flatMap { f -> f.readLines().collect { line -> line.tokenize('\t') } }
    versions = channel.topic('versions')
        .mix(version_files)
        .map { entry -> tuple(entry[0].toString().tokenize(':').last(), "    ${entry[1]}: ${entry[2]}") }
        .unique()
        .groupBy()
        .map { proc, tools -> "${proc}:\n" + tools.toList().toSorted().join('\n') }
        .collect()
        .map { blocks ->
            blocks.toList().toSorted().join('\n') +
                "\nWorkflow:\n    ${workflow.manifest.name}: ${workflow.manifest.version}\n    Nextflow: ${nextflow.version}"
        }

    multiqc_reports = pre.fastqc.flatMap { r -> r.zip }
        .mix(pre.fastp.map { r -> r.json })
        .mix(pre.host_removal.map { r -> r.log })
        .mix(summary.flatMap { r -> r.multiqc })

    multiqc = MULTIQC(
        multiqc_reports.collect().map { fs -> fs.toSet() },
        file("${projectDir}/assets/multiqc_config.yml"),
        params.multiqc_config,
        params.multiqc_logo,
        versions
    )

    publish:
    fastqc       = pre.fastqc
    fastp        = pre.fastp.map { r -> record(id: r.id, json: r.json, html: r.html, log: r.log) }
    host_removal = pre.host_removal.map { r -> record(id: r.id, log: r.log) }
    clean_reads  = params.save_clean_reads ? pre.reads : channel.empty()
    read_stats   = pre.read_stats
    gtpro        = geno.gtpro.map { r -> record(id: r.id, stats: r.stats) }
    gtpro_raw    = params.save_gtpro_raw ? geno.gtpro.map { r -> record(id: r.id, raw: r.raw) } : channel.empty()
    genotypes    = geno.genotypes
    kinship      = geno.kinship
    abundance    = abundance
    summary      = summary
    multiqc      = multiqc

    onComplete:
    log.info(
        workflow.success
            ? "Pipeline completed successfully. Results in: ${workflow.outputDir}"
            : "Pipeline completed with errors"
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    fastqc: Channel<Record> {
        path 'qc/fastqc'
    }
    fastp: Channel<Record> {
        path 'qc/fastp'
    }
    host_removal: Channel<Record> {
        path 'qc/host_removal'
    }
    clean_reads: Channel<Record> {
        path 'clean_reads'
    }
    read_stats: Channel<Record> {
        path 'qc/read_stats'
    }
    gtpro: Channel<Record> {
        path 'gtpro/stats'
    }
    gtpro_raw: Channel<Record> {
        path { r -> "gtpro/raw/${r.id}" }
    }
    genotypes: Channel<Record> {
        path { r -> "genotypes/${r.species}" }
        index {
            path 'genotypes/index.csv'
            header true
        }
    }
    kinship: Channel<Record> {
        path { r -> "kinship/${r.species}" }
        index {
            path 'kinship/index.csv'
            header true
        }
    }
    abundance: Channel<Record> {
        path 'abundance/coverm'
    }
    summary: Record {
        path 'summary'
    }
    multiqc: Record {
        path 'multiqc'
    }
}
