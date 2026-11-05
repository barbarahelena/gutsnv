nextflow.enable.types = true

include { Sample ; RunGroup ; Reads ; toReads } from '../../../types/records'
include { CAT_FASTQ           } from '../../../modules/local/cat_fastq'
include { FASTQC              } from '../../../modules/local/fastqc'
include { FASTP               } from '../../../modules/local/fastp'
include { BOWTIE2_BUILD       } from '../../../modules/local/bowtie2_build'
include { BOWTIE2_HOSTREMOVAL } from '../../../modules/local/bowtie2_hostremoval'
include { SEQKIT_STATS        } from '../../../modules/local/seqkit_stats'

/*
 * Read preprocessing: merge runs per sample, raw-read QC (FastQC), adapter/quality
 * trimming (fastp), host-read removal (Bowtie2) and read counting (seqkit).
 * Trimming and host removal can each be skipped, e.g. for reads that are already cleaned.
 */
workflow PREPROCESSING {
    take:
    samples: Channel<Sample>
    skip_fastqc: Boolean
    skip_qc: Boolean
    skip_host_removal: Boolean
    host_fasta: Path?
    host_index: Path?

    main:
    // One record per sample; runs of the same sample are concatenated
    grouped = samples
        .map { s -> tuple(s.sample, s) }
        .groupBy()
        .map { id, runs -> toRunGroup(id, runs.toList()) }

    one_run = grouped
        .filter { r -> r.reads_1.size() == 1 }
        .map { r -> toReads(r.id, r.single_end, r.reads_1 + r.reads_2) }
    merged = CAT_FASTQ(grouped.filter { r -> r.reads_1.size() > 1 })
        .map { r -> toReads(r.id, r.single_end, r.reads) }
    raw = one_run.mix(merged)

    fastqc = skip_fastqc ? channel.empty() : FASTQC(raw)

    fastp = skip_qc ? channel.empty() : FASTP(raw)
    trimmed = skip_qc
        ? raw
        : fastp.map { r -> toReads(r.id, r.single_end, r.reads) }

    host_removal = channel.empty()
    clean = trimmed
    if (!skip_host_removal) {
        index = host_index ? channel.value(host_index) : BOWTIE2_BUILD(host_fasta)
        host_removal = BOWTIE2_HOSTREMOVAL(trimmed, index)
        clean = host_removal.map { r -> toReads(r.id, r.single_end, r.reads) }
    }

    read_stats = SEQKIT_STATS(clean)

    emit:
    reads: Channel<Reads> = clean
    fastqc: Channel<Record> = fastqc
    fastp: Channel<Record> = fastp
    host_removal: Channel<Record> = host_removal
    read_stats: Channel<Record> = read_stats
}

def toRunGroup(id: String, rows: List<Sample>) -> RunGroup {
    def single_end = rows[0].fastq_2 == null
    if (rows.any { r -> (r.fastq_2 == null) != single_end }) {
        error("Sample '${id}' mixes single-end and paired-end runs")
    }
    // in paired-end mode every row has fastq_2; the fallback only satisfies the type checker
    def reads_1 = rows.collect { r -> r.fastq_1 }.toList()
    def reads_2 = single_end ? [] : rows.collect { r -> r.fastq_2 ?: r.fastq_1 }.toList()
    return record(id: id, single_end: single_end, reads_1: reads_1, reads_2: reads_2)
}
