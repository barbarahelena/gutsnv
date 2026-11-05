nextflow.enable.types = true

/*
 * One row of the input samplesheet (--input). Columns are matched by name.
 * Samples sequenced over several runs/lanes have one row per run with the same `sample`.
 */
record Sample {
    sample: String
    fastq_1: Path
    fastq_2: Path?
    cohort: String?
    subject: String?
    timepoint: String?
}

/*
 * All runs of one sample (after grouping samplesheet rows by `sample`).
 */
record RunGroup {
    id: String
    single_end: Boolean
    reads_1: List<Path>
    reads_2: List<Path>
}

/*
 * The reads of one sample: [R1, R2] for paired-end, [R] for single-end data.
 */
record Reads {
    id: String
    single_end: Boolean
    reads: List<Path>
}

def toReads(id: String, single_end: Boolean, reads: Iterable<Path>) -> Reads {
    return record(id: id, single_end: single_end, reads: reads.toList().toSorted())
}
