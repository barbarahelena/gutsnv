#!/usr/bin/env bash
#
# Nextflow template for GTPRO_GENOTYPE.
#
# 1. Genotype each read file with GT-Pro.
# 2. Decode the raw GT-Pro output into per-site allele counts, split into species bins.
#
# Raw GT-Pro lines are "<snp_id> <count>", where snp_id = 6-digit species ID,
# 1 digit allele type (0 = reference, 1 = alternative), then the SNP's global position.
# Counts of all read files (R1 and R2) are summed per site.
#
# Output (no header), one gzipped file per bin (species_id modulo n_bins):
#   <prefix>.binNNN.tsv.gz : species  global_pos  sample  ref_count  alt_count
# plus <prefix>.gtpro_stats.tsv : sample  n_sites  n_species  total_count

set -euo pipefail

prefix="${task.ext.prefix ?: id}"

mkdir -p raw
i=0
for r in ${reads.join(' ')}; do
    i=\$((i + 1))
    GT_Pro genotype \\
        -d ${db}/${db_name} \\
        -t ${task.cpus} \\
        -o "raw/\${prefix}_R\${i}" \\
        ${task.ext.args ?: ''} \\
        "\$r"
done

decompress() {
    case "\$1" in
        *.lz4) lz4 -dc "\$1" ;;
        *.gz)  gzip -dc "\$1" ;;
        *.bz2) bzip2 -dc "\$1" ;;
        *)     cat "\$1" ;;
    esac
}

for f in raw/*; do decompress "\$f"; done \\
| awk -v sample="\${prefix}" -v nb="${n_bins}" '
    BEGIN { OFS = "\\t" }
    \$1 ~ /^[0-9]+\$/ && NF >= 2 {
        sp  = substr(\$1, 1, 6)
        typ = substr(\$1, 7, 1)
        pos = substr(\$1, 8) + 0
        key = sp "\\t" pos
        if (typ == "0") ref[key] += \$2; else alt[key] += \$2
        seen[key] = 1
        total += \$2
    }
    END {
        n = 0
        for (key in seen) {
            split(key, k, "\\t")
            sp_seen[k[1]] = 1
            out = sprintf("%s.bin%03d.tsv.gz", sample, (k[1] + 0) % nb)
            print k[1], k[2], sample, ref[key] + 0, alt[key] + 0 | ("gzip -c > " out)
            n++
        }
        ns = 0
        for (s in sp_seen) ns++
        printf "sample\\tn_sites\\tn_species\\ttotal_count\\n%s\\t%d\\t%d\\t%d\\n", sample, n, ns, total > (sample ".gtpro_stats.tsv")
    }
'
