#!/usr/bin/env bash
set -euo pipefail

# Usage: bash make_test_reads.sh consensus_100.fa sample1.bam ... sample4.bam
if (( $# != 5 )); then
    echo "Usage: $0 consensus_100.fa sample1.bam sample2.bam sample3.bam sample4.bam" >&2
    exit 1
fi

reference=$1
shift

outdir=test_100_loci
mkdir -p "$outdir"

# The .fai provides each locus name and length.
samtools faidx "$reference"
awk 'BEGIN { OFS = "\t" } { print $1, 0, $2 }' \
    "${reference}.fai" > "$outdir/loci.bed"

printf 'sample_name,read1,read2\n' > "$outdir/samplesheet.csv"

for bam in "$@"; do
    [[ -f "$bam" ]] || { echo "Missing BAM: $bam" >&2; exit 1; }

    sample_name=$(basename "$bam" .bam)
    names="$outdir/${sample_name}.read_names.txt"
    read1="$outdir/${sample_name}_R1.fastq.gz"
    read2="$outdir/${sample_name}_R2.fastq.gz"

    # Find templates with an alignment overlapping any test locus.
    samtools view -M -L "$outdir/loci.bed" -F 0x900 "$bam" |
        cut -f1 |
        LC_ALL=C sort -u > "$names"

    if [[ ! -s "$names" ]]; then
        echo "No reads found on the test loci in: $bam" >&2
        exit 1
    fi

    # Retrieve both mates by read name, including mates that mapped elsewhere
    # or were unmapped. Collate before writing synchronised paired FASTQs.
    samtools view -b -N "$names" "$bam" |
        samtools collate -u -O - |
        samtools fastq \
            -1 "$read1" \
            -2 "$read2" \
            -0 /dev/null \
            -s /dev/null \
            -n -

    printf '%s,%s,%s\n' \
        "$sample_name" \
        "$(readlink -f "$read1")" \
        "$(readlink -f "$read2")" \
        >> "$outdir/samplesheet.csv"

    echo "$sample_name: $(wc -l < "$names") selected read names"
done

echo "Test samplesheet: $outdir/samplesheet.csv"