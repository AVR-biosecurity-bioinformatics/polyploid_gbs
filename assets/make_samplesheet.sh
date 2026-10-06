#!/usr/bin/env bash
set -euo pipefail

# Usage: bash make_samplesheet.sh samplesheet.csv DIR [DIR ...]
if (( $# < 2 )); then
    echo "Usage: $0 samplesheet.csv FASTQ_DIR [FASTQ_DIR ...]" >&2
    exit 1
fi

output=$1
shift

for dir in "$@"; do
    [[ -d "$dir" ]] || { echo "Not a directory: $dir" >&2; exit 1; }
done

# Read all FASTQs from the supplied directories, including subdirectories.
mapfile -d '' -t fastqs < <(
    find "$@" -type f \
        \( -name '*_1.fastq.gz' -o -name '*_2.fastq.gz' \) \
        -print0 | sort -z
)

(( ${#fastqs[@]} > 0 )) || {
    echo "No paired-end FASTQs found" >&2
    exit 1
}

declare -A seen_library
printf 'sample_name,read1,read2\n' > "$output"

for fastq in "${fastqs[@]}"; do
    [[ "$fastq" == *_1.fastq.gz ]] || continue

    read1=$(readlink -f -- "$fastq")
    read2="${read1%_1.fastq.gz}_2.fastq.gz"

    [[ -f "$read2" ]] || {
        echo "Missing R2 for: $read1" >&2
        exit 1
    }

    filename=${read1##*/}
    library=${filename%_1.fastq.gz}
    sample_name=${library%_*}  # Remove the final barcode, retain the date.

    # Your current Nextflow parser derives library_id from the R1 filename.
    # Identical R1 filenames in different directories would collide.
    if [[ -v seen_library["$library"] ]]; then
        echo "Duplicate library filename across directories: $filename" >&2
        exit 1
    fi
    seen_library["$library"]=1

    printf '%s,%s,%s\n' "$sample_name" "$read1" "$read2" >> "$output"
done

# Also catch R2 files that have no corresponding R1.
for fastq in "${fastqs[@]}"; do
    [[ "$fastq" == *_2.fastq.gz ]] || continue
    [[ -f "${fastq%_2.fastq.gz}_1.fastq.gz" ]] || {
        echo "Missing R1 for: $fastq" >&2
        exit 1
    }
done

echo "Wrote $(($(wc -l < "$output") - 1)) pairs to $output"
