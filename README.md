# polyploid ddRAD variant-calling pipeline

This Nextflow pipeline aligns paired-end ddRAD reads to a FASTA of target loci and calls variants jointly across samples. It trims reads with fastp, aligns them with Minibwa, merges libraries that share a `sample_name`, calls variants with FreeBayes across parallel groups of whole loci, and generates Riker, bcftools and MultiQC statistics.

## Requirements

Run the pipeline on BASC with Nextflow, SLURM and Conda available. The pipeline uses process-specific Conda packages defined in `main.nf`.

## 1. Prepare a samplesheet

The samplesheet is a CSV with exactly these columns:

```csv
sample_name,read1,read2
12122025_RK_001_02-1E,/path/to/sample_R1.fastq.gz,/path/to/sample_R2.fastq.gz
```

Each row represents one paired-end read set. Use the same `sample_name` on multiple rows when those read sets belong to the same biological sample and should be merged after alignment.

A samplesheet can be generate from one or more FASTQ directories using the included helper `assets/make_samplesheet.sh`

```bash
bash assets/make_samplesheet.sh gorse_samplesheet.csv \
    /group/pathogens/IAWS/Projects/Biocontrol/gorse/AGRF_NXGSQCAGRF25080393-2_23MVYYLT3_gbs-2/fastq \
    /path/to/another/fastq_directory
```

Check the generated sample names and paths before running. The current pipeline parser resolves relative `read1` and `read2` paths against the **pipeline project directory**; absolute paths also work.

## Run a small test

The pipleine includes a small test dataset that can be run using:

```bash
nextflow run . -profile test --slurm_account 'fruitfly' -resume
```

## Run the full dataset (Diploid calling)

From the pipeline directory:

```bash
REFERENCE="/group/pathogens/IAWS/Projects/Biocontrol/gorse/AGRF_NXGSQCAGRF25080393-2_23MVYYLT3_gbs-2/NXGSQCAGRF25080393-2_consensus.fa"

nextflow run . \
    --reference "$REFERENCE" \
    --samplesheet gorse_samplesheet.csv \
    --slurm_account 'fruitfly' \
    --nchunks 1000 \
    --ploidy 2 \
    -output-dir gorse_diploid \ 
    -resume
```

Replace `fruitfly` with the SLURM account you are authorised to use. The number of parallel FreeBayes chunks defaults to the value in `main.nf`; override it with `--nchunks N` if needed.

## Changing pipeline parameters

Pipeline parameters are declared in the `params` block at the top of `main.nf`. `reference`, `samplesheet` and  `slurm_account` are required. All other parameters have defaults and can be overridden for an individual run using
`--parameter_name value` or by directly editing the `main.nf` file

For example, to call 50 chunks as diploids with a different minimum alternate-allele fraction:

```bash
nextflow run . \
    --reference "$REFERENCE" \
    --samplesheet gorse_samplesheet.csv \
    --slurm_account fruitfly \
    --nchunks 50 \
    --ploidy 2 \
    --min_alternate_fraction 0.10 \
    -resume
```

## Outputs

The default output directory is `output/results/`, as configured in `nextflow.config`:

- `bam/`: final per-sample BAMs and indexes.
- `vcf/`: joint cohort VCF and index.
- `qc/bam_stats/`: raw per-sample Riker files.
- `qc/vcf_stats/`: bcftools VCF statistics.
- `qc/`: MultiQC report.

Use `-output-dir /path/to/results` to override the output directory for a run. Intermediate files remain in Nextflow's `work/` directory; `-resume` reuses eligible completed tasks.