# polyploid ddRAD variant-calling pipeline

This Nextflow pipeline aligns paired-end ddRAD reads to a FASTA of target loci and calls variants jointly across samples. It trims reads with fastp, aligns them with Minibwa, merges libraries that share a `sample_name`, calls variants with FreeBayes across parallel groups of whole loci, and generates Riker, bcftools and MultiQC statistics.

## Requirements

Run the pipeline on BASC with Nextflow, SLURM and Conda available. The pipeline uses process-specific Conda packages defined in `main.nf`.

## 1. Install nextflow

Install a [self-install `nextflow` package](https://www.nextflow.io/docs/latest/install.html#self-install) 

# 2. clone the repository

```
analysis_dir=/path/to/my/analysis/dir

# clone Freyr repository
git clone https://github.com/AVR-biosecurity-bioinformatics/polyploid_gbs.git $analysis_dir 
```

## 3. Prepare a samplesheet

The samplesheet is a CSV with exactly these columns:

```csv
sample_name,read1,read2
12122025_RK_001_02-1E,/path/to/sample_R1.fastq.gz,/path/to/sample_R2.fastq.gz
```

Each row represents one paired-end read set. Use the same `sample_name` on multiple rows when those read sets belong to the same biological sample and should be merged after alignment.

A samplesheet can be generate from one or more FASTQ directories using the included helper `assets/make_samplesheet.sh` or generated manually

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
    -output-dir /group/pathogens/IAWS/Projects/Biocontrol/gorse/gorse_diploid \ 
    -resume
```

## Run the full dataset (Hexaploid calling)

From the pipeline directory:

```bash
REFERENCE="/group/pathogens/IAWS/Projects/Biocontrol/gorse/AGRF_NXGSQCAGRF25080393-2_23MVYYLT3_gbs-2/NXGSQCAGRF25080393-2_consensus.fa"

nextflow run . \
    --reference "$REFERENCE" \
    --samplesheet gorse_samplesheet.csv \
    --slurm_account 'fruitfly' \
    --nchunks 1000 \
    --ploidy 6 \
    -output-dir /group/pathogens/IAWS/Projects/Biocontrol/gorse/gorse_hexaploid \ 
    -resume
```

Replace `fruitfly` with the SLURM account you are authorised to use. The number of parallel FreeBayes chunks defaults to the value in `main.nf`; override it with `--nchunks N` if needed.

## Changing pipeline parameters

Pipeline parameters are declared in the `params` block at the top of `main.nf`. `reference`, `samplesheet` and  `slurm_account` are required. All other parameters have defaults and can be overridden for an individual run using
`--parameter_name value` or by directly editing the `main.nf` file

For example, to increase the mapping quality to 30 and a min_alernat_fraction of 10:
```bash
nextflow run . \
    --reference "$REFERENCE" \
    --samplesheet gorse_samplesheet.csv \
    --slurm_account fruitfly \
    --ploidy 2 \
    --min_mapping_quality 30 \
    --min_alternate_fraction 0.10 \
    -resume
```

## Full parameter list
| Parameter | Type | Required | Description |
|------------|------|----------|-------------|
| `reference` | `Path` | Yes | FASTA containing the ddRAD target loci. |
| `samplesheet` | `String` | Yes | Path to the sample sheet used to locate input reads. |
| `slurm_account` | `String` | Yes | SLURM account for submitted jobs. |
| `nchunks` | `Integer` | `200` | Number of groups of whole loci called in parallel. |
| `lower_percentile` | `Integer` | `1` | Exclude loci below this percentile. |
| `upper_percentile` | `Integer` | `99` | Exclude loci above this percentile. |
| `ploidy` | `Integer` | `2` | Number of chromosome copies per sample. |
| `theta` | `Double` | `0.001` | Expected population diversity used by the population prior. |
| `pooled_discrete` | `Boolean` | `false` | Model samples as discrete pools rather than individuals. |
| `hwe_priors_off` | `Boolean` | `true` | Disable the Hardy-Weinberg equilibrium prior. |
| `min_mapping_quality` | `Integer` | `20` | Ignore reads with mapping quality below this threshold. |
| `min_base_quality` | `Integer` | `20` | Ignore allele observations with base quality below this threshold. |
| `min_supporting_allele_qsum` | `Integer` | `0` | Minimum summed quality supporting an allele. |
| `min_coverage` | `Integer` | `10` | Minimum site coverage required for FreeBayes to process a site. |
| `limit_coverage` | `Integer` | `100000` | Coverage limit applied during variant calling. |
| `min_alternate_fraction` | `Double` | `0.05` | Minimum alternate allele fraction within a sample to evaluate a site. |
| `min_alternate_count` | `Integer` | `2` | Minimum alternate-supporting observations within a sample. |
| `use_best_n_alleles` | `Integer` | `4` | Evaluate at most the four best-supported SNP alleles. |
| `no_partial_observations` | `Boolean` | `false` | Exclude reads that do not span the detection window. |


## Outputs

The default output directory is `output/results/`, as configured in `nextflow.config`:

- `bam/`: final per-sample BAMs and indexes.
- `vcf/`: joint cohort VCF and index.
- `qc/bam_stats/`: raw per-sample Riker files.
- `qc/vcf_stats/`: bcftools VCF statistics.
- `qc/`: MultiQC report.

Use `-output-dir /path/to/results` to override the output directory for a run. Intermediate files remain in Nextflow's `work/` directory; `-resume` reuses eligible completed tasks.