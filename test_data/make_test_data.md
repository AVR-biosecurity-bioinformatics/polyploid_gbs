

# Make test loci list

```
REFERENCE="/group/pathogens/IAWS/Projects/Biocontrol/gorse/AGRF_NXGSQCAGRF25080393-2_23MVYYLT3_gbs-2/NXGSQCAGRF25080393-2_consensus.fa"

awk '/^>/ { n++; if (n > 100) exit } n > 0 { print }' \
    "$REFERENCE" > consensus_100.fa

grep -c '^>' consensus_100.fa

```

# Make test fastqs

Extract just reads aligned to these loci for 4 samples from bams
```
module load SAMtools/1.23.1-GCC-13.3.0
bash test_data/make_test_reads.sh consensus_100.fa \
    output/results/bam/12122025_RK_001_03-1E.bam \
    output/results/bam/20260122_JE__001_03-1E.bam \
    output/results/bam/12122025_RK_001_02-1E.bam \
    output/results/bam/251120_GL_001_01-1E.bam

```