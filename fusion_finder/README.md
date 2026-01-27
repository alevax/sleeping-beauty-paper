# Sleeping Beauty Transposon RNA-seq Analysis

RNA-seq alignment pipeline for detecting Sleeping Beauty (T2/Onc2) transposon insertion sites in mouse tumor samples.

## Methods

Based on Temiz et al. (Genome Research, 2013):
- mm10 genome with En2 gene masked (chr5:28,148,329-28,157,139)
- T2/Onc2 transposon sequence (2163bp) appended as chromosome `chrSB`
- HISAT2 alignment with splice-site awareness from GENCODE annotations
- FUSION_FINDER (Sarver, 2013) for identifying transposon-genome junction sites
- **Paired-end data only**: FUSION_FINDER relies on discordant read pairs to detect fusions

## Package Versions

| Package | Version |
|---------|---------|
| HISAT2 | 2.2.1 |
| samtools | 1.17 |
| bedtools | 2.29.2 |
| GENCODE | vM25 (mm10/GRCm38) |
| Perl | 5.30.0 |

## Directory Structure

```
fusion_finder/
├── README.md
├── input.txt                # FUSION_FINDER input config
├── loc_SB.txt               # chrSB regions for fusion detection
├── data/
│   └── fastq/
│       ├── samples.txt      # Sample metadata (ID, paired)
│       └── *.fastq.gz       # Raw FASTQ files (gitignored)
├── reference/
│   ├── chrSB.fa             # T2/Onc2 transposon sequence (committed)
│   ├── mm10.fa              # Mouse genome (gitignored)
│   ├── mm10_chrSB.fa        # Combined reference (gitignored)
│   ├── gencode.vM25.annotation.gtf (gitignored)
│   └── hisat2_index/        # HISAT2 index files (gitignored)
├── results/
│   ├── bams/                # Aligned BAM files (gitignored)
│   └── fusions/
│       └── fusions_all.txt  # Fusion detection results (input for downstream analysis)
├── scripts/                 # Preprocessing pipeline
│   ├── build_reference.sh
│   ├── align_samples.sh
│   ├── slurm_build_reference.sh
│   ├── slurm_align_array.sh
│   ├── fusion_finder.pl
│   ├── run_fusion_finder.sh
│   └── slurm_fusion_finder.sh
└── logs/                    # SLURM job logs (gitignored)

# CIS data moved to central location
../data/cis-results/
├── all_insertions_master.csv   # CIS insertions (all genes)
└── key_genes_insertions.csv    # CIS insertions (key genes)

# Downstream analysis scripts
../libs/fusion_finder_analysis.R  # Consolidated analysis script

# Analysis outputs
../experiments/fusion-finder/
├── fusion-finder-analysis-log.txt
├── figures/
│   ├── fusions_by_tumor_type.pdf
│   └── fusions_by_transposon_region.pdf
└── tables/
    ├── all_genes_cis_fusions.csv
    └── key_genes_cis_fusions.csv
```

## Samples

29 paired-end samples (CMZ327-CMZ357).

See `data/fastq/samples.txt` for sample IDs.

## Usage

### Step 1: Build Reference
```bash
# Local
bash scripts/build_reference.sh

# SLURM (recommended)
sbatch scripts/slurm_build_reference.sh
```

### Step 2: Align Samples
```bash
# Local
bash scripts/align_samples.sh

# SLURM (recommended - runs 29 samples in parallel)
sbatch scripts/slurm_align_array.sh
```

### Step 3: Find Fusions
```bash
# Local
bash scripts/run_fusion_finder.sh

# SLURM
sbatch scripts/slurm_fusion_finder.sh
```

## Downstream Analysis

The fusion finder analysis is integrated into the main paper analysis pipeline. From the repository root:

```bash
# Run full pipeline including fusion analysis
bash main.sh

# Or run fusion analysis only
Rscript libs/fusion_finder_analysis.R
```

**Outputs** (written to `experiments/fusion-finder/`):
- `tables/all_genes_cis_fusions.csv` - CIS-fusion overlap for all genes
- `tables/key_genes_cis_fusions.csv` - CIS-fusion overlap for key cancer genes
- `figures/fusions_by_tumor_type.pdf` - Fusion distribution by tumor type
- `figures/fusions_by_transposon_region.pdf` - Fusion distribution by transposon region

**Requirements**:
- Completed preprocessing (Steps 1-3 above)
- CIS data files in `data/cis-results/`
- Fusion detection results (`results/fusions/fusions_all.txt`)
