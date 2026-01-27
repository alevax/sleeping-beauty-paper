# Isoform Switching Analysis

This directory contains the isoform switching analysis pipeline using **IsoformSwitchAnalyzeR**.

## Overview

Isoform switching occurs when the relative expression of isoforms from the same gene changes between conditions. This analysis identifies genes where isoform usage significantly differs between NEPC and nonNE conditions, and predicts the functional consequences of these switches.

## Pipeline Workflow

The IsoformSwitchAnalyzeR workflow is divided into two parts with 11 steps total:

### Part 1: Extract Isoform Switches and Sequences

| Step | Function | Description |
|------|----------|-------------|
| 1 | `importRdata()` | Import Salmon quantifications, GTF annotation, and FASTA sequences |
| 2 | `preFilter()` | Filter lowly expressed genes/isoforms |
| 3 | `isoformSwitchTestDEXSeq()` | Statistical test for differential isoform usage |
| 4 | `analyzeORF()` | Predict/annotate Open Reading Frames |

### Part 2: Annotate and Analyze Consequences

| Step | Function | Description |
|------|----------|-------------|
| 5 | **`analyzeAlternativeSplicing()`** | **CRITICAL** - Classifies alternative splicing events (IR, ES, A5, A3, ATSS, ATTS) |
| 6 | `analyzeSwitchConsequences()` | Predict functional consequences of isoform switches |
| 7 | `extractSwitchSummary()` | Extract significant isoform switches |
| 8 | Save results | Export RDS objects and CSV tables |
| 9 | Individual visualizations | Generate switch plots for top genes |
| 10 | Global splicing analysis | `extractSplicingSummary()`, `extractSplicingEnrichment()`, `extractSplicingGenomeWide()` |
| 11 | Global consequence analysis | `extractConsequenceSummary()`, `extractConsequenceEnrichment()`, `extractConsequenceGenomeWide()` |

### Critical Dependency

**`analyzeAlternativeSplicing()` MUST be called BEFORE `analyzeSwitchConsequences()`** if you want to analyze:
- `intron_retention`
- `intron_structure`

This is because intron retention events need to be classified before they can be evaluated as consequences.

## Comparisons

This analysis performs two comparisons:

1. **NEPC vs nonNE** - Compares neuroendocrine prostate cancer vs non-neuroendocrine within SBT2+ samples
2. **SBT2+ vs SBT2-** - Compares all SBT2-positive (NEPC + nonNE) vs SBT2-negative samples

## Directory Structure

```
isoform_switching/
├── R/                                   # Upstream analysis scripts (documentation)
│   ├── isoform_switch_analysis.R        # Main IsoformSwitchAnalyzeR pipeline
│   ├── 01_cis_isoform_switch_overlap.R  # CIS integration Part 1
│   └── 02_cis_switch_consequences.R     # CIS integration Part 2
├── scripts/
│   ├── build_salmon_index.sh            # Build Salmon index
│   ├── fetch_hmmer_results.py           # Fetch HMMER API results
│   └── slurm_*.sh                       # SLURM job scripts
├── data/                                # Input data (gitignored)
├── reference/                           # Reference files (gitignored)
├── results/                             # Pre-computed results (gitignored)
│   └── isoformswitch/cis_overlap/       # CIS-integrated analysis results
│       ├── checkpoint_cis_overlap.rds   # Used by publication script
│       ├── pfam_results.txt             # Used by publication script
│       ├── cis_switch_genes_detail.csv  # Used by publication script
│       └── cis_insertions_in_switching_genes.csv
├── docs/
│   └── IsoformSwitchAnalyzeR_vignette.html
└── README.md                            # This file

Note: data/, results/, and logs/ are gitignored but contain the upstream analysis
that generates pre-computed files used by the main paper analysis pipeline.
```

## Running the Analysis

### Submit SLURM Job

```bash
cd /lab/barcheese01/mdiberna/sleeping_beauty/isoform_switching
sbatch slurm/run_isoform_switch.sh
```

### Monitor Progress

```bash
# Check job status
squeue -u $USER

# View log file
tail -f logs/isoform_switch_*.log
```

## Checkpointing (Auto-Recovery)

The pipeline automatically saves checkpoints after expensive operations. If a job fails or is interrupted, re-running will automatically skip completed steps.

### Checkpoint Locations

| Checkpoint | After Step | Time Saved | File |
|------------|------------|------------|------|
| 1 | Import data | 5-15 min | `checkpoint_*_01_imported.rds` |
| 2 | DEXSeq testing | 10-30 min | `checkpoint_*_03_dexseq.rds` |
| 3 | Consequence analysis | 20-50 min | `checkpoint_*_06_consequences.rds` |

**Location:** `results/isoformswitch/checkpoints/`

### How It Works

- Before each expensive step, the script checks if a checkpoint exists
- If checkpoint exists → load it and skip to next step
- If no checkpoint → run the step and save checkpoint on success
- Checkpoints are comparison-specific (e.g., `NEPC_vs_nonNE`, `SBT2pos_vs_SBT2neg`)

### Force Fresh Re-run

To re-run from scratch, delete the checkpoint files:

```bash
rm -rf results/isoformswitch/checkpoints/
```

Or delete specific checkpoints:

```bash
# Re-run only DEXSeq and later steps for NEPC comparison
rm results/isoformswitch/checkpoints/checkpoint_NEPC_vs_nonNE_03_dexseq.rds
rm results/isoformswitch/checkpoints/checkpoint_NEPC_vs_nonNE_06_consequences.rds
```

## Input Requirements

1. **Salmon quantification** - `results/salmon/{sample}/quant.sf` for each sample
2. **Sample metadata** - CSV with `sample` and `condition` columns
3. **GTF annotation** - GENCODE vM25 annotation file
4. **Transcript FASTA** - GENCODE vM25 transcript sequences

## Output Files

### Per Comparison - Data Tables

| File | Description |
|------|-------------|
| `switchList_*.rds` | Complete R object with all analysis data |
| `top_switches_*.csv` | Significant isoform switches with statistics |
| `switching_genes_*.csv` | Gene-level summary of switches |
| `splicing_summary_*.csv` | Summary of alternative splicing events |
| `splicing_summary_global_*.csv` | Global splicing summary across all genes |
| `splicing_enrichment_*.csv` | Splicing type enrichment analysis |
| `splicing_genomewide_*.csv` | Genome-wide splicing patterns |
| `consequence_enrichment_*.csv` | Enrichment of functional consequences |
| `consequence_summary_global_*.csv` | Global consequence summary |
| `consequence_enrichment_global_*.csv` | Consequence gain/loss enrichment |
| `consequence_genomewide_*.csv` | Genome-wide consequence patterns |

### Per Comparison - Plots

| File | Description |
|------|-------------|
| `switch_*_GENE.pdf` | Individual gene switch visualization |

**Note**: Global summary plots are generated by the main IsoformSwitchAnalyzeR pipeline if needed.

## Interpretation

### Key Metrics

- **dIF (delta Isoform Fraction)**: Change in isoform usage between conditions (default cutoff: 0.1 = 10%)
- **q-value**: FDR-adjusted p-value from DEXSeq (default cutoff: 0.05)
- **IF (Isoform Fraction)**: Proportion of gene expression from each isoform

### Alternative Splicing Types

| Type | Abbreviation | Description |
|------|--------------|-------------|
| Intron Retention | IR | An intron is retained in the mature transcript |
| Exon Skipping | ES | One or more exons are skipped |
| Mutually Exclusive Exons | MEE | One of two exons is included, but not both |
| Alternative 5' Splice Site | A5 | Alternative donor site usage |
| Alternative 3' Splice Site | A3 | Alternative acceptor site usage |
| Alternative TSS | ATSS | Different transcription start site |
| Alternative TTS | ATTS | Different transcription termination site |

### Consequence Types

| Consequence | Description |
|-------------|-------------|
| `intron_retention` | Isoform retains an intron (often leads to NMD) |
| `coding_potential` | Change in predicted coding potential |
| `NMD_status` | Nonsense-mediated decay sensitivity |
| `ORF_seq_similarity` | Similarity of predicted protein sequences |
| `domains_identified` | Gain/loss of protein domains |
| `signal_peptide_identified` | Change in signal peptide prediction |

## Troubleshooting

### Error: "To test for intron retention alternative splicing must first be classified"

**Cause**: `analyzeSwitchConsequences()` was called with `intron_retention` in `consequencesToAnalyze` before running `analyzeAlternativeSplicing()`.

**Solution**: Add `analyzeAlternativeSplicing()` call before `analyzeSwitchConsequences()`:

```r
# Correct order:
switchList <- analyzeAlternativeSplicing(switchList)
switchList <- analyzeSwitchConsequences(switchList, consequencesToAnalyze = c('intron_retention', ...))
```

### Memory Issues

IsoformSwitchAnalyzeR can be memory-intensive. Request at least 64GB RAM for full analysis:

```bash
#SBATCH --mem=64G
```

## CIS-Integrated Analysis (Transposon Insertion Sites)

This analysis integrates Common Insertion Site (CIS) data from DNA-based transposon mapping with RNA isoform switching results. This allows identification of genes where transposon insertions may be causing functional changes through altered isoform usage.

### Why This Matters

Sleeping Beauty transposon insertions can disrupt genes in multiple ways:
- **Gene disruption**: Insertions that truncate or inactivate transcripts
- **Isoform switching**: Insertions that alter splicing or promoter usage
- **Domain loss/gain**: Changes in protein domain composition

By integrating CIS data with isoform switching analysis, we can identify genes where insertions correlate with functional isoform changes.

### CIS Integration Workflow

The analysis is split into two parts to allow for external Pfam analysis:

#### Part 1: Identify Overlapping Genes & Extract Sequences

```bash
Rscript R/01_cis_isoform_switch_overlap.R
```

This script:
1. Loads CIS insertion data from `../fusion_finder/R/all_insertions_master.csv`
2. Loads isoform switching results from checkpoints
3. Identifies genes with **both** CIS insertions AND significant isoform switches
4. Analyzes sample-level patterns (do samples with insertions show different isoform usage?)
5. Extracts amino acid sequences for overlapping genes
6. Outputs FASTA file for Pfam/HMMER analysis

**Output files:**
- `results/isoformswitch/cis_overlap/cis_switch_overlap_summary.csv`
- `results/isoformswitch/cis_overlap/cis_switch_genes_detail.csv`
- `results/isoformswitch/cis_overlap/sequences/cis_genes_AA.fasta`

#### External Step: Pfam Domain Analysis via HMMER Webserver

Upload the AA sequences to HMMER for protein domain prediction:

1. Go to: https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan
2. Switch to "Upload a File" tab
3. Upload `results/isoformswitch/cis_overlap/sequences/cis_genes_AA.fasta`
4. Select **Pfam** as target database
5. Submit job

**Important**: The webserver processes each sequence separately, creating individual job IDs. To download all results:

**Option A: Use the API fetch script (recommended)**
```bash
# After submission, note your job ID from the URL (e.g., d44e5578-851a-465f-ba4e-a820d8e10efc)
# First, get the list of sub-job IDs:
curl -s "https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/YOUR_JOB_ID" \
    -H "Accept: application/json" > results/isoformswitch/cis_overlap/hmmer_jobs.json

# Then fetch all results and convert to IsoformSwitchAnalyzeR format:
python3 scripts/fetch_hmmer_results.py \
    results/isoformswitch/cis_overlap/hmmer_jobs.json \
    results/isoformswitch/cis_overlap/pfam_results.txt
```

**Option B: Run from project directory with defaults**
```bash
# If hmmer_jobs.json already exists in the default location:
cd /lab/barcheese01/mdiberna/sleeping_beauty/isoform_switching
python3 scripts/fetch_hmmer_results.py
```

**Option C: Email results** (for small jobs)
- Provide your email when submitting
- Copy/paste results from email into a text file
- Save as `pfam_results.txt`

The fetch script (`scripts/fetch_hmmer_results.py`) converts HMMER API JSON to the tab-delimited format expected by IsoformSwitchAnalyzeR's `analyzePFAM()` function.

#### Part 2: Consequence Analysis with Domain Information

```bash
Rscript R/02_cis_switch_consequences.R results/isoformswitch/cis_overlap/pfam_results.txt
```

This script:
1. Imports Pfam domain predictions
2. Runs `analyzeSwitchConsequences()` with domain analysis enabled
3. Generates switch plots for each CIS-overlapping gene
4. Creates summary tables of functional consequences

**Output files:**
- `results/isoformswitch/cis_overlap/switchList_with_pfam.rds` - R object with Pfam domain data
- `results/isoformswitch/cis_overlap/domain_switches.csv` - Domains in switching isoforms
- `results/isoformswitch/cis_overlap/plots/*.pdf` - Switch plots for each CIS-overlapping gene
- `results/isoformswitch/cis_overlap/cis_switch_consequences.csv` - Consequence analysis results
- `results/isoformswitch/cis_overlap/cis_switch_summary_with_consequences.csv` - Gene summary with consequences

**Note**: The switch plots show isoform usage comparison between NEPC (right bars) and nonNE (left bars). Positive dIF indicates increased usage in NEPC.

### CIS Data Requirements

The CIS integration requires these files from the fusion_finder analysis:

| File | Description |
|------|-------------|
| `../fusion_finder/R/all_insertions_master.csv` | All CIS insertions with gene annotations |
| `../fusion_finder/R/key_genes_insertions.csv` | Key driver gene insertions (optional) |

These files should contain columns: `gene_name`, `sample`, `rnaseq_id`, `is_nepc`, `classification`, `sb_mechanism`

### Interpretation

The CIS-integrated analysis helps answer:

1. **Which CIS genes show isoform switching?** - Genes where transposon insertion correlates with altered isoform usage
2. **What are the functional consequences?** - Domain loss/gain, NMD sensitivity, coding potential changes
3. **Is the switching sample-specific?** - Do samples with insertions show different patterns than those without?

## Integration with Main Analysis Pipeline

The isoform switching analysis is now integrated into the main paper analysis pipeline at:
`/lab/barcheese01/mdiberna/sleeping-beauty-paper/main.sh`

### Two-Stage Workflow

**Stage 1: Upstream Analysis** (this directory)
- Scripts in `R/` perform the full IsoformSwitchAnalyzeR pipeline
- Outputs saved to `results/isoformswitch/`
- CIS integration scripts (`01_cis_isoform_switch_overlap.R`, `02_cis_switch_consequences.R`) run separately
- All results saved as checkpoints for downstream use

**Stage 2: Publication Outputs** (main paper repo)
- Script: `libs/isoform_switching_analysis.R` in the main paper repo
- Reads pre-computed results from this directory's `results/`
- Generates publication-quality plots and tables
- Outputs to: `experiments/isoform-switching/` in the main paper repo
- Called automatically by `main.sh`

### File Locations

```
isoform_switching/                               # THIS DIRECTORY (source data)
├── R/                                           # Upstream analysis scripts
│   ├── isoform_switch_analysis.R                # Main IsoformSwitchAnalyzeR pipeline
│   ├── 01_cis_isoform_switch_overlap.R          # CIS integration Part 1
│   └── 02_cis_switch_consequences.R             # CIS integration Part 2
└── results/isoformswitch/cis_overlap/           # Pre-computed results (READ by main paper)
    ├── checkpoint_cis_overlap.rds
    ├── pfam_results.txt
    ├── cis_switch_genes_detail.csv
    └── cis_insertions_in_switching_genes.csv

../sleeping-beauty-paper/                        # MAIN PAPER REPO
├── libs/isoform_switching_analysis.R            # Publication script (reads from above)
└── experiments/isoform-switching/               # Publication outputs (WRITE by main paper)
    ├── tables/
    │   ├── cis_isoform_switch_summary.csv
    │   ├── cis_isoform_switch_summary_full.csv
    │   └── *.csv
    └── plots/
        └── switch_*.pdf
```

### Running the Integrated Analysis

From the main paper repository:
```bash
cd /lab/barcheese01/mdiberna/sleeping-beauty-paper
bash main.sh  # Runs all analyses including isoform switching
```

Or run just the publication script:
```bash
Rscript libs/isoform_switching_analysis.R
```

**Note**: The publication script expects pre-computed results in `isoform_switching/results/`.
If you need to re-run the upstream analysis, use the scripts in this directory's `R/` folder.

## References

- [IsoformSwitchAnalyzeR Bioconductor Page](https://bioconductor.org/packages/IsoformSwitchAnalyzeR/)
- [IsoformSwitchAnalyzeR Vignette](docs/IsoformSwitchAnalyzeR_vignette.html)
- [HMMER Web Server](https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan) - For Pfam domain analysis
- Vitting-Seerup K, Sandelin A. IsoformSwitchAnalyzeR: analysis of changes in genome-wide patterns of alternative splicing and its functional consequences. Bioinformatics. 2019.

## Version Information

- R: 4.2.1
- IsoformSwitchAnalyzeR: Check with `packageVersion('IsoformSwitchAnalyzeR')`
- GENCODE: vM25 (mouse)
- HMMER: 3.3.1 (for local Pfam analysis, optional)
