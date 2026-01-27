# Sleeping Beauty Paper

This repository contains all analysis code and data for the paper "A forward genetic screen identifies Sirtuin1 as a driver of neuroendocrine prostate cancer."

## TAPDANCE Analysis (Perl)

The TAPDANCE (Transposon Analysis Pipeline for Driver And Cancer Element identification) implementation for this project is in the `TAPDANCE_NEPC` subdirectory.

Full instructions for the TAPDANCE pipeline, including setup, configuration, and execution, can be found in the [TAPDANCE README](TAPDANCE_NEPC/README.md).

The large input and intermediate files not included in this repository are available at [Zenodo](https://zenodo.org/records/15098998). Final common insertion site level data outputs are contained in the `TAPDANCE_NEPC` subdirectory.

## Downstream Analysis (R)

The main.sh file contains a bash pipeline that will recreate the entire analysis.

To run this pipeline on a Mac, you can do the following:

```bash
bash PATH_TO_DIR/sleeping-beauty-paper/main.sh
```

This will create an `experiments/` directory within `sleeping-beauty-paper/` containing the results.

### Requirements

Note: pandoc is required to generate the network figures. When running the code within RStudio, it is loaded through the application. However, when running it with an RScript (which is what our bash pipeline does), you may encounter the following error:

```
Error in htmlwidgets::saveWidget(graph, file, selfcontained = selfcontained, :
    Saving a widget with selfcontained = TRUE requires pandoc.
```

To address this, include pandoc in your environment's PATH variable. In Mac, this can be done by adding the following lines of code into your shell source (e.g. .zshrc) as shown below:

```bash
PATH="${PATH}:/Applications/RStudio.app/Contents/MacOS/quarto/bin/tools/pandoc"
export PATH
```

After reloading your shell, you should now be able to run the script without encountering any errors.

## RNA-Seq Analyses

Three additional RNA-seq analyses have been integrated into the repository, examining fusion events, isoform switching, and chimeric transcripts in NEPC samples with transposon insertions.

### Fusion Finder Analysis

**Location:** `fusion_finder/`
**Purpose:** Detect RNA fusion events from transposon-genome junctions

The fusion finder analysis identifies RNA-seq reads spanning transposon insertion sites and host genome sequences. It cross-references detected fusion events with CIS (Common Insertion Site) data to identify genes where transposon insertions correlate with aberrant transcription.

**Key scripts:**
- `libs/fusion-finder/overlap_cis_fusions.R` - Cross-reference CIS data with fusion events
- `libs/fusion-finder/figure_S4BC_fusions.R` - Generate Figure S4B-C (fusion counts by tumor type and transposon region)
- `libs/fusion-finder/figure_S4B_fusions.R` - Generate Figure S4B (box plot)

**Preprocessing:** Alignment and fusion detection scripts are in `fusion_finder/scripts/`

**Outputs:**
- `fusion_finder/results/cis_fusion_overlap/` - CIS-fusion overlap tables
- `fusion_finder/results/cis_fusion_overlap/figures/` - Publication figures

### Isoform Switching Analysis

**Location:** `isoform_switching/`
**Purpose:** Identify isoform switches in genes with transposon insertions

This analysis uses IsoformSwitchAnalyzeR to detect differential isoform usage between NEPC and non-NEPC samples, focusing on genes with CIS insertions. It identifies cases where transposon insertions may drive alternative splicing or transcript isoform switches.

**Key scripts:**
- `libs/isoform-switching/isoform_switch_analysis.R` - Main isoform switching detection
- `libs/isoform-switching/01_cis_isoform_switch_overlap.R` - Identify CIS genes with isoform switches
- `libs/isoform-switching/02_cis_switch_consequences.R` - Functional consequence analysis with Pfam domains
- `libs/isoform-switching/plot_sirt1_coverage.R` - SIRT1 locus coverage analysis

**Preprocessing:** Salmon quantification and alignment scripts in `isoform_switching/scripts/`

**Outputs:**
- `isoform_switching/results/isoformswitch/` - Switch analysis results
- `isoform_switching/results/isoformswitch/cis_overlap/` - CIS-isoform overlap tables
- `isoform_switching/results/isoformswitch/plots/` - Switch plots

**Note:** The isoform switching analysis requires Pfam domain annotation. After running `01_cis_isoform_switch_overlap.R`, upload the generated FASTA file to the HMMER webserver (https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan) and then run `02_cis_switch_consequences.R` with the Pfam results.

### Chimeric Isoforms Analysis

**Location:** `chimeric_isoforms/`
**Purpose:** Analyze RNA-seq coverage at CIS insertion sites to detect chimeric transcripts

This analysis examines read coverage across gene bodies in samples with transposon insertions, comparing NEPC samples to control groups (SB-negative and SB-positive non-NEPC). It identifies regions where transposon-driven transcription may create chimeric or aberrant isoforms.

**Key scripts:**
- `libs/chimeric-isoforms/00_prepare_data.R` - Filter CIS data and prepare sample groups
- `libs/chimeric-isoforms/01_coverage_analysis.R` - Calculate RNA-seq coverage at insertion sites
- `libs/chimeric-isoforms/02_publication_plots.R` - Generate publication-quality coverage plots
- `libs/chimeric-isoforms/02_trackviewer_plots.R` - Alternative plotting with trackViewer

**Outputs:**
- `chimeric_isoforms/results/coverage_data/` - Coverage data by gene
- `chimeric_isoforms/results/plots/` - Coverage visualization plots
- `chimeric_isoforms/results/tables/` - Summary statistics

### Shared CIS Data

Common Insertion Site (CIS) results used across analyses are stored in:
- `data/cis-results/all_insertions_master.csv` - All CIS insertions with gene annotations
- `data/cis-results/key_genes_insertions.csv` - Key gene CIS insertions

### Running RNA-Seq Analyses

The RNA-seq analyses are integrated into `main.sh` and will run as part of the full pipeline. To run individual analyses:

```bash
# Fusion finder
Rscript libs/fusion-finder/overlap_cis_fusions.R

# Isoform switching
Rscript libs/isoform-switching/01_cis_isoform_switch_overlap.R

# Chimeric isoforms
Rscript libs/chimeric-isoforms/00_prepare_data.R
Rscript libs/chimeric-isoforms/01_coverage_analysis.R
Rscript libs/chimeric-isoforms/02_publication_plots.R
```

**Note:** These analyses require preprocessed data (aligned BAM files, Salmon quantifications) which are not included in this repository due to size constraints. Preprocessing scripts are available in each analysis directory under `scripts/`.
