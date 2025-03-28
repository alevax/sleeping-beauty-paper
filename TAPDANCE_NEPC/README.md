# TAPDANCE_NEPC

Transposon Analysis Pipeline for Driver And Cancer Element identification, modified for Neuroendocrine Prostate Cancer Sleeping Beauty analysis. This version is adapted from Tim Starr's TAPDANCE [implementation](https://sourceforge.net/p/tapdancebio/home/Home/).

## Pre-Analyzed Data

This repository includes the complete TAPDANCE analysis results for the NEPC study. All output files and data are provided, with two exceptions:

1. **`seqs.tab`**: This large input file is not included in the repository but is available at [Zenodo](https://zenodo.org/records/15098998) 

2. **`mapping/` folder**: The intermediate mapping files produced by TAPDANCE.pl are also available at the same [Zenodo](https://zenodo.org/records/15098998) link

Additional analysis and results beyond the TAPDANCE pipeline are available in a companion repository: [https://github.com/mat10d/TAPDANCE_NEPC](https://github.com/mat10d/TAPDANCE_NEPC)

## Setup Instructions

### Prerequisites
- Perl 5 with cpanm
- Homebrew (for macOS users)
- Bowtie

### MySQL Setup
1. Install MySQL:
   ```bash
   brew install mysql
   ```

2. Configure MySQL:
   ```bash
   # Set root password
   mysqladmin -u root password 'YOUR_SECURE_PASSWORD'

   # Log in to MySQL
   mysql -u root -p
   
   # Configure MySQL settings
   SET GLOBAL sql_mode=(SELECT REPLACE(@@sql_mode,'ONLY_FULL_GROUP_BY',''));
   SET GLOBAL local_infile = 'ON';
   
   # Create database
   CREATE DATABASE tapdance;
   QUIT;
   ```

### Perl Dependencies
Install required Perl modules:
```bash
cpanm DBD::mysql
cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
cpan DBD::mysql
```

### Bowtie Setup
1. Install [Bowtie](https://bowtie-bio.sourceforge.net/index.shtml). This TAPDANCE analysis used Bowtie 1.2.3. 
2. Add the mm10 genome folder to the `indexes` subfolder of your Bowtie installation
3. Update your PATH to include Bowtie:
   ```bash
   PATH=$PATH:/path/to/your/bowtie
   ```

### Annotation Setup
The pipeline requires a genomic annotation file in BED format (`lib/mm10.bed`) that contains gene annotations.

#### Obtaining the mm10.bed file:
The annotation file can be generated using [UCSC's Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables?command=start):
- **Assembly**: Dec. 2011 (GRCm38/mm10)
- **Track**: NCBI RefSeq
- **Table**: UCSC RefSeq (refGene)
- **Output format**: BED

Example of mm10.bed format:
```
chr1	3214481	3671498	Xkr4	0	-	3216021	3671348	0	3	2487,200,947,	0,207220,456070,
chr1	4290845	4409241	Rp1	0	-	4292980	4409187	0	4	2167,172,636,72,	0,61064,61356,118324,
```

The annotation version should match the Bowtie genome version. The mm10.bed file used for this TAPDANCE implementation is included in the repository under `lib/mm10.bed`.

## Configuration

### Modifying `config.pl`
Before running the pipeline, you must modify the `config.pl` file to match your environment and experiment. This is a **critical step** that requires several adjustments:

1. **Database Settings**
   ```perl
   $db_name     = 'tapdance_database';  # Change to your database name
   $db_username = 'tapdance_username';  # Your MySQL username
   $db_password = 'tapdance_password';  # Change to your MySQL password
   $db_host     = 'localhost;mysql_local_infile=1';
   ```

2. **Project Settings**
   ```perl
   $proj = 'prostate';                 # Project name (used as suffix for all tables)
   $annotation_file = 'lib/mm10.bed';  # Path to genome annotation file
   $bowtie_path = '/path/to/your/bowtie/indexes/mm10/mm10';  # Update with your Bowtie path
   ```

3. **Analysis Parameters**
   ```perl
   $library_percent = '0.01';     # Threshold for including insertions (% of total reads)
   $CIS_total_pvalue = '0.05';    # p-value threshold for CIS total count
   $CIS_library_pvalue = '0.05';  # p-value threshold for CIS library count
   $CIS_region_pvalue = '0.05';   # p-value threshold for CIS region count
   $cocis_threshold = '0.001';    # p-value threshold for co-occurring CIS
   ```

4. **Barcode and IRDR Processing**
   
   The `resolve_barcodes` and `resolve_IRDR` subroutines must be adjusted based on your experiment:
   
   ```perl
   sub resolve_barcodes {
     # Adjust the number (13) to match your barcode length + 1
     $sth = $dbh->prepare("create table illumina_decoded_$proj select library, id, substring(sequence,13) as decoded_sequence from barcode_$proj,illumina_raw_$proj where sequence like concat(seq,'%')");
     $sth->execute;
   }
   
   sub resolve_IRDR {
     # Adjust the sequence pattern and substring length (30) to match your IRDR sequence
     $sth = $dbh->prepare("create table illumina_without_IRDR_$proj select library,id,substring(decoded_sequence,30) as insertion_sequence, 'good' as type from illumina_decoded_$proj where decoded_sequence like '___TGTATGTAAACTTCCGACTTCAACTG%'");
     $sth->execute;
   }
   ```
   
   **Important notes for adjusting these values:**
   - For `resolve_barcodes`: The number (13) should be your barcode length + 1
   - For `resolve_IRDR`: The sequence pattern should match your transposon's IRDR sequence
   - The substring length (30) must match the length of the search string, including wildcards (`_`) and the `%` sign

## Input Data Files

The pipeline requires several input files to be placed in the `data` directory:

### 1. `seqs.tab`
Contains raw sequencing reads that include barcodes and IRDR (inverted repeat/direct repeat) sequences.

Example:
```
AACGGACAACTT            AACGGACAACTTAAGTGTATGTAAACTTCCGACTTCAACTGTATATAAAGAACTCAAGAAGGTAGACTCCAGAAAATCAAATAACCCCATTAAAAAATGGG
AACGGACAACTT            AACGGACAACTTAAGTGTATGTAAACTTCCGACTTCAACTGTAGCTCAGCTGGGATGTGAACAAAAGTTTCCGGGATTGTGTGTTACTTCCATTCAGTATG
AACGGACAACTT            AACGGACAACTTAAGTGTATGTAAACTTCCGACTTCAACTGTATACTTAAGTGGTGTTAAGTATATGCACACTGTTGTACAACCATTTCTACCATTTTCTC
```

**Used in**: `TAPDANCE.pl` - Primary input file for identifying and mapping transposon insertions.

**Note**: Due to its large size, this file is not included in the repository but is available at the [Zenodo](https://zenodo.org/records/15098998) link.

### 2. `barcode2lib.txt`
Maps barcodes to specific libraries and indicates whether they're from the left or right end of the transposon.

Example:
```
AACGGACAACTT    B-A01-L    Left
AACTACGACTGT    N-A01-R    Right
AACTCTAATGGC    B-A02-L    Left
```

**Used in**: `TAPDANCE.pl` - Used for identifying the source (library) of each sequence in seqs.tab.

### 3. `chromo.tab`
Lists chromosomes to be excluded from the analysis.

Example (to exclude chromosome 1):
```
chr1
```

**Used in**: `TAP2.pl` - Specifies which chromosomes to exclude when identifying common insertion sites.

### 4. `metadata.tab`
Defines groups of libraries for comparison and CIS analysis.

Example:
```
B3-1    metastatic    cis
B3-4    metastatic    cis
B3-8    metastatic    cis
```

**Used in**: `TAP2.pl` - Groups libraries for the identification of common insertion sites (CIS).

## Running TAPDANCE

### Execution
Navigate to the TAPDANCE_NEPC directory:
```bash
cd TAPDANCE_NEPC
```

#### Step 1: Map fragments to reads
```bash
perl lib/TAPDANCE.pl
```
This step:
- Processes raw sequence data from `seqs.tab`
- Identifies barcode sequences using `barcode2lib.txt`
- Removes IRDR sequences and maps to the genome using Bowtie
- Creates tables of mapped insertions

#### Step 2: Identify Common Insertion Sites (CIS)
```bash
perl lib/TAP2.pl
```
This step:
- Uses the libraries defined in `metadata.tab`
- Excludes chromosomes listed in `chromo.tab`
- Identifies genomic regions with statistically significant insertion densities
- Associates insertions with nearby genes

**Important**: You can run `TAP2.pl` multiple times with different settings without having to rerun `TAPDANCE.pl`. This allows you to:
- Try different `$library_percent` thresholds in `config.pl`
- Modify `metadata.tab` to analyze different groups of libraries
- Adjust p-value thresholds for CIS detection
- Test different chromosome exclusions in `chromo.tab`

After changing any of these settings, simply run `TAP2.pl` again to generate a new analysis.

### Results
Analysis results will be stored in `results/[metadata_descriptor]`, where `metadata_descriptor` is specified in the metadata.tab file.

**Note**: This repository already contains the complete analysis results from running TAPDANCE on the NEPC dataset. The R code for higher-level analysis is not included here but is available in the companion repository [https://github.com/mat10d/TAPDANCE_NEPC](https://github.com/mat10d/TAPDANCE_NEPC).

## Reference
If you use this pipeline in your research, please cite the original TAPDANCE publication:
- Sarver AL, et al. (2012). TAPDANCE: An automated tool to identify and annotate transposon insertion CISs and associations between CISs from next generation sequence data. BMC Bioinformatics, 13:154.