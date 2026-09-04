# Pre-GWAS Quality Control with PLINK

This Bash pipeline performs quality control (QC) on genetic data before a genome-wide association study (GWAS). It applies consecutive variant- and sample-level filters and generates diagnostic plots with R scripts.

The main script is `01_QC_TARGET_DATA_PIPELINE_PRC_v4.1.sh`, version 4.1, authored by Pol Ramon Cañellas.

## Pipeline overview

| Step | Analysis | Main action |
|---|---|---|
| 0 | Data preparation | Converts PED/MAP files to binary format, optionally updates phenotypes, and removes indels and non-ACGT alleles. |
| 1 | Missingness | Examines and filters variants and individuals according to missing genotype rates. |
| 2 | Sex discrepancy | Compares recorded sex with X-chromosome estimates and removes discrepancies. |
| 3 | Minor allele frequency | Keeps autosomal SNPs and removes variants with a low MAF. |
| 4 | Hardy-Weinberg equilibrium | Filters variants by HWE in controls and then across the complete dataset. |
| 5 | Heterozygosity | Performs linkage disequilibrium pruning and removes individuals more than three standard deviations from the mean. |
| 6 | Relatedness | Evaluates identity by descent, keeps founders, and removes individuals above the relatedness threshold. |

## Requirements

- Linux or an HPC environment with Bash and GNU utilities (`awk`, `sed`, `grep`, `find`, `sort`, and `uniq`).
- [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/).
- R 4.x.
- The R packages required by the auxiliary scripts.
- The `rename` command if required by any auxiliary script.

The script was prepared for an HPC environment using these modules:

```bash
module load PLINK/1.9b_6.21-x86_64
module load R/4.3.2-gfbf-2023a
```

On systems without environment modules, PLINK and `Rscript` must be available in `PATH`. Alternatively, set `runplink` to the absolute path of the PLINK executable.

## Input files

All input and auxiliary files must be in the directory from which the pipeline is run.

### Genetic data

The code uses the fixed prefix `Target_data` and therefore expects:

```text
Target_data.ped
Target_data.map
```

The initial conversion generates and uses:

```text
Target_data.bed
Target_data.bim
Target_data.fam
```

If the data are already in binary format, they must use this same prefix. Although the original comments mention names such as `Target_data_*****`, the script does not detect suffixes automatically.

### Auxiliary scripts and files

```text
QC0_Target_data_phenotype.R          # Only when phenotype adaptation is needed
QC1_Hist_miss.R
QC2_Gender_check.R
QC3_MAF_check.R
QC4_HWE.R
QC5_Check_heterozygosity_rate.R
QC5_Inversion.txt                    # High-LD regions
QC6_Relatedness.R
```

Additional requirements:

- `Target_data_pheno_sex.txt` is required when step 2 is enabled. It must contain the columns `FID IID SEX`, without a header, using PLINK's accepted sex coding.
- If `step0fam="YES"`, the phenotype file expected by `QC0_Target_data_phenotype.R` is also required. That R script must first be adapted to the input data and must generate `0A_phenotypes.txt` containing `FID IID PHENO`.
- PLINK normally codes phenotypes as `1` for controls, `2` for cases, and `-9` for missing values.

## Configuration

Edit the relevant variables near the beginning of the main script:

```bash
DATA="Target_data"
runplink="plink"

step0fam="NO"
step1="YES"
step2="YES"
step3="YES"
step4="YES"
step5="YES"
step6="YES"
```

Use uppercase `YES` to enable a step and `NO` to skip it.

The default thresholds are:

| Variable | Value | Purpose |
|---|---:|---|
| `mymissigness1` | 0.20 | Initial missingness filter for SNPs and individuals. |
| `mymissigness2` | 0.02 | More stringent missingness filter. |
| `mymaf` | 0.01 | Minimum minor allele frequency. |
| `myhwe1` | 1e-6 | HWE threshold initially applied by PLINK to controls. |
| `myhwe2` | 1e-10 | HWE threshold applied with `--hwe-all`. |
| `mywindowsize` | 50 | Window size for LD pruning. |
| `myshiftwindow` | 5 | Number of SNPs by which the pruning window is shifted. |
| `mypairwiser2` | 0.20 | Correlation threshold for `--indep-pairwise`. |
| `myrelatedness` | 0.125 | Minimum PI_HAT threshold for detecting relatedness. |

These thresholds should be reviewed according to the sample size, population, and study design.

## Usage

Change to the directory containing the input data, main script, and auxiliary files, then run:

```bash
chmod +x 01_QC_TARGET_DATA_PIPELINE_PRC_v4.1.sh
bash 01_QC_TARGET_DATA_PIPELINE_PRC_v4.1.sh
```

The pipeline repeatedly updates `Target_data.bed`, `Target_data.bim`, and `Target_data.fam`, so every filter is applied to the result of the preceding filter.

## Output

The final quality-controlled dataset is saved as:

```text
QCtargetdata.bed
QCtargetdata.bim
QCtargetdata.fam
```

The parameters used for the current run are recorded in:

```text
01_QC_TARGET_DATA_PIPELINE_PRC_v4.1.log
```

Intermediate reports, exclusion lists, PLINK files, and plots are renamed with their step prefix (`0A_`, `1B_`, `3A_`, and so on). At the end of the run, they are moved into `FilesCreated/`.

## Dependencies between steps

- Step 6 uses `4B_indepSNP.prune.in`, which is generated during step 5. Do not run step 6 while skipping step 5 unless a compatible copy of this file already exists.
- Later steps operate on the cumulative output of earlier steps. Skipping a step may change how subsequent checks should be interpreted.
- Step 2 invokes `plink` directly once instead of using `$runplink`. Therefore, `plink` must be available in `PATH` even when `runplink` points to another location.

## Important warnings

- Always work on a copy of the original data. The pipeline overwrites the `Target_data.*` files throughout the analysis.
- At the end of the pipeline, `rm -r FilesCreated` is run before the directory is recreated. If the directory already exists, all its contents are deleted without confirmation.
- Files ending in `~` are also deleted from the working directory.
- The patterns used to move results (`0*`, `1*`, ..., `6*`) are broad. Run the pipeline in a dedicated directory containing no unrelated files with those prefixes.
- The relatedness stage removes every unique identifier appearing as the first individual in each pair in the `.genome` file. Review `6B_Target_data.genome` and `6B_Target_data.txt` before interpreting the final dataset.
- The script does not stop automatically when a command fails. Review the terminal output and generated logs before using the results.

## Recommended directory structure

```text
project_qc/
├── 01_QC_TARGET_DATA_PIPELINE_PRC_v4.1.sh
├── Target_data.ped
├── Target_data.map
├── Target_data_pheno_sex.txt
├── QC0_Target_data_phenotype.R
├── QC1_Hist_miss.R
├── QC2_Gender_check.R
├── QC3_MAF_check.R
├── QC4_HWE.R
├── QC5_Check_heterozygosity_rate.R
├── QC5_Inversion.txt
└── QC6_Relatedness.R
```

After a successful run, the `QCtargetdata.*` dataset, parameter log, and `FilesCreated/` directory will also be present.

## Limitations

This pipeline is designed for PLINK 1.9, biallelic data, and a case-control study design. It does not replace manual review of diagnostic plots, exclusion lists, or cohort-specific criteria. Further checks may be required before a GWAS, including ancestry, duplicates, allele concordance, batch effects, and comparison against a reference panel.
