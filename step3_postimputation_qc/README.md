# STEP 3 - Post-Imputation Quality Control Pipeline

This Bash pipeline performs post-imputation quality control on genotyped data. It prepares imputed variants for downstream polygenic risk score (PRS) analysis by removing ambiguous or duplicated variants, applying minor allele frequency and missingness filters, and generating diagnostic plots.

The main script is `03_POSTIMPUTATION_DATA_PIPELINE_PRC_v4.0.sh`, version 4.0, authored by Pol Ramon Cañellas.

## Pipeline overview

| Step | Main action |
|--- | ---|
| 1  | Downloads Michigan Imputation Server results, extracts chromosomes 1–22, and concatenates them into one VCF. |
| 2  | Filters variants by imputation quality (`INFO/R2 > 0.9`) and converts the VCF to PLINK binary format. |
| 3  | Removes palindromic variants, indels, duplicated positions, and duplicated IDs; assigns missing variant IDs. |
| 4  | Plots the MAF distribution and removes variants below the configured MAF threshold. |
| 5  | Plots variant/sample missingness and applies two consecutive SNP missingness filters. |

The pipeline therefore runs the complete workflow, from the imputation-server results through the final post-imputation PLINK dataset.

## Requirements

- Linux or an HPC environment with Bash.
- GNU command-line utilities, including `awk`, `cut`, `find`, `sort`, and `uniq`.
- [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/).
- R 4.x.
- BCFtools.
- `curl` and `unzip` for result download and extraction.
- The R packages imported by the plotting scripts.

The original HPC configuration loads:

```bash
module load BCFtools/1.21-GCC-12.3.0
module load PLINK/1.9b_6.21-x86_64
module load R/4.3.2-gfbf-2023a
```

On systems without environment modules, the corresponding executables must be available in `PATH`.

## Input files

Run the pipeline from a dedicated directory containing the main script and auxiliary R scripts:

```text
03_POSTIMPUTATION_DATA_PIPELINE_PRC_v4.0.sh
POS4_MAF_check.R
POS5_Hist_miss.R
```

### Imputation-server results

Step 1 expects password-protected chromosome archives from the Michigan Imputation Server:

```text
chr_1.zip
chr_2.zip
...
chr_22.zip
```

After extraction, the script expects:

```text
chr1.dose.vcf.gz
chr2.dose.vcf.gz
...
chr22.dose.vcf.gz
```

These files are concatenated into `1_AllChromosomes.vcf.gz`, filtered in step 2, and converted to PLINK format.

Configure the download URL and archive password without committing either secret to version control:

```bash
impdata="<Michigan Imputation Server download URL>"
imppassword="<archive password>"
```

There is a naming mismatch in the supplied commands: step 2 writes `Target.imputated.*`, while step 3 reads `2_TOPMED_Target.imputated.*`. These prefixes must be made identical for the full pipeline to continue successfully.

## Configuration

Set the PLINK executable and filtering thresholds near the beginning of the script:

```bash
runplink="plink"

s4mymaf=0.001
s5mymissingness1=0.2
s5mymissingness2=0.02
```

| Variable | Default | Purpose |
|---|---:|---|
| `s4mymaf` | 0.001 | Minimum minor allele frequency retained after imputation. |
| `s5mymissingness1` | 0.20 | Initial maximum missing genotype rate per SNP. |
| `s5mymissingness2` | 0.02 | Final, more stringent maximum missing genotype rate per SNP. |

Thresholds should be chosen according to sample size, imputation panel, downstream analysis, and study protocol.

## Usage

From the directory containing all required files:

```bash
chmod +x 03_POSTIMPUTATION_DATA_PIPELINE_PRC_v4.0.sh
bash 03_POSTIMPUTATION_DATA_PIPELINE_PRC_v4.0.sh
```

The parameters for the current run are written to:

```text
03_POSTIMPUTATION_DATA_PIPELINE_PRC_v4.0.log
```

## Processing details

### Step 1: download and merge imputed data

The pipeline downloads the Michigan Imputation Server results, extracts the 22 password-protected chromosome archives, and concatenates `chr1.dose.vcf.gz` through `chr22.dose.vcf.gz` into:

```text
1_AllChromosomes.vcf.gz
```

The chromosome files must be ordered numerically, as shown in the script, to produce a correctly ordered combined VCF.

### Step 2: imputation-quality filtering

BCFtools retains variants satisfying:

```text
INFO/R2 > 0.9
```

PLINK then converts the filtered VCF into a binary dataset using `--const-fid`. The large filtered intermediate VCF is removed after conversion.

### Step 3: data adaptation

The script performs the following operations in order:

1. Removes palindromic A/T, T/A, C/G, and G/C variants.
2. Removes indels and alleles outside A/C/G/T using `--snps-only just-acgt`.
3. Assigns coordinate- and allele-based IDs to variants with missing identifiers.
4. Detects and removes duplicated variant positions.
5. Detects and removes duplicated variant IDs.

The working dataset prefix becomes `TOPMED_Target.imputated`.

### Step 4: minor allele frequency

PLINK calculates allele frequencies before and after filtering. `POS4_MAF_check.R` generates the corresponding MAF distribution plots. Variants with a MAF below `s4mymaf` are removed.

The resulting working prefix is:

```text
TOPMED_Target.imputated_maf
```

### Step 5: SNP missingness

PLINK calculates missingness per individual and per variant. `POS5_Hist_miss.R` generates diagnostic histograms. Two consecutive `--geno` filters are then applied:

1. Missingness greater than 0.20 is removed.
2. Missingness greater than 0.02 is removed from the remaining variants.

The script only filters SNP missingness in this stage; it does not remove individuals with `--mind`.

## Output files

The final post-imputation PLINK dataset is:

```text
TOPMED_Target.POSTimputed.bed
TOPMED_Target.POSTimputed.bim
TOPMED_Target.POSTimputed.fam
```

Intermediate reports, exclusion lists, PLINK logs, and plots are moved to:

```text
FilesCreated/
```

Examples include MAF distributions, missingness histograms, duplicate lists, and stage-specific PLINK logs.

## Important warnings

- Always work on copies of the original data. PLINK repeatedly overwrites the working dataset prefixes.
- Step 2 generates `Target.imputated.*`, while step 3 expects `2_TOPMED_Target.imputated.*`. Align these prefixes before running the full workflow; otherwise, step 3 will fail even when steps 1 and 2 complete successfully.
- Verify that `INFO/R2` is the correct imputation-quality field for the selected server and reference panel before enabling the step 2 filter.
- Do not store an imputation-server download URL or password in a tracked script. Supply secrets securely at runtime and rotate any credentials that may have been exposed.
- `FilesCreated/` is deleted recursively without confirmation before being recreated.
- Files ending in `~` are deleted from the working directory.
- Broad patterns (`1*` through `5*`) are moved into `FilesCreated/`. Use a dedicated directory containing no unrelated files that begin with those digits.
- The script does not enable `set -e`, so it may continue after a command fails. Check the terminal output and confirm every expected output file before downstream use.
- Converting dosage VCFs to PLINK 1 binary files may turn probabilistic imputed genotypes into hard calls. Confirm PLINK's genotype-call behavior and any required dosage threshold for the intended downstream analysis.

## Recommended directory structure

For the complete workflow:

```text
step3_postimputation_qc/
├── 03_POSTIMPUTATION_DATA_PIPELINE_PRC_v4.0.sh
├── chr_1.zip
├── chr_2.zip
├── ...
├── chr_22.zip
├── POS4_MAF_check.R
└── POS5_Hist_miss.R
```

After a successful run, the final `TOPMED_Target.POSTimputed.*` dataset, the parameter log, and `FilesCreated/` directory will also be present.

## Limitations

This pipeline assumes autosomal data already imputed through a compatible workflow and is tailored to PLINK 1.9 binary output. It does not perform sample-level post-imputation QC, ancestry reassessment, phenotype validation, or association testing. Diagnostic plots and exclusion counts should be reviewed before the final dataset is used for PRS or GWAS analyses.

