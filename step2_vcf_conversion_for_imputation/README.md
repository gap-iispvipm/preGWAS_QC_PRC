# STEP 2 - Pre-Imputation Pipeline for Genotyped Data

This Bash pipeline prepares quality-controlled genotype data for imputation. It harmonizes a target dataset with a reference dataset, performs principal component analysis (PCA), checks and updates alleles against the HRC reference, and creates compressed, indexed VCF files for chromosomes 1–22.

The main script is `02_PREIMPUTATION_DATA_PIPELINE_PRC_v5.0.sh`, version 5.0, authored by Pol Ramon Cañellas.

## Pipeline overview

| Step | Analysis | Main action |
|---|---|---|
| 0 | Dataset preparation | Copies input PLINK files, removes ambiguous and duplicate variants, assigns coordinate-based IDs, and optionally performs LiftOver. |
| 1 | Dataset harmonization | Sets reference alleles, detects strand differences, flips target variants, and removes unresolved variants. |
| 2 | Dataset merge | Merges target and reference datasets for population structure analysis. |
| 3 | Merged PCA | Performs LD pruning and PCA on the merged target/reference dataset. |
| 4 | Target PCA | Calculates the first 10 principal components for use as covariates. |
| 5 | HRC allele check | Downloads the HRC checking utility and reference sites, then generates and runs allele-correction commands. |
| 6 | VCF preparation | Sorts and indexes one VCF file per autosomal chromosome. |

## Requirements

- Linux or an HPC environment with Bash.
- GNU command-line utilities, including `awk`, `cut`, `find`, `sort`, `uniq`, `wget`, `unzip`, and `gunzip`.
- [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/).
- R 4.x.
- Perl with BioPerl.
- BCFtools.
- Internet access during step 5.
- UCSC `liftOver` and an appropriate chain file if genome-build conversion is enabled.

The original HPC configuration loads:

```bash
module load PLINK/1.9b_6.21-x86_64
module load R/4.3.2-gfbf-2023a
module load BioPerl/1.7.8-GCCcore-12.3.0
module load BCFtools/1.21-GCC-12.3.0
```

On systems without environment modules, the corresponding executables must be available in `PATH`.

## Input files

Run the pipeline from a dedicated directory containing the main script, both datasets, and all auxiliary files.

### Target dataset

```text
QCtargetdata_<name>.bed
QCtargetdata_<name>.bim
QCtargetdata_<name>.fam
```

### Reference dataset

```text
QCreferencedata_<name>.bed
QCreferencedata_<name>.bim
QCreferencedata_<name>.fam
```

The wildcard patterns must match exactly one BED/BIM/FAM triplet for each dataset. At the start of the run, the files are copied to the fixed prefixes `QCtargetdata` and `QCreferencedata`.

Both datasets must use the same reference genome build, either GRCh37/hg19 or GRCh38/hg38. If they do not, configure LiftOver for one dataset before running the pipeline.

### Auxiliary R scripts

```text
PRE0_FID_IID.R
PRE0D_Duplicates.R
PRE3_PCA.R
```

Their expected roles are:

- `PRE0_FID_IID.R`: creates the ID update file used after LiftOver.
- `PRE0D_Duplicates.R`: chooses duplicate variants to exclude.
- `PRE3_PCA.R`: plots and interprets PCA results and may generate an ancestry-outlier list.

Install all R packages imported by these scripts before starting the pipeline.

## Configuration

Edit these variables near the beginning of the script:

```bash
runplink="plink"

s0liftoverTarget="NO"
s0liftoverReference="NO"

mywindowsize=50
myshiftwindow=5
mypairwiser2=0.2
```

If LiftOver is enabled, uncomment and configure:

```bash
runliftover="/path/to/liftOver"
liftoverChain="/path/to/hg38ToHg19.over.chain"
```

The chain file must represent the intended source-to-destination conversion. Do not infer the conversion direction from the example filename without first confirming the builds of both input datasets.

Only one LiftOver branch can run: the script uses `if ... elif`, so if both variables are set to `YES`, only the target dataset is lifted.

## Usage

From the directory containing all input files:

```bash
chmod +x 02_PREIMPUTATION_DATA_PIPELINE_PRC_v5.0.sh
bash 02_PREIMPUTATION_DATA_PIPELINE_PRC_v5.0.sh
```

The parameters used for the run are written to:

```text
02_PREIMPUTATION_DATA_PIPELINE_PRC_v5.0.log
```

## Processing details

### Variant cleanup

The pipeline:

1. Removes palindromic A/T, T/A, C/G, and G/C variants from the target dataset.
2. Assigns IDs to variants with missing identifiers.
3. Removes duplicated variant IDs.
4. Uses frequency information and `PRE0D_Duplicates.R` to choose duplicate positions for removal.
5. Renames variants as `CHR:BP`.
6. Removes variants that still share the same coordinate-based ID.

### Target/reference harmonization

The reference allele is taken from column 5 of `QCreferencedata.bim`. Target variants with apparent strand differences are flipped. Variants that remain inconsistent are excluded from both datasets before merging.

### Principal component analysis

The merged dataset is LD-pruned with:

```text
--indep-pairwise 50 5 0.2
```

PCA is then run on the pruned merged dataset. A second PCA calculates 10 components from the target dataset, producing:

```text
PREtarget_PCA.results.covariate
```

This file can be used as a covariate input in downstream analyses such as PRSice.

### HRC preparation

Step 5 downloads and runs `HRC-1000G-check-bim.pl` with the GRCh37 HRC r1.1 sites file. The utility creates `Run-plink.sh`, which the pipeline executes to update variant names, strands, and reference alleles and to export chromosome-level VCF files.

## Output files

The updated PLINK target dataset is renamed to:

```text
PREIMPUTEDtargetdata.bed
PREIMPUTEDtargetdata.bim
PREIMPUTEDtargetdata.fam
```

The target PCA covariate file is:

```text
PREtarget_PCA.results.covariate
```

The chromosome-level files expected from steps 5 and 6 are:

```text
QCtargetdata-updated-chr1.vcf.gz
QCtargetdata-updated-chr1.vcf.gz.tbi
...
QCtargetdata-updated-chr22.vcf.gz
QCtargetdata-updated-chr22.vcf.gz.tbi
```

Intermediate reports, plots, and lists are moved to `FilesCreated/`. Intermediate BED/BIM/FAM files are subsequently deleted from that directory to reduce disk usage.

## Important warnings

- Work on copies of the original data. The pipeline repeatedly overwrites working datasets and deletes several intermediate files.
- `FilesCreated/` is deleted recursively without confirmation before being recreated.
- Files ending in `~` are deleted from the working directory.
- Broad patterns such as `0*`, `1*`, ..., `6*` are moved into `FilesCreated/`. Use a dedicated directory with no unrelated files matching those patterns.
- Step 0D explicitly reads `QCtargetdata.bed~`, `QCtargetdata.bim~`, and `QCtargetdata.fam~`. These backup files are not guaranteed to be generated by every PLINK installation. Verify that they exist and represent the intended dataset before running this step.
- Step 1 uses `0F_QCtargetdata` as input. Confirm that this dataset was created successfully and, when LiftOver is enabled, that the intended post-LiftOver dataset is actually used.
- The LiftOver branches are marked as untested in the source script and require validation before production use.
- The script does not enable `set -e`, so it may continue after a failed command. Check command output and confirm that every expected file exists.
- Downloads use legacy HTTP and FTP endpoints. Their availability and integrity should be verified before use.
- The script comments refer to both TOPMed and HRC, but the implemented check downloads the HRC r1.1 GRCh37 reference and recommends the HRC imputation panel. This workflow should not be described as a TOPMed preparation pipeline without changing and validating the reference-specific steps.
- The introductory comment claims that the final files are named `2_PreImputed.split.XX.vcf.gz`, but the implemented commands operate on `QCtargetdata-updated-chrXX.vcf.gz`.

## Imputation-server settings stated by the script

The final messages recommend uploading the chromosome-level VCF files to the Michigan Imputation Server with:

| Setting | Value |
|---|---|
| Reference panel | HRC r1.1 2016 |
| Array build | GRCh37/hg19 |
| R² filter | 0.3 |
| Phasing | Eagle v2.4, phased output |
| Population | EUR |
| Mode | Quality Control & Imputation |

These values are study-specific defaults from the script, not universal recommendations. Confirm that the selected panel, population, genome build, and server options are appropriate for the cohort and are still supported before submitting data.

## Recommended directory structure

```text
preimputation_project/
├── 02_PREIMPUTATION_DATA_PIPELINE_PRC_v5.0.sh
├── QCtargetdata_<name>.bed
├── QCtargetdata_<name>.bim
├── QCtargetdata_<name>.fam
├── QCreferencedata_<name>.bed
├── QCreferencedata_<name>.bim
├── QCreferencedata_<name>.fam
├── PRE0_FID_IID.R
├── PRE0D_Duplicates.R
└── PRE3_PCA.R
```

## Limitations

This pipeline is written for PLINK 1.9 and autosomal, biallelic genotype data. Population assignment, PCA outlier removal, genome-build conversion, allele harmonization, and reference-panel selection require cohort-specific review. Always inspect the PLINK logs, PCA plots, exclusion lists, HRC check report, and generated VCF files before submitting data for imputation.
