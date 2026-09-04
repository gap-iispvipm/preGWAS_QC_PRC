###################################################################
######################### START ANALISIS ##########################
###################################################################
###### Script: 03 POSTIMPUTATION PIPELINE FOR GENOTYPED DATA ######
################### Author: POL RAMON CAÑELLAS ####################
########################## Version: v4.0 ##########################
###################################################################


###################################################################
### IMPORTANT INSTRUCTIONS TO RUN THE PIPELINE OF THIS SCRIPT ###

# Run the present script into the folder containing all files from your data.

# FOLDER MUST CONTAINING FILES:
	#### 03_POSTIMPUTATION_DATA_PIPELINE_PRC_v4.0.sh
	#### POS4_MAF_check.R
	#### POS5_Hist_miss.R
	#### chr_(1-22).zip

# ===================================
# MODULES NEEDED
module load BCFtools/1.21-GCC-12.3.0
module load PLINK/1.9b_6.21-x86_64
module load R/4.3.2-gfbf-2023a
# ===================================


# In this script, THE ONLY variables that need to be modified are the followings: 
#### Specify pathways for tools (folder and executable):
runplink="plink" #path to run plink

# ===================================
# IMPUTATION SERVER DATA
#impdata="https://imputationserver.sph.umich.edu/get/{YOUR_LINK}" #URL to download the imputed data from web server
#imppassword="{YOUR_PASSWORD}"
# ===================================

##### The current variables for PLINK are going to be the following:
s4mymaf=0.001 #step4
s5mymissingness1=0.2 #step5
s5mymissingness2=0.02 #step5

# At this moment, in order to run the script, you need to be in the current folder from a terminal and type the following command: "bash 03_POSTIMPUTATION_DATA_PIPELINE_PRC_v4.0.sh"

# And once it finish, your final post-imputed data is going to be under the name "POS....".

# PRINT THE PARAMETRES OF THE CURRENT RUN:
printf "### Script: 03 POSTIMPUTATION PIPELINE FOR GENOTYPED DATA
### Author: POL RAMON CAÑELLAS
### Version: v4.0

### PARAMETERS OF THE CURRENT RUN:
Step 4: maf = %s
Step 5: missigness 1 = %s
Step 5: missigness 2 = %s" "$s4mymaf" "$s5mymissingness1" "$s5mymissingness2" > 03_POSTIMPUTATION_DATA_PIPELINE_PRC_v4.0.log 

echo "### Script running...enjoy!"


###################################################################
### Step 1: IMPUTED DATA MERGING ###
#echo "### Starting Step 1: imputed data merging ###"

# Rearrange imputated data manually using terminal commands.
# FIRST, download manually from Michigan all the files resulting from the run: chromosome files (.zip and .log) as well as other run info.
	curl -sL $impdata | bash

# SECOND, unzip all the chromosome files using the Michigan password from the run:
for chr in {1..22}; do unzip -P $imppassword chr_${chr}.zip; done

# THIRD, merge all chromosme files into one.
bcftools concat -O z -o 1_AllChromosomes.vcf.gz chr1.dose.vcf.gz chr2.dose.vcf.gz chr3.dose.vcf.gz chr4.dose.vcf.gz chr5.dose.vcf.gz chr6.dose.vcf.gz chr7.dose.vcf.gz chr8.dose.vcf.gz chr9.dose.vcf.gz chr10.dose.vcf.gz chr11.dose.vcf.gz chr12.dose.vcf.gz chr13.dose.vcf.gz chr14.dose.vcf.gz chr15.dose.vcf.gz chr16.dose.vcf.gz chr17.dose.vcf.gz chr18.dose.vcf.gz chr19.dose.vcf.gz chr20.dose.vcf.gz chr21.dose.vcf.gz chr22.dose.vcf.gz

#echo "### End of Step 1: imputed data merging ###"


###################################################################
### Step 2: QC TARGET - INFO SCORE ###
#echo "### Starting Step 2: qc target - info score ###"

# Filter according info (r2) results from merged files and keep compressed.
bcftools filter -Oz -i 'INFO/R2>0.9' 1_AllChromosomes.vcf.gz > 2_Target.imputated.vcf.gz
	
# Create data into bfiles.
$runplink --vcf 2_Target.imputated.vcf.gz --make-bed --const-fid --out Target.imputated
for i in $( find . -name 'Target.imputated.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "2_$i" ; done
	                                                                                                   
# Delete big intermedied files (TB!):
rm 2_Target.imputated.vcf.gz

echo "### End of Step 2: qc target - info score ###"


###################################################################
### Step 3: DATA ADAPTATION ###
echo "### Starting Step 3: data adaptation ###"

# For SNPS in .bim file: a) remove palindromic, b) remove indels, c) assign missing IDs, d) remove duplicated position, and e) remove duplicated IDs:

	# a) Obtain and remove palindromic SNPs.
	awk '($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "C" && $6 == "G") || ($5 == "G" && $6 == "C")' 2_TOPMED_Target.imputated.bim | cut -f2 > 3A_palindromic_snps.txt
	$runplink --bfile 2_TOPMED_Target.imputated --exclude 3A_palindromic_snps.txt --make-bed --out TOPMED_Target.imputated
	for i in $( find . -name 'TOPMED_Target.imputated.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "3A_$i" ; done


	# b) Remove those snps with insertions and deletions in the alleles A1 and A2 (I/D, or -), in order to work only with two alternative bases per variant.
	### --snps-only excludes all variants with one or more multi-character allele codes. With 'just-acgt', variants with single-character allele codes outside of {'A', 'C', 'G', 'T', 'a', 'c', 'g','t', <missing code>} are also excluded.
	$runplink --bfile TOPMED_Target.imputated --snps-only 'just-acgt' --make-bed --out TOPMED_Target.imputated
	for i in $( find . -name 'TOPMED_Target.imputated.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "3B_$i" ; done


	# c) Missing IDs.
	### If the .bim file contains SNPs without an identifier (these SNPs may be are indicated as "."), we will assign unique indentifiers to those SNPs.
	### Identify those SNPs that has no IDs.
	$runplink --bfile TOPMED_Target.imputated --set-missing-var-ids @:#:\$1:\$2 --make-bed --out TOPMED_Target.imputated
	for i in $( find . -name 'TOPMED_Target.imputated.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "3C_$i" ; done
                                                                                                          

	# d) Remove SNPS that have same positon.
	### Identify IDs which have same position.
	$runplink --bfile TOPMED_Target.imputated --list-duplicate-vars --out 3D_IDduplicated
	
	### Exclude IDs which position is duplicated.
	$runplink --bfile TOPMED_Target.imputated --exclude 3D_IDduplicated.dupvar --make-bed --out TOPMED_Target.imputated
	for i in $( find . -name 'TOPMED_Target.imputated.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "3D_$i" ; done


	# e) Remove SNPS that have same ID.
	### Identify duplicated IDs.
	$runplink --bfile TOPMED_Target.imputated --write-snplist --out 3D_ID
	awk 'NR==FNR{a[$1]++;next;}{ if (a[$1] > 1)print;}' 3D_ID.snplist 3D_ID.snplist | sort | uniq > 3D_IDduplicated.txt
		
	### Exclude IDs which ID is duplicated.
	$runplink --bfile TOPMED_Target.imputated --exclude 3D_IDduplicated.txt --make-bed --out TOPMED_Target.imputated
	for i in $( find . -name 'TOPMED_Target.imputated.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "3E_$i" ; done

echo "### End of Step 3: data adaptation ###"


###################################################################
### Step 4: QC TARGET - MINOR ALLELE FREQUENCY ###
echo "### Starting Step 4: qc target - minor allelle frequency ###"

# Generate a plot of the MAF distribution.
$runplink --bfile TOPMED_Target.imputated --freq --out 3E_MAF_check
Rscript --no-save POS4_MAF_check.R
mv MAF_distribution.pdf 3E_MAF_distribution.pdf

# Remove SNPs with a low MAF frequency ($s4mymaf).
$runplink --bfile TOPMED_Target.imputated --maf $s4mymaf --make-bed --out TOPMED_Target.imputated_maf
for i in $( find . -name 'TOPMED_Target.imputated_maf.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "4_$i" ; done
	
	# Visualise the remaining MAF distribution.
	$runplink --bfile TOPMED_Target.imputated_maf --freq --out 4_MAF_check
	Rscript --no-save POS4_MAF_check.R
	mv MAF_distribution.pdf 4_MAF_distribution.pdf

echo "### End of Step 4: qc target - minor allelle frequency ###"


###################################################################
### Step 5:  QC TARGET - SNPs MISSINGNESS ###
echo "### Starting Step 5: qc target - SNPs missingness ###"

# Investigate missingness per individual and per SNP and make histograms.
$runplink --bfile TOPMED_Target.imputated_maf --missing
	# Output: plink.imiss and plink.lmiss, these files show respectively the proportion of missing SNPs per individual and the proportion of missing individuals per SNP.
	# Generate plots to visualize the missingness results.
	Rscript --no-save POS5_Hist_miss.R
	for f in plink.* hist*; do mv "$f" "4_$f" ; done

# Delete SNPs with high levels of missingness.
# The following two QC commands will not remove any SNPs or individuals. However, it is good practice to start the QC with these non-stringent thresholds.

# Delete SNPs with missingness >0.2 ($s5mymissingness1).
$runplink --bfile TOPMED_Target.imputated_maf --geno $s5mymissingness1 --make-bed --out TOPMED_Target.imputated_maf
for i in $( find . -name 'TOPMED_Target.imputated_maf.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "5A_$i" ; done

	# Generate plots to visualize the missingness results.
	$runplink --bfile TOPMED_Target.imputated_maf --missing
	Rscript --no-save POS5_Hist_miss.R
	for f in plink.* hist*; do mv "$f" "5A_$f" ; done


# Delete SNPs with missingness >0.02 ($s5mymissingness2).
$runplink --bfile TOPMED_Target.imputated_maf --geno $s5mymissingness2 --make-bed --out TOPMED_Target.imputated_maf
for i in $( find . -name 'TOPMED_Target.imputated_maf.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "5B_$i" ; done

	# Generate plots to visualize the missingness results.
	$runplink --bfile TOPMED_Target.imputated_maf --missing
	Rscript --no-save POS5_Hist_miss.R
    for f in plink.* hist*; do mv "$f" "5B_$f" ; done

echo "### End of Step 5: qc target - SNPs missingness ###"

###################################################################
### ENDING ###

# Rename Bfiles to POST:
mv TOPMED_Target.imputated_maf.bed TOPMED_Target.POSTimputed.bed
mv TOPMED_Target.imputated_maf.bim TOPMED_Target.POSTimputed.bim
mv TOPMED_Target.imputated_maf.fam TOPMED_Target.POSTimputed.fam

# Move all files created into another subdirectory:
rm -r FilesCreated; mkdir FilesCreated
mv 1* 2* 3* 4* 5* ./FilesCreated/

# Remove unnecesary input Bfiles created cretated by plink (ended with ~)
rm *~

echo "### CONGRATULATIONS! You've just succesfully completed the postimputation tutorial! You now must conduct the PRS process."

###################################################################
