###################################################################
######################### START ANALISIS ##########################
###################################################################
############### Script: 01 QC TARGET DATA PIPELINE ################
################### Author: POL RAMON CAÑELLAS ####################
########################## Version: v4.1 ##########################
###################################################################


###################################################################
### IMPORTANT INSTRUCTIONS TO RUN THE PIPELINE OF THIS SCRIPT ###

# Run the present script into the folder containing all files from your data.

# RENAME ORIGINAL FILES: original files needs to be renamed adding at the begining of the name "Target_data" (for .bed/.bim/.fam files) keeping also the original name.
	#### According to the original database info, in case it is needed to take in consideration Phenotypes and the .fam file DO NOT contain only and specifically in this order the information (FID | IID | PHENO) with PHENO values 1(for Healthy Controls), 2(for Patients), and -9 (for NA), this file needs to be transformed using the information of another file with the renamed name "Target_data_pheno_*****" (for .txt file), and the "QC0_Target_data_phenotype.R" script for the STEP 0 needs to be adapted according that file.

# FOLDER MUST CONTAINING FILES:
	####01_QC_TARGET_DATA_PIPELINE_PRC_v4.1.sh
	####Target_data_*****.ped 			# <- plink .ped file 
	####Target_data_*****.map 			# <- plink .map file
	####Target_data_*****.bed 			# <- plink .bed file
	####Target_data_*****.bim 			# <- plink .bim file
	####Target_data_*****.fam 			# <- plink .fam file
	####Target_data_*****.txt *(MUST be present if .fam file is NOT correct, in order to perform step 0 and rearrange it)
	####QC0_Target_data_phenotype.R *(MUST be rearranged according to the needs of the .fam file)
	####QC1_Hist_miss.R
	####QC2_Gender_check.R
	####QC3_MAF_check.R
	####QC4_HWE.R
	####QC5_Check_heterozygosity_rate.R
	####QC5_Inversion.txt
	####QC6_Relatedness.R

# ===================================
# VARIABLES TO DEFINE
DATA="Target_data"

# MODULES NEEDED
module load PLINK/1.9b_6.21-x86_64
module load R/4.3.2-gfbf-2023a
# ===================================


# In this script, THE ONLY variables that need to be modified are the followings: 
	#### Specify plink executable (folder and executable):
	runplink="plink"

	##### Do you need to perform the following steps? (Answer: "YES" or "NO"):
		##### Step 0: RAW FILES SCRIPT ADAPTATION.
			##### .fam file transformation of the data?
			step0fam="NO" #step0

		##### Step 1: MISSIGNESS.
			step1="YES" #step1
		##### Step 2: SEX DISCREPANCY.
			step2="YES" #step2
		##### Step 3: MINOR ALLELE FREQUENCY
			step3="YES" #step3
		##### Step 4: HARDY-WEINBERG EQUILIBRIUM.
			step4="YES" #step4
		##### Step 5: HETEROZYGOSITY.
			step5="YES" #step5
		##### Step 6: RELATEDNESS.
			step6="YES" #step6

	##### The current variables for PLINK are going to be the following:
	mymissigness1=0.2 #step1
	mymissigness2=0.02 #step1
	mymaf=0.01 #step3
	myhwe1=1e-6 #step4
	myhwe2=1e-10 #step4
	mywindowsize=50 #step5
	myshiftwindow=5 #step5
	mypairwiser2=0.2 #step5
	myrelatedness=0.125 #step6

	#### NOTE: 'rename' package is needed. Thus, it is going to be checked and installed if needed.


# At this moment, in order to run the script, you need to be in the current folder from a terminal and type the following command: "bash 01_QC_TARGET_DATA_PIPELINE_PRC_v4.0.sh"

# And once it finish, your final QC data is going to be under the name "QCtargetdata.bed/bim/fam".

# PRINT AND STORE THE PARAMETRES OF THE CURRENT RUN:
printf "### Script: 01 QC TARGET DATA PIPELINE
### Author: POL RAMON CAÑELLAS
### Version: v4.1

### PARAMETERS OF THE CURRENT RUN:
Step 0: .fam adaptation = %s
Step 1: MISSIGNESS = %s
Step 1: missigness 1 = %s
Step 1: missigness 2 = %s
Step 2: SEX DISCREPANCY = %s
Step 3: MINOR ALLELE FREQUENCY = %s
Step 3: maf = %s
Step 4: HARDY-WEINBERG EQUILIBRIUM = %s
Step 4: hwe 1 = %s
Step 4: hwe 2 = %s
Step 5: HETEROZYGOSITY = %s
Step 5: window size = %s
Step 5: shift window = %s
Step 5: pairwise r2 = %s
Step 6: RELATEDNESS = %s
Step 6: relatedness = %s" "$step0fam" "$step1" "$mymissigness1" "$mymissigness2" "$step2" "$step3" "$mymaf" "$step4" "$myhwe1" "$myhwe2" "$step5" "$mywindowsize" "$myshiftwindow" "$mypairwiser2" "$step6" "$myrelatedness" > 01_QC_TARGET_DATA_PIPELINE_PRC_v4.1.log


echo "### Script running...enjoy!"


###################################################################
### Step 0: RAW FILES SCRIPT ADAPTATION ###
echo "### Starting Step 0: raw files script adaptation ###"

#0.1# Convert .ped and .map files to binary files
$runplink --file ${DATA} --make-bed --out Target_data


#0.2# FOR .fam file: 
	##### If it is needed, adapt first the script "0_Target_data_phenotype.R" and exectue it.
	##### This script replace values of Phenotype column according to standard .fam file: Control (=1), Paciente (=2) and NA (=-9).
	if [ "$step0fam" = "YES" ];
	  then
	     echo "# Running Step 0 phenotypes"
	     # Run script to create proper raw phenotype file, as well as summary of FID_IID_Pheno (0A_phenotypes.txt) used next.
	     	Rscript --no-save QC0_Target_data_phenotype.R
	     
	    # Update phenotypes in bfiles.
		$runplink --bfile Target_data --pheno 0A_phenotypes.txt --allow-no-sex --make-bed --out Target_data
		for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0A_$i" ; done
	
	     # Remove individuals with no phenotype.
		awk '$3==-9 {print $1,$2}' 0A_phenotypes.txt > 0B_no.pheno.iid.txt
		$runplink --bfile Target_data --remove 0B_no.pheno.iid.txt --allow-no-sex --make-bed --out Target_data
		for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0B_$i" ; done
	     
	  else
	     echo "# NO performed Step 0 phenotypes"
	fi

	
#0.3# FOR .bim file: remove indels to work with snps.

	##### Remove those snps with insertions and deletions in the alleles A1 and A2 (I/D, or -), in order to work only with two alternative bases per variant.
	##### --snps-only excludes all variants with one or more multi-character allele codes. With 'just-acgt', variants with single-character allele codes outside of {'A', 'C', 'G', 'T', 'a', 'c', 'g','t', <missing code>} are also excluded.
	$runplink --bfile Target_data --snps-only 'just-acgt' --make-bed --out Target_data
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0C_$i" ; done


echo "### End of Step 0: raw files script adaptation ###"


###################################################################
### Step 1: MISSIGNESS ###

if [ "$step1" = "NO" ];
   then
	echo "### NO performed Step 1: Missigness ###"   
   else
	echo "### Starting Step 1: Missigness ###"

	# Investigate missingness per individual and per SNP and make histograms.
	$runplink --bfile Target_data --missing
		# Output: plink.imiss and plink.lmiss, these files show respectively the proportion of missing SNPs per individual and the proportion of missing individuals per SNP.
		# Generate plots to visualize the missingness results.
		Rscript --no-save QC1_Hist_miss.R
		for f in plink.* hist*; do mv "$f" "0C_$f"; done

	# Delete SNPs and individuals with high levels of missingness, explanation of this and all following steps can be found in box 1 and table 1 of the article mentioned in the comments of this script.
	# The following two QC commands will not remove any SNPs or individuals. However, it is good practice to start the QC with these non-stringent thresholds.

	# Delete SNPs with missingness >0.2 ($mymissigness1).
	$runplink --bfile Target_data --geno $mymissigness1 --make-bed --out Target_data
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "1A_$i" ; done
		# Generate plots to visualize the missingness results.
		$runplink --bfile Target_data --missing
		Rscript --no-save QC1_Hist_miss.R
		for f in plink.* hist*; do mv "$f" "1A_$f" ; done

	# Delete individuals with missingness >0.2 ($mymissigness1).
	$runplink --bfile Target_data --mind $mymissigness1 --make-bed --out Target_data
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "1B_$i" ; done
		# Generate plots to visualize the missingness results.
		$runplink --bfile Target_data --missing
		Rscript --no-save QC1_Hist_miss.R
		for f in plink.* hist*; do mv "$f" "1B_$f" ; done
		

	# Delete SNPs with missingness >0.02 ($mymissigness2).
	$runplink --bfile Target_data --geno $mymissigness2 --make-bed --out Target_data
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "1C_$i" ; done
		# Generate plots to visualize the missingness results.
		$runplink --bfile Target_data --missing
		Rscript --no-save QC1_Hist_miss.R
		for f in plink.* hist*; do mv "$f" "1C_$f" ; done

	# Delete individuals with missingness >0.02 ($mymissigness2).
	$runplink --bfile Target_data --mind $mymissigness2 --make-bed --out Target_data
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "1D_$i" ; done
		# Generate plots to visualize the missingness results.
		$runplink --bfile Target_data --missing
		Rscript --no-save QC1_Hist_miss.R
		for f in plink.* hist*; do mv "$f" "1D_$f" ; done

	echo "### End of Step 1: Missigness ###"
fi


###################################################################
### Step 2: SEX DISCREPANCY ####

if [ "$step2" = "NO" ];
   then
	echo "### NO performed Step 2: Sex discrepancy ###"   
   else
	echo "### Starting Step 2: Sex discrepancy ###"
	#### Afegir warning per si no esta disponible el txt amb el phenotype 
	#### Amb aquesta linia actualitzem amb un fitxer .txt el phenotip de sexes [FID, IID, SEX] (sense capçalera) 
    plink --bfile Target_data --update-sex Target_data_pheno_sex.txt --make-bed --out Target_data

	# Check for sex discrepancy.
	# Subjects who were a priori determined as females must have a F value of <0.2, and subjects who were a priori determined as males must have a F value >0.8. This F value is based on the X chromosome inbreeding (homozygosity) estimate.
	# Subjects who do not fulfil these requirements are flagged "PROBLEM" by PLINK.

	$runplink --bfile Target_data --check-sex

		# Generate plots to visualize the sex-check results.
		Rscript --no-save QC2_Gender_check.R
		for f in plink.* Gender_* Men_* Women*; do mv "$f" "1D_$f" ; done
		# These checks indicate that there is one woman with a sex discrepancy, F value of 0.99. (When using other datasets often a few discrepancies will be found).


	# Delete individuals with sex discrepancy.
	grep "PROBLEM" 1D_plink.sexcheck | awk '{print$1,$2}'> 2_sex_discrepancy.txt
	# This command generates a list of individuals with the status PROBLEM.
	$runplink --bfile Target_data  --remove 2_sex_discrepancy.txt --make-bed --out Target_data 
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "2_$i" ; done
	# This command removes the list of individuals with the status PROBLEM.
		
		# Check again for sex discrepancy.
		$runplink --bfile Target_data --check-sex
		
		# Generate plots to visualize the sex-check results.
		Rscript --no-save QC2_Gender_check.R
		for f in plink.* Gender_* Men_* Women*; do mv "$f" "2_$f" ; done


	echo "### End of Step 2: Sex discrepancy ###"
fi


###################################################################
### Step 3: MINOR ALLELE FREQUENCY ###

if [ "$step3" = "NO" ];
   then
	echo "### NO performed Step 3: minor allelle frequency ###"   
   else
	echo "### Starting Step 3: minor allelle frequency ###"

	# Generate a bfile with autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF).

	# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).
	awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' Target_data.bim > 2_snp_1_22.txt
	$runplink --bfile Target_data --extract 2_snp_1_22.txt --make-bed --out Target_data
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "3A_$i" ; done

	# Generate a plot of the MAF distribution.
	$runplink --bfile Target_data --freq --out 3A_MAF_check
	Rscript --no-save QC3_MAF_check.R
	mv "MAF_distribution.pdf" "3A_MAF_distribution.pdf"

	# Remove SNPs with a low MAF frequency ($mymaf).
	# A conventional MAF threshold for a regular GWAS is between 0.01 or 0.05, depending on sample size.
	$runplink --bfile Target_data --maf $mymaf --make-bed --out Target_data
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "3B_$i" ; done
		
		# Visualise the remaining MAF distribution.
		$runplink --bfile Target_data --freq --out 3B_MAF_check
		Rscript --no-save QC3_MAF_check.R
		mv "MAF_distribution.pdf" "3B_MAF_distribution.pdf"

	echo "### End of Step 3: minor allelle frequency ###"
fi


###################################################################
### Step 4: HARDY-WEINBERG EQUILIBRIUM ###

if [ "$step4" = "NO" ];
   then
	echo "### NO performed Step 4: hardy-weinberg equilibrium ###"   
   else
	echo "### Starting Step 4: hardy-weinberg equilibrium ###"

	# Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).

	# Check the distribution of HWE p-values of all SNPs.
	$runplink --bfile Target_data --hardy
	for f in plink.*; do mv "$f" "3B_$f" ; done

	# Selecting SNPs with HWE p-value below 0.00001, required for one of the two plot generated by the next Rscript, allows to zoom in on strongly deviating SNPs.
	Rscript --no-save QC4_HWE.R

	# By default the --hwe option in plink only filters for controls.
	# Therefore, we use two steps, first we use a stringent HWE threshold for controls, followed by a less stringent threshold for the case data.
	$runplink --bfile Target_data --hwe $myhwe1 --make-bed --out Target_data
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "4A_$i" ; done

	# The HWE threshold for the cases filters out only SNPs which deviate extremely from HWE.
	# This second HWE step only focusses on cases because in the controls all SNPs with a HWE p-value < hwe 1e-6 were already removed.
	$runplink --bfile Target_data --hwe $myhwe2 --hwe-all --make-bed --out Target_data
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "4B_$i" ; done

	echo "### End of Step 4: hardy-weinberg equilibrium ###"
fi


###################################################################
### Step 5: HETEROZYGOSITY ###

if [ "$step5" = "NO" ];
   then
	echo "### NO performed Step 5: heterozygosity ###"   
   else
	echo "### Starting Step 5: heterozygosity ###"

	# Generate a plot of the distribution of the heterozygosity rate of your subjects.
	# And remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.

	# Checks for heterozygosity are performed on a set of SNPs which are not highly correlated.
	# Therefore, to generate a list of non-(highly)correlated SNPs, we exclude high inversion regions (QC5_Inversion.txt [High LD regions]) and prune the SNPs using the command --indep-pairwise.
	# The parameters 50 5 0.2 stand respectively for: the window size ($mywindowsize), the number of SNPs to shift the window at each step ($myshiftwindow), and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously ($mypairwiser2).
	$runplink --bfile Target_data --exclude QC5_Inversion.txt --range --indep-pairwise $mywindowsize $myshiftwindow $mypairwiser2 --out 4B_indepSNP

	$runplink --bfile Target_data --extract 4B_indepSNP.prune.in --het --out 5_pruned.data
	# The output file contains your pruned data set.

	#1st Plot of the heterozygosity rate distribution with the +-3SD thresholds.
	#2nd Generates a list of individuals (outliers) who deviate more than 3 standard deviations (+-3SD) from the heterozygosity rate mean.
	Rscript --no-save QC5_Check_heterozygosity_rate.R

	# Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns.
	sed 's/"// g' 5_fail-het-qc.txt | awk '{print$1, $2}'> 5_het_fail_ind.txt

	# Remove heterozygosity rate outliers.
	$runplink --bfile Target_data --remove 5_het_fail_ind.txt --make-bed --out Target_data
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "5_$i" ; done

	echo "### End of Step 5: heterozygosity ###"
fi


###################################################################
### Step 6: RELATEDNESS ###

if [ "$step6" = "NO" ];
   then
	echo "### NO performed Step 6: relatedness ###"   
   else
	echo "### Starting Step 6: relatedness ###"

	# It is essential to check datasets you analyse for cryptic relatedness.
	# Assuming a random population sample we are going to exclude all individuals above the pihat threshold of 0.125 in this tutorial.

	# Check for relationships between individuals with a pihat > 0.125 ($myrelatedness).
	$runplink --bfile Target_data --extract 4B_indepSNP.prune.in --genome --min $myrelatedness --out 5_pihat_min0.125

	# If your dataset contain parent-offspring relations, the following commands will visualize specifically these parent-offspring relations, using the z values.
	awk '{ if ($8 >0.9) print $0 }' 5_pihat_min0.125.genome > 5_zoom_pihat.genome

	# Generate a plot to assess the type of relationship.
	Rscript --no-save QC6_Relatedness.R

	# Explanation of the generated plot according to the individuals relationship fo the dataset: PO = parent-offspring, UN = unrelated individuals, OT = others.
	# Normally, family based data should be analyzed using specific family based methods. In this tutorial, for demonstrative purposes, we treat the relatedness as cryptic relatedness in a random population sample.
	# In this tutorial, we aim to remove all 'relatedness' from our dataset.
	# To demonstrate that the majority of the relatedness was due to parent-offspring we only include founders (individuals without parents in the dataset).
	$runplink --bfile Target_data --filter-founders --make-bed --out Target_data
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "6A_$i" ; done

	# Now we will look again for individuals with a pihat >0.125.
	$runplink --bfile Target_data --extract 4B_indepSNP.prune.in --genome --min $myrelatedness --out Target_data
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "6B_$i" ; done
	# The file '6B_Target_data.genome' shows that, after exclusion of all non-founders, LOTS OF individuals pair with a pihat greater than 0.125 remains in the Cibersam data.

	# Get the Delete the founders individuals with a pihat >0.125.
	awk '{print$1,$2}' 6B_Target_data.genome | sort | uniq > 6B_Target_data.txt

	# Delete the founders individuals with a pihat >0.125.
	$runplink --bfile Target_data --remove 6B_Target_data.txt --make-bed --out Target_data
	for i in $( find . -name 'Target_data.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "6C_$i" ; done

	echo "### End of Step 6: relatedness ###"
fi


###################################################################
### ENDING ###

# Rename Bfiles to QC:
mv Target_data.bed QCtargetdata.bed
mv Target_data.bim QCtargetdata.bim
mv Target_data.fam QCtargetdata.fam

# Move all files created into another subdirectory:
rm -r FilesCreated; mkdir FilesCreated
mv 0[ABCDEFGHI]* 1* 2* 3* 4* 5* 6* ./FilesCreated/

# Remove unnecesary input Bfiles created cretated by plink (ended with ~)
rm *~

echo "### CONGRATULATIONS! You've just succesfully completed the QC tutorial! You are now able to conduct a proper genetic QC."


###################################################################
