###################################################################
######################### START ANALISIS ##########################
###################################################################
###### Script: 02 PREIMPUTATION PIPELINE FOR GENOTYPED DATA #######
################### Author: POL RAMON CAÑELLAS ####################
########################## Version: v5.0 ##########################
###################################################################


###################################################################
### IMPORTANT INSTRUCTIONS TO RUN THE PIPELINE OF THIS SCRIPT ###

# =========================================
# MODULES NEEDED
module load PLINK/1.9b_6.21-x86_64
module load R/4.3.2-gfbf-2023a
module load BioPerl/1.7.8-GCCcore-12.3.0
module load BCFtools/1.21-GCC-12.3.0
# =========================================

# Run the present script into the folder containing all files from your data.

# RENAME QC FILES (.bed/.bim/.fam files) AND TARGET PEDIGREE FILES (.ped): this files needs to be renamed adding at the begining of the name the following:
	#### FOR TARGET DATA: "QCtargetdata_*****".
	#### FOR REFERENCE DATA: "QCreferencedata_*****".
	#### Additionally, check the reference genome build of both datasets (GRCh37/hg19 or GRCh38/hg38). If they are NOT under the same, then a LiftOver is needed and thus, the variable "qcliftover" for Step0 is needed to set as "YES".

# FOLDER MUST CONTAINING FILES:
	####02_PREIMPUTATION_DATA_PIPELINE_PRC_v5.0.sh
	####QCtargetdata_*****.bed
	####QCtargetdata_*****.bim
	####QCtargetdata_*****.fam
	####QCreferencedata_*****.bed
	####QCreferencedata_*****.bim
	####QCreferencedata_*****.fam
	####PRE0_FID_IID.R
	####PRE0D_Duplicates.R
	####PRE3_PCA.R


# In this script, THE ONLY variables that need to be modified are the followings: 
#### Specify plink and liftover executable (folder and executable):
runplink="plink"
#runliftover="/home/pol/bin/liftOver/liftOver" #path to UCSC liftOver tool
#liftoverChain="/home/pol/bin/liftOver/hg38ToHg19.over.chain" #chain file to conversion hg38 to hg19

##### Do you need to perform this steps? Answer: "YES" or "NO"
s0liftoverTarget="NO" #step0
s0liftoverReference="NO" #step0

##### The current variables for PLINK are going to be the following:
mywindowsize=50 #step3
myshiftwindow=5 #step3
mypairwiser2=0.2 #step3

# At this moment, in order to run the script, you need to be in the current folder from a terminal and type the following command: "bash 02_PREIMPUTATION_DATA_PIPELINE_PRC_v5.0_MGR.sh"

# And once it is finished, your final pre-imputed data is going to be under the name "2_PreImputed.split.XX.vcf.gz".

# PRINT AND STORE THE PARAMETERS OF THE CURRENT RUN:
printf "### Script: 02 PREIMPUTATION PIPELINE FOR GENOTYPED DATA
### Author: POL RAMON CAÑELLAS
### Version: v5.0

### PARAMETERS OF THE CURRENT RUN:
Step 0: LiftOver Target = %s
Step 0: LiftOver Reference = %s
Step 3: window size = %s
Step 3: shift window = %s
Step 3: pairwise r2 = %s" "$s0liftoverTarget" "$s0liftoverReference" "$mywindowsize" "$myshiftwindow" "$mypairwiser2" > 02_PREIMPUTATION_DATA_PIPELINE_PRC_v5.0.log

echo "### Script running...enjoy!"


###################################################################
### Step 0: QC BFILES ADAPTATION  ###
echo "### Starting Step 0: qc bfiles adaptation ###"

#0.1# Rename .bim/.bed/.fam files as "QC*.*" for the pipeline:
	#"QCtargetdata_*.*" as "QCtargetdata.*":
	cp QCtargetdata_*.bed QCtargetdata.bed
	cp QCtargetdata_*.bim QCtargetdata.bim
	cp QCtargetdata_*.fam QCtargetdata.fam
	
	#"QCreferencedata_*.*" as "QCreferencedata.*":
	cp QCreferencedata_*.bed QCreferencedata.bed
	cp QCreferencedata_*.bim QCreferencedata.bim
	cp QCreferencedata_*.fam QCreferencedata.fam


#0.2# Arrange SNPs in target and reference data:

	#A) Palindromic SNPs (ambiguous).
	#### Obtain and remove palindromic SNPs.
	
		##### ONLY for Targetdata
		awk '($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "C" && $6 == "G") || ($5 == "G" && $6 == "C")' QCtargetdata.bim | cut -f2 > 0A_palindromic_snps.txt

		$runplink --bfile QCtargetdata --exclude 0A_palindromic_snps.txt --make-bed --out QCtargetdata
		for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0A_$i" ; done

	#B) Missing rsIDs.
	#### If the .bim file contains SNPs without an identifier (these SNPs may be are indicated as "."), we will assign unique indentifiers to those SNPs.
	#### Identify those SNPs that have no IDs and rename them with unique identifiers.
	
		##### For Targetdata
		$runplink --bfile QCtargetdata --set-missing-var-ids @:#:\$1:\$2 --make-bed --out QCtargetdata
		for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0B_$i" ; done
				
		##### For Referencedata
		$runplink --bfile QCreferencedata --set-missing-var-ids @:#:\$1:\$2 --make-bed --out QCreferencedata
		for i in $( find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0B_$i" ; done
			

	#C) Duplicated rsIDs.
	#### Identify SNPS that have same rsID and exclude them as duplicated IDs.
				
		##### For Targetdata
		$runplink --bfile QCtargetdata --write-snplist --out 0B_targetID.snp
		awk 'NR==FNR{a[$1]++;next;}{ if (a[$1] > 1)print;}' 0B_targetID.snp.snplist 0B_targetID.snp.snplist | sort | uniq > 0B_targetID.duplicated.txt
					
		$runplink --bfile QCtargetdata --exclude 0B_targetID.duplicated.txt --make-bed --out QCtargetdata
		for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0C_$i" ; done
					
		##### For Referencedata
		$runplink --bfile QCreferencedata --write-snplist --out 0B_referenceID.snp
		awk 'NR==FNR{a[$1]++;next;}{ if (a[$1] > 1)print;}' 0B_referenceID.snp.snplist 0B_referenceID.snp.snplist | sort | uniq > 0B_referenceID.duplicated.txt
					
		$runplink --bfile QCreferencedata --exclude 0B_referenceID.duplicated.txt --make-bed --out QCreferencedata
		for i in $( find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0C_$i" ; done		
	
			
	#D) Same position (CH, POS, A1, A2).		
	#### Remove SNPS that have same positon (this means chromosome, position, A1, and A2). NOTE: THOSE THAT HAVE DIFFERENT A1 OR A2 WILL NOT BE REMOVED!! THIS IS GOING TO BE DEAL IN NEXT STEP 'E'.
	#### Check frequencies. Identify IDs which have same position. Compare the frequencies of the same positions according the IDs. Exclude IDs which chromosome, position, A1, and A2, is duplicated and have lower NCHROBS (this means lower frequency) in order to keep only the ones with higer frecueny in this position, independently of the name.
	
		##### Get the list of duplicates and the frequencies.
		$runplink --bfile QCtargetdata --list-duplicate-vars --freq --out 0C_targetID.duplicated
		$runplink --bfile QCreferencedata --list-duplicate-vars --freq --out 0C_referenceID.duplicated

		##### Get the list of ids with the lowest frecuencies of each duplicated.
		Rscript --no-save PRE0D_Duplicates.R
	
		##### Exclude the ids ####<=== HE HAGUT D'AFEGIR 3 FILES PER SEPARAT PERQUE SOLAMENT ESTAVEN AMB LA ~.
		$runplink --bed QCtargetdata.bed~ --bim QCtargetdata.bim~ --fam QCtargetdata.fam~ --exclude 0C_targetID.duplicated.exclude.txt --make-bed --out QCtargetdata
		for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0D_$i" ; done
		
		$runplink --bfile QCreferencedata --exclude 0C_referenceID.duplicated.exclude.txt --make-bed --out QCreferencedata
		for i in $( find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0D_$i" ; done
		
			
	#E) Same position (CH, POS).
	#### Set new ID to SNPs as Chr:BP. 
	
		##### For Targetdata
		awk  '{print $2, $1":"$4}' QCtargetdata.bim > 0D_targetIDcoordinates.txt				
		$runplink --bfile QCtargetdata --update-name 0D_targetIDcoordinates.txt --make-bed --out QCtargetdata
		for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0E_$i" ; done
				
		##### For Referencedata
		awk  '{print $2, $1":"$4}' QCreferencedata.bim > 0D_referenceIDcoordinates.txt 
		$runplink --bfile QCreferencedata --update-name 0D_referenceIDcoordinates.txt --make-bed --out QCreferencedata
		for i in $( find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0E_$i" ; done		
						
			
	#F) Remove same position (CH, POS).
	#### Remove SNPS that have same rsID (this time with the new IDs meaning that they have different A1 or A2). This step is necessary to avoid dupplicated issues in setting reference Allele in step 1.
	### Identify duplicated rsIDs. Exclude rsIDs which rsID is duplicated.
		
		##### For Targetdata
		$runplink --bfile QCtargetdata --write-snplist --out 0E_targetID.snp
		awk 'NR==FNR{a[$1]++;next;}{ if (a[$1] > 1)print;}' 0E_targetID.snp.snplist 0E_targetID.snp.snplist | sort | uniq > 0E_targetID.duplicated.txt
					
		$runplink --bfile QCtargetdata --exclude 0E_targetID.duplicated.txt --make-bed --out 0F_QCtargetdata
		for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0F_$i" ; done
					
		##### For Referencedata
		$runplink --bfile QCreferencedata --write-snplist --out 0E_referenceID.snp
		awk 'NR==FNR{a[$1]++;next;}{ if (a[$1] > 1)print;}' 0E_referenceID.snp.snplist 0E_referenceID.snp.snplist | sort | uniq > 0E_referenceID.duplicated.txt
					
		$runplink --bfile QCreferencedata --exclude 0E_referenceID.duplicated.txt --make-bed --out QCreferencedata
		for i in $( find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0F_$i" ; done	


#0.3# If both data sets are not under the same Reference Genome, then a flipover is needed for step 
###### (THIS STEP IS NOT TESTED HERE!!)
	if [ "$s0liftoverTarget" = "YES" ];
	   then
	      echo "### Starting Step 0: LiftOver TARGET ###"
	     	# Lift PLINK format (PLINK format usually refers to .ped and .map files).
		# Obtain PLINK.map file:
		$runplink --bfile QCtargetdata --recode --tab --make-bed --out QCtargetdata
		for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0A_$i" ; done
		
		# Convert .map to liftover .bed file
		#### By rearrange columns of .map file, we obtain a standard BED format file for liftover.
		awk '{print "chr"$1,$4,$4+1,$2}' 0A_QCtargetdata.map > 0B_QCtargetdata.bed
		
		# LiftOver .bed file
		chmod +x $runliftover
		chmod +x $liftoverChain
		$runliftover 0B_QCtargetdata.bed $liftoverChain 0C_QCtargetdata.bed 0C_unlifted.bed

		# Convert lifted .bed file back to .map file
		### By rearrange column of .bed file to obtain .map file in the new build.
		awk 'BEGIN {OFS="\t"} ; {print substr($1,4),$4,0,$2}' 0C_QCtargetdata.bed > 0D_QCtargetdata.map	
			
		# Modify .ped file
		### .ped file have many column files. By convention, the first six columns are family_id, person_id, father_id, mother_id, sex, and phenotype. From the 7th column, there are two letters/digits representing a genotype at the certain marker. In step "liftOver", as some genome positions cannot be lifted to the new version, we need to drop their corresponding columns from .ped file to keep consistency. You can use PLINK --exclude those snps, see Remove a subset of SNPs.
		
		### Select snps which contain only CHROMOSOME NUMBERS from the .map file
		awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5 || $1==6 || $1==7 || $1==8 || $1==9 || $1==10 || $1==11 || $1==12 || $1==13 || $1==14 || $1==15 || $1==16 || $1==17 || $1==18 || $1==19 || $1==20 || $1==21 || $1==22' 0D_QCtargetdata.map > 0E_QCtargetdata.map
		awk '{print $2,$4}' 0E_QCtargetdata.map > 0E_ID-POS.txt
		awk '{print $2}' 0E_QCtargetdata.map > 0E_ID.txt
		
		# KEEP ONLY LIFTED SNPS FROM ORIGINAL FILES.
		$runplink --bfile QCtargetdata --extract 0E_ID.txt --make-bed --out QCtargetdata
		for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0F_$i" ; done

		# UPDATE MAP.
		$runplink --bfile QCtargetdata --update-map 0E_ID-POS.txt --make-bed --out QCtargetdata
		for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0G_$i" ; done

		# Set new ID to SNPs as Chr:BP
		awk  '{print $2, $1":"$4}' QCtargetdata.bim > 0G_IDcoordinates.txt
		$runplink --bfile QCtargetdata --update-name 0G_IDcoordinates.txt --make-bed --out QCtargetdata
		for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0H_$i" ; done

		# Rearrange individual IDs
		Rscript --no-save PRE0_FID_IID.R
		$runplink --bfile QCtargetdata --update-ids 0H_FID.IID.txt --make-bed --out QCtargetdata
		for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0I_$i" ; done
		
	      echo "### End of Step 0: LiftOver TARGET ###"
		
		
	elif [ "$s0liftoverReference" = "YES" ];
	   then
	      echo "### Starting Step 0: LiftOver REFERENCE ###"
	     	# Lift PLINK format (PLINK format usually refers to .ped and .map files).
		# Obtain PLINK.map file:
		$runplink --bfile QCreferencedata --recode --tab --make-bed --out QCreferencedata
		for i in $( find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0A_$i" ; done
		
		# Convert plink .map to .bed liftover file
		#### By rearrange columns of .map file, we obtain a standard BED format file for liftover.
		awk '{print "chr"$1,$4,$4+1,$2}' 0A_QCreferencedata.map > 0B_QCreferencedata.bed
		
		# LiftOver .bed file
		chmod +x $runliftover
		chmod +x $liftoverChain
		$runliftover 0B_QCreferencedata.bed $liftoverChain 0C_QCreferencedata.bed 0C_unlifted.bed

		# Convert lifted .bed file back to .map file
		### By rearrange column of .bed file to obtain .map file in the new build.
		awk 'BEGIN {OFS="\t"} ; {print substr($1,4),$4,0,$2}' 0C_QCreferencedata.bed > 0D_QCreferencedata.map	
			
		# Modify .bed file
		### .ped file have many column files. By convention, the first six columns are family_id, person_id, father_id, mother_id, sex, and phenotype. From the 7th column, there are two letters/digits representing a genotype at the certain marker. In step "liftOver", as some genome positions cannot be lifted to the new version, we need to drop their corresponding columns from .ped file to keep consistency. You can use PLINK --exclude those snps, see Remove a subset of SNPs.
		
		### Select snps which contain only CHROMOSOME NUMBERS from the .map file
		awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5 || $1==6 || $1==7 || $1==8 || $1==9 || $1==10 || $1==11 || $1==12 || $1==13 || $1==14 || $1==15 || $1==16 || $1==17 || $1==18 || $1==19 || $1==20 || $1==21 || $1==22' 0D_QCreferencedata.map > 0E_QCreferencedata.map
		awk '{print $2,$4}' 0E_QCreferencedata.map > 0E_ID-POS.txt
		awk '{print $2}' 0E_QCreferencedata.map > 0E_ID.txt
		
		# KEEP ONLY LIFTED SNPS FROM ORIGINAL FILES.
		$runplink --bfile QCreferencedata --extract 0E_ID.txt --make-bed --out QCreferencedata
		for i in $( find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0F_$i" ; done

		# UPDATE MAP.
		$runplink --bfile QCreferencedata --update-map 0E_ID-POS.txt --make-bed --out QCreferencedata
		for i in $( find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0G_$i" ; done

		# Set new ID to SNPs as Chr:BP
		awk  '{print $2, $1":"$4}' QCreferencedata.bim > 0G_IDcoordinates.txt
		$runplink --bfile QCreferencedata --update-name 0G_IDcoordinates.txt --make-bed --out QCreferencedata
		for i in $( find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0H_$i" ; done

		# Rearrange individual IDs
		Rscript --no-save PRE0_FID_IID.R
		$runplink --bfile QCreferencedata --update-ids 0H_FID.IID.txt --make-bed --out QCreferencedata
		for i in $( find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "0I_$i" ; done
		
	      echo "### End of Step 0: LiftOver REFERENCE ###"
			
	fi

echo "### End of Step 0: qc bfiles adaptation ###"


###################################################################
### Step 1: PREPARE TARGET AND REFERENCE DATASET FOR MERGE ###
echo "### Starting Step 1: prepare target and reference dataset for merge ###"

# Check target and reference dataset merge:
### Prior to merging target data with reference data make sure that the files are mergeable conducting the following 3 steps:


# The following steps are maybe quite technical in terms of commands, but we just compare the two data sets and make sure they correspond.

	# 1.1) Make sure the reference genome is similar in both datasets (GRCh37/hg19 or GRCh38/hg38). Otherwise, step0 Liftover is needed to be perform inicially.
	#Then, set reference genome.
	awk '{print$2,$5}' QCreferencedata.bim > 1A_QCreferencedata-list.txt
	$runplink --bfile 0F_QCtargetdata --reference-allele 1A_QCreferencedata-list.txt --make-bed --out QCtargetdata
	for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "1A_$i" ; done


	# 1.2) Resolve strand issues.
	# Check for potential strand issues.
	awk '{print$2,$5,$6}' QCreferencedata.bim > 1A_QCreferencedata_strand.txt
	awk '{print$2,$5,$6}' QCtargetdata.bim > 1A_QCtargetdata_strand.txt
	sort 1A_QCreferencedata_strand.txt 1A_QCtargetdata_strand.txt | uniq -u > 1A_strand_differences.txt

	# Flip SNPs for resolving strand issues.
	#### Print SNP-identifier and remove duplicates.
	awk '{print$1}' 1A_strand_differences.txt | sort -u > 1A_flip_list.txt
	#### Generates a file of SNPs which are the non-corresponded SNPs between the two files. 
	#### Flip the non-corresponding SNPs. 
	$runplink --bfile QCtargetdata --flip 1A_flip_list.txt --reference-allele 1A_QCreferencedata-list.txt --make-bed --out QCtargetdata 
	for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "1B_$i" ; done

	# Check for SNPs which are still problematic after they have been flipped.
	awk '{print$2,$5,$6}' QCtargetdata.bim > 1B_QCtargetdata_problematics.txt
	sort 1A_QCreferencedata_strand.txt 1B_QCtargetdata_problematics.txt |uniq -u  > 1B_QCtargetdata_uncorresponded.txt
	#### This file demonstrates that there are differences between the files.


	# 1.3) Remove the SNPs which after the previous two steps still differ between datasets.
	# Remove problematic SNPs from both data sets.
	awk '{print$1}' 1B_QCtargetdata_uncorresponded.txt | sort -u > 1B_QCtargetdata_excluded.txt
	# The command above generates a list of SNPs which caused the differences between the Target and the Reference data sets after flipping and setting of the reference genome.

	# Remove the problematic SNPs from both datasets.
	$runplink --bfile QCtargetdata --exclude 1B_QCtargetdata_excluded.txt --make-bed --out QCtargetdata
	for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "1C_$i" ; done

	$runplink --bfile QCreferencedata --exclude 1B_QCtargetdata_excluded.txt --make-bed --out QCreferencedata
	for i in $( find . -name 'QCreferencedata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "1C_$i" ; done


echo "### End of Step 1: prepare target and reference dataset for merge ###"


###################################################################
### Step 2: MERGE BOTH DATASETS ###
echo "### Starting Step 2: merge both databases ###"

# First, get the individuals IDs per each dataset that is needed for running the R PCA script in step 5:
	# For target:
	awk '{print$1,$2}' QCtargetdata.fam > 1C_QCtargetdata_individualsID.txt
	
	# For reference:
	awk '{print$1,$2}' QCreferencedata.fam > 1C_QCreferencedata_individualsID.txt

# Second, merge target with reference dataset.
$runplink --bfile QCtargetdata --bmerge QCreferencedata.bed QCreferencedata.bim QCreferencedata.fam --allow-no-sex --make-bed --out PCA_datasets
for i in $( find . -name 'PCA_datasets.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "2_$i" ; done

echo "### End of Step 2: merge both dataset ###"


###################################################################
### Step 3: PRINCIPAL COMPONENT ANALYSIS - MERGED ###
echo "### Starting Step 3: principal component analysis - MERGED ###"

#Pruning.
	#.prune.in file
	$runplink --bfile PCA_datasets --indep-pairwise $mywindowsize $myshiftwindow $mypairwiser2 --out 3A_PCA_datasets
	# create a pruned file
	$runplink --bfile PCA_datasets --extract 3A_PCA_datasets.prune.in --make-bed --out 3B_PCA_datasets

#Run PCA.
$runplink --bfile 3B_PCA_datasets --pca --out 3_PCA_results

#Plot PCA.
Rscript --no-save PRE3_PCA.R
#### Here the ethnical outliers file (individuals excluded) is going to be used in later steps.

echo "### End of Step 3: principal component analysis - MERGED ###"


###################################################################
### Step 4: PRINCIPAL COMPONENT ANALYSIS - QC TARGET  ###
echo "### Starting Step 4: principal component analysis - TARGET ###"

#Pruning.
	#.prune.in file
	$runplink --bfile 0F_QCtargetdata --indep-pairwise $mywindowsize $myshiftwindow $mypairwiser2 --out 4_QCtargetdata
	# create a pruned file
	$runplink --bfile 0F_QCtargetdata --extract 4_QCtargetdata.prune.in --make-bed --out 4_QCtargetdata_PCA
	
# Calculate the first 10 PCs of the target data to obtain the covariate file that is going to be used in the PRS (PRSice). 
$runplink --bfile 4_QCtargetdata_PCA --pca 10 --out 4_PCA_results

# Change name of the file .eigenvec to .covariate
cp 4_PCA_results.eigenvec 4_PCA_results.covariate

echo "### End of Step 4: principal component analysis - TARGET ###"


###################################################################
### Step 5: VCF COOKBOOK FOR IMPUTATION ("https://genepi.github.io/michigan-imputationserver/prepare-your-data/") ###
echo "### Starting Step 5: allele update according Topmed ###"

	# Get the frequencies.
	$runplink --bfile 0F_QCtargetdata --freq --out QCtargetdata
	for i in $( find . -name 'QCtargetdata.*' -printf '%f\n' | awk '!/.bed/ && !/.bim/ && !/.fam/' ); do mv "$i" "5A_$i" ; done


	# Download files needed.
	wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip
	unzip HRC-1000G-check-bim-v4.3.0.zip
	rm -r HRC-1000G-check-bim-v4.3.0.zip
	
	wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz  # Latest versión from: https://www.chg.ox.ac.uk/~wrayner/tools/
	gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

	
	# Reference bfile to HRC-1000G
	perl HRC-1000G-check-bim.pl -b QCtargetdata.bim -f 5A_QCtargetdata.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
	# After running perl script, Run-plink.sh file will be created, with the corresponding QC files (SNPs to rename, SNPs to flip, and SNPs to change the reference allele)
	
	sh Run-plink.sh

echo "### End of Step 5: allele update according Topmed ###"


###################################################################
### Step 6: QC TARGET SPLIT INTO CHROMOSOMES ###
echo "### Starting Step 6: qc target divide into chromosomes ###"

	# Create vcf:
	for chr in {1..22}; do bcftools sort QCtargetdata-updated-chr${chr}.vcf -Oz -o QCtargetdata-updated-chr${chr}.vcf.gz; done
	for chr in {1..22}; do bcftools index --tbi QCtargetdata-updated-chr${chr}.vcf.gz; done

echo "### End of Step 6: qc target divide into chromosomes ###"


###################################################################
### ENDING ###

# Copy PRE files and give an appropiated name:
cp 4_PCA_results.covariate PREtarget_PCA.results.covariate

# Remove unnecesary input Bfiles created cretated by plink (ended with ~)
rm *~

# Move all files created into another subdirectory:
mv QCtargetdata.bed PREIMPUTEDtargetdata.bed; mv QCtargetdata.bim PREIMPUTEDtargetdata.bim; mv QCtargetdata.fam PREIMPUTEDtargetdata.fam
rm -r FilesCreated; mkdir FilesCreated
mv 0[ABCDEFG]* 1* 2* 3* 4* 5* 6* HRC-1000G-check-bim.pl ./FilesCreated/


# To save local memory space, remove the intermedied .bed .bim .fam files created
rm -r ./FilesCreated/*.bed; rm -r ./FilesCreated/*.bim; rm -r ./FilesCreated/*.fam
rm -r PCA_datasets.bed; rm -r PCA_datasets.bim; rm -r PCA_datasets.fam
rm -r HRC-1000G-check-bim-v4.3.0.zip
rm -r QCreferencedata.bed; rm -r QCreferencedata.bim; rm -r QCreferencedata.fam

echo "### CONGRATULATIONS! You've just succesfully completed the preimputation tutorial! You now must conduct the TOPMed imputation process."


###################################################################
### MICHIGAN IMPUTATION ###
echo "### MICHIGAN IMPUTATION ###"

echo "Load the 02_PreImputed_TargetData_.split.X.vcf.gz files to MICHIGAN Imputation Server using the following parameters:"
echo "# Reference panel: HRC r1.1 2016 (GRCh37/hg19)"
echo "# Array Build: GRCh37/hg19 (depending on the target base)"
echo "# rsq filter: 0.3"
echo "# Phasing: Eagle v2.4 (phased output)"
echo "# Population: EUR"
echo "# Mode: Quality Control & Imputation"

echo "### After the data is imputed, continue with the pipeline with the postimputation."
echo "### NOTE: The output of MICHIGAN can be decides to be under GRCh37. Thus, no liftover step is needed."
echo "### WARNING! The output of TOPMed is under GRCh38. Thus, target data need to be liftover."


###################################################################