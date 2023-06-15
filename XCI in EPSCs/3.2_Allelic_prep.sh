#!/bin/bash/

### A FAZER:
# OK Conseguir esse arquivo do gatk: resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.Fixed.vcf.gz


#Now you run 
# Example: sh 03.2_Allelic_Imbalance_prep.sh Female_tri_Vol12_D12_IVC7_E1 

#Arguments --------------------------------------------------------
#Embryo name
embryo=$1
#embryo="Female_tri_Vol12_D12_IVC7_E1"
thr=16

#create output dirs
###mkdir -p /home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/${embryo}/tmp
tmp_output="/home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/"${embryo}/"tmp/"
embryo_output="/home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/"${embryo}

#set up external files---------------------------------------------
vcf=/home/anaorsi/genomas/new.vcf.gz
fasta=/home/anaorsi/genomas/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa


echo "--> ${embryo} process started on:"
date

#1. merge single-cell alignments from the same embryo
###echo "1. merging bam files and sorting with samtools"
###ls "/home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/"${embryo}/${embryo}*_split.bam > ${tmp_output}${embryo}_bam_path.txt
#merge all cells from the same embryo
###/usr/local/bin/samtools merge \
###  -O BAM \
###  ${tmp_output}/${embryo}_merged.bam \
### -b ${tmp_output}${embryo}_bam_path.txt \
### -@ $thr

###/usr/local/bin/samtools sort \
### -O BAM \
###  -o ${tmp_output}/${embryo}_sorted.bam \
###  ${tmp_output}/${embryo}_merged.bam \
###  -@ $thr

###/usr/local/bin/samtools index ${tmp_output}/${embryo}_sorted.bam

###rm ${tmp_output}/${embryo}_merged.bam


###echo "3. Recalibrate base quality scores"
#3.1. Recalibrate base quality scores in order to correct sequencing errors and other experimental artifacts.
###echo "	3.1. BaseRecalibrator" 
###/home/anaorsi/tools/gatk-4.2.5.0/gatk BaseRecalibrator \
### -R $fasta \
### -I ${tmp_output}/${embryo}_sorted.bam \
### --known-sites $vcf \
### -O ${tmp_output}/${embryo}_recal.grp \
### -verbosity ERROR

#3.2. ApplyBQSR: BaseRecalibrator with -BQSR option
###echo "	3.2. Apply base quality score" 
###/home/anaorsi/tools/gatk-4.2.5.0/gatk ApplyBQSR \
### -R $fasta \
### -I ${tmp_output}/${embryo}_sorted.bam \
### --bqsr-recal-file ${tmp_output}/${embryo}_recal.grp \
### -O ${tmp_output}/${embryo}.realn_Recal.bam \
### -verbosity ERROR

##rm ${tmp_output}/${embryo}_recal.grp
##rm ${tmp_output}/${embryo}_realignedBam.bai
##rm ${tmp_output}/${embryo}_realignedBam.bam

#4. HaplotypeCaller 
##echo "4. HaplotypeCaller" 
##/home/anaorsi/tools/gatk-4.2.5.0/gatk HaplotypeCaller \
##  -R $fasta \
##  -I ${tmp_output}/${embryo}.realn_Recal.bam \
##  --dont-use-soft-clipped-bases true \
##  --standard-min-confidence-threshold-for-calling 20.0 \
##  -L X \
##  --native-pair-hmm-threads $thr \
##  -O ${tmp_output}/${embryo}_X.vcf \
##  -verbosity ERROR

bgzip ${tmp_output}/${embryo}_all.vcf
tabix -p vcf ${tmp_output}/${embryo}_all.vcf.gz

#9. limit to heterozygous sites
echo "5. limit to heterozygous sites"
/home/anaorsi/tools/gatk-4.2.5.0/gatk SelectVariants \
  -R $fasta \
  -V ${tmp_output}/${embryo}_all.vcf.gz \
  -O ${tmp_output}/${embryo}_HTR_SNP.vcf \
  -select-type SNP \
  -select "vc.getGenotype(\""${embryo}"\").isHet()"  \
  -verbosity ERROR

bgzip ${tmp_output}/${embryo}_HTR_SNP.vcf
tabix -p vcf ${tmp_output}/${embryo}_HTR_SNP.vcf.gz

/home/anaorsi/tools/gatk-4.2.5.0/gatk \
  VariantFiltration \
  -R $fasta \
  -V ${tmp_output}/${embryo}_HTR_SNP.vcf.gz \
  -O ${tmp_output}/${embryo}_HTR_filtered.vcf \
  -window 35 \
  -cluster 3 \
  --filter-name FS -filter "FS > 30.0" \
  --filter-name QD -filter "QD < 2.0" \
  -verbosity ERROR
  
bgzip ${tmp_output}/${embryo}_HTR_filtered.vcf 
tabix -p vcf ${tmp_output}/${embryo}_HTR_filtered.vcf.gz

##if test -f "${tmp_output}/${embryo}.realn_Recal.bam" ; then

##   mv ${tmp_output}${embryo}*vcf* ${embryo_output} 
##    mv ${tmp_output}${embryo}*realn_Recal* ${embryo_output}   
#    rm -r ${tmp_output}    
##    echo "--> process completed on:"
##    date
##else 
##    echo "Something went wrong"
##    echo "--> process completed on:"
##    date
##fi


