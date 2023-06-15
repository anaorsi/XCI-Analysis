#!/bin/bash/
# step for SplitNCigarReads its not going when I merge all bam files from the same embryo (~20 cells...) 
# Example: sh 03.3_Allelic_Imbalance.sh Female_tri_Vol12_D12_IVC7_E1 HRR052798
##HRR052798 HRR052799 HRR052800

#Arguments --------------------------------------------------------
#Embryo name
embryo=$1
#embryo="Female_tri_Vol12_D12_IVC7_E1"

cell_accession=$2
#HRR052798 HRR052799 HRR052800

#create output dirs
#mkdir -p /home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/${embryo}/${cell_accession}

##MOD
cell_output="/home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/"${embryo}/${cell_accession}
embryo_output="/home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/"${embryo}
##cell_output="/home/anaorsi/data/new_ffEPSCs/00.Teste/"${embryo}/${cell_accession}
##embryo_output="/home/anaorsi/data/new_ffEPSCs/00.Teste/"${embryo}


## MOD
vcf_input=${embryo_output}/tmp/${embryo}_HTR_filtered.vcf.gz
##vcf_input=/home/anaorsi/genomas/X7f.vcf.gz


#set up external files---------------------------------------------
fasta=/home/anaorsi/genomas/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa

echo "--> ${embryo}: ${cell_accession} process started on:"
date

#1. add read group information per cell
# echo "1. adding read group information with Picard" 
#picard AddOrReplaceReadGroups \
#  -I ${embryo_output}/${embryo}_${cell_accession}_split.bam \
#  -O ${tmp_output}/${cell_accession}_reheader.bam \
# -RGLB ${cell_accession} \
# -RGPL illumina \
# -RGPU ${cell_accession} \
# -RGSM ${cell_accession} \
# -VERBOSITY ERROR -QUIET true
 
# then run the ASE read counter on individual cell BAMs at those heterozygous sites
# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php

#2. Read count
echo "2. Read count HTR sites" 
# -I ${tmp_output}${cell_accession}_reheader.bam \
/home/anaorsi/tools/gatk-4.2.5.0/gatk ASEReadCounter \
 -R $fasta \
 -I ${embryo_output}/${embryo}_${cell_accession}_split.bam \
 -V ${vcf_input} \
 -O ${cell_output}/${cell_accession}_HTR_ASEReadCount.txt \
 -verbosity ERROR

#Clean
#apagar so no final se tudo tiver dado certo
if test -f "${cell_output}/${cell_accession}_HTR_ASEReadCount.txt" ; then
    
  #  rm -r ${tmp_output}
    echo "--> process completed on:"
    date
else 
    echo "Something went wrong"
    echo "--> process completed on:"
    date

fi
 
 

