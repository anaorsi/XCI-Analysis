#!/usr/bin/bash

# A FAZER:
# OK conseguir esse arquivo do GATK: resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.Fixed.vcf
# aparentemente OKconseguir esse arquivo: Homo_sapiens.GRCh38.dna.primary_assembly.chrsizes.bed
# OK seguir as instruções: https://gatk.broadinstitute.org/hc/en-us/articles/360041320571--How-to-Install-all-software-packages-required-to-follow-the-GATK-Best-Practices
# picard, gatk, gatk3??

# first prepare all cells from the same embryo
# step for SplitNCigarReads its not going when I merge all bam files from the same embryo (~20 cells...) 
# Example: sh 03.1_Allelic_Imbalance_prep.sh Female_tri_Vol12_D12_IVC7_E1 HRR052798
##HRR052798 HRR052799 HRR052800

#Arguments --------------------------------------------------------
#Embryo name
embryo=$1
#embryo="Female_tri_Vol12_D12_IVC7_E1"

#cell_accession=$2
cell_accession=$2
#HRR052798 HRR052799 HRR052800

#create output dirs
#mkdir -p "/home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/"${embryo}
mkdir -p "/home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/"${embryo}/${cell_accession}"/tmp"
rm -r "/home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/"${embryo}/${cell_accession}"/tmp"
mkdir "/home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/"${embryo}/${cell_accession}"/tmp"


#set up variables
##input
bam_input="/home/anaorsi/data/new_ffEPSCs/02.HISAT2/"${cell_accession}"/"${cell_accession}"_aln.bam"
cell_input="/home/anaorsi/data/new_ffEPSCs/02.HISAT2/"${cell_accession}""
vcf_input="/home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/"${embryo}"_HTR_filtered.vcf.gz"

#output
embryo_output="/home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/"${embryo}
cell_output="/home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/"${embryo}/${cell_accession}/
tmp_output="/home/anaorsi/data/new_ffEPSCs/04.SNP_Calling/"${embryo}/${cell_accession}/"tmp/"

#set up external files---------------------------------------------
fasta=/home/anaorsi/genomas/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
vcf=/home/anaorsi/genomas/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf

thr=16

echo "cell_accession" ${cell_accession}
echo "cell_input" ${cell_input}
echo "tmp_output" ${tmp_output}
echo "cell_output" ${cell_output}
echo "bam_input" ${bam_input}

echo "--> ${cell_accession} process started on:"
date

#1. remove reads not mapped on main chromosomes
echo "1. remove reads not mapped on main chromosomes" 
/usr/local/bin/samtools view -bh \
-L ~/genomas/hisat2/Homo_sapiens.GRCh38.84.gtf.bed12 \
${bam_input} > ${tmp_output}${cell_accession}_temp.bam

#/usr/local/bin/samtools view -H \
#${tmp_output}${cell_accession}_temp.bam |  grep -v "SN:G" | grep -v "SN:K" > ${tmp_output}/newHeader.txt

####/usr/local/bin/samtools reheader \
####${tmp_output}/newHeader.txt ${tmp_output}${cell_accession}_temp.bam > ${tmp_output}${cell_accession}_MainChrs.bam

#rm ${tmp_output}${cell_accession}_temp.bam
#rm ${tmp_output}/newHeader.txt


echo "2. samtools sort"
/usr/local/bin/samtools sort \
  -O BAM \
  -o ${tmp_output}${cell_accession}_sorted.bam \
  ${tmp_output}${cell_accession}_temp.bam
 #### ${tmp_output}${cell_accession}_MainChrs.bam
  
#3. add read group information per embryo
echo "3. adding read group information with Picard"
#TODO: A Rui coloca RGID=$sp
java -jar /home/anaorsi/tools/picard.jar AddOrReplaceReadGroups \
  -I ${tmp_output}${cell_accession}_sorted.bam \
  -O ${tmp_output}${embryo}_${cell_accession}_reheader.bam \
  -RGLB ${embryo} \
  -RGPL illumina \
  -RGPU ${embryo} \
  -RGSM ${embryo} \
  -VERBOSITY ERROR -QUIET true
  
#rm ${tmp_output}${cell_accession}_MainChrs.bam

#4. mark duplicates and create index 
echo "4. marking duplicates with Picard" 
java -jar /home/anaorsi/tools/picard.jar MarkDuplicates \
  -I ${tmp_output}${embryo}_${cell_accession}_reheader.bam  \
  -O ${tmp_output}${embryo}_${cell_accession}_dedupped.bam  \
  -CREATE_INDEX true \
  -VALIDATION_STRINGENCY SILENT \
  -M ${tmp_output}${embryo}_${cell_accession}.metrics \
  -VERBOSITY ERROR -QUIET true

/usr/local/bin/samtools index ${tmp_output}${embryo}_${cell_accession}_dedupped.bam
#rm ${tmp_output}${embryo}_${cell_accession}_reheader.bam

#5. "split'n'trim" the reads and adjust mapping quality (no need to adjust quality)
echo "5. split'n'trim the reads" 

/home/anaorsi/tools/gatk-4.2.5.0/gatk  SplitNCigarReads \
 -R $fasta \
 -I ${tmp_output}${embryo}_${cell_accession}_dedupped.bam \
 -O ${tmp_output}${embryo}_${cell_accession}_split.bam

if test -f "${tmp_output}${embryo}_${cell_accession}_split.bam" ; then
    mv ${tmp_output}${cell_accession}_sorted.bam ${cell_output}
    mv ${tmp_output}${embryo}_${cell_accession}_split.* ${embryo_output}
    rm -r ${tmp_output}    
    echo "--> process completed on:"
    date
else 
    echo "Something went wrong"
    echo "--> process completed on:"
    date

fi

