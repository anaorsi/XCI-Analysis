#!/usr/bin/bash

echo "eu funciono"

# 2021.01.06
# Quality check and aligment
# This script first makes use of Rui Wang quality control pipeline (https://github.com/WRui/Post_Implantation) to remove low-quality reads and adapter contaminated reads; 
# The filtered fastq is then aligend with HISAT2. 
# HISAT2 is the aligner used, and the outputted SAM file is converted to BAM, sorted by genomic coordinates and has duplicate reads removed using sambamba.
# USAGE: sh 02_HISAT2.sh <path to the fastq file> <path to the genome index dir>
# You have to give the whole path for the genome index
# Example: sh 02_HISAT2.sh HRR052691 ~/Documents/References/Human/GRCh38/HISAT2/grch38_snp_tran/genome

#Arguments --------------------------------------------------------
# sample name
fq=$1
# directory with genome index
genome=$2

# make all of the output directories ------------------------------
##Raw files
#mkdir -p ~/data/new_ffEPSCs/00.Raw_files/${fq}
##log and quality files
mkdir -p ~/logs/ffepscs_logs/${fq}
#trimmed data
mkdir -p ~/data/new_ffEPSCs/01.Clean_data/${fq}
#bam files
mkdir -p ~/data/new_ffEPSCs/02.HISAT2/${fq}
#count files
mkdir -p ~/data/new_ffEPSCs/03.HTSeq/${fq}

# Set up input required scripts locations--------------------------
QC=~/scripts/ffepscs_scripts/QC.pl


# Set up file locations--------------------------------------------
##log dir
log_dir=~/logs/ffepscs_logs/${fq}
##Raw fq to clean
raw_fq_in=~/data/new_ffEPSCs/00.Raw_files/${fq}
##Cleaned fq
clean_fq_out=~/data/new_ffEPSCs/01.Clean_data
##Cleaned fq to align
fastq_in=~/data/new_ffEPSCs/01.Clean_data/${fq}
##aligned bam files
aln_results=~/data/new_ffEPSCs/02.HISAT2/${fq}
##counted reads
read_count=~/data/new_ffEPSCs/03.HTSeq/${fq} 


#set up external files---------------------------------------------
## GTF for read count
gtf=~/genomas/hisat2/Homo_sapiens.GRCh38.84.gtf


## set up output file names----------------------------------------
align_log=${log_dir}/${fq}_HISAT.log

align_out=${aln_results}/${fq}_unsorted.sam
align_bam=${aln_results}/${fq}_unsorted.bam
align_sorted=${aln_results}/${fq}_sorted.bam
align_filtered=${aln_results}/${fq}_aln.bam
htseq_bam_out=${aln_results}/${fq}_aln.htseq.sort.bam
htseq_sam_out=${aln_results}/${fq}_aln.htseq.sort.sam

htseq_out=${read_count}/${fq}_htseq.out 

#echo "QC script: ${QC}"
#echo "log_dir: ${log_dir}"
echo "raw_fq_in: ${raw_fq_in}"
echo "fastq_in: ${fastq_in}"
#echo "aln_results: ${aln_results}"
#echo "read_count: ${read_count}"
#echo "gtf file path: ${gtf}"
#echo "align_log: ${align_log}"
#echo "align_out: $align_out"
#echo "align_bam: $align_bam"
#echo "align_sorted: $align_sorted"
#echo "align_filtered: $align_filtered"
#echo "htseq_bam_out: ${htseq_bam_out}"
#echo "htseq_out: ${htseq_out}"

thr=15

#----------------------- Start process ---------------------------#
echo "--> ${fq} process started on:"
date

# I have to download the raw files ........................
#echo "Download started"
#tools/sra-tools/sratoolkit.2.9.6-ubuntu64/bin/fastq-dump --outdir ${raw_fq_in} --skip-technical --dumpbase --split-3 --clip ${fq}
date

#Quality and trimming----------------------------------------------
echo "running QC_plus_rm_primer_polyA_T_trimTSO"
perl ${QC} --indir ${raw_fq_in} --outdir ${clean_fq_out} --sample ${fq} --end 2 --scRNA 1 
mv -t ${log_dir} ${clean_fq_out}/${fq}/*.log ${clean_fq_out}/${fq}/*.pdf ${clean_fq_out}/${fq}/*.ATGC ${clean_fq_out}/${fq}/*.mean_quality 
#rm ${raw_fq_in}/${fq}/*fastq
date


#Aligment----------------------------------------------------------
#State HISAT2 version
echo " 1. Alignment"
#echo "--> HISAT2 version:"
#/usr/local/bin/hisat2 --version

#Message indicating the chosen genome
echo "--> HISAT2 using $genome as genome index"

# Run HISAT2
#echo "--> command line: \n hisat2 \ \n -x $genome \ \n -U $fastq_in/${fq}.R1.clean.fq.gz,$fastq_in/${fq}.R2.clean.fq.gz \ \n -S $align_out  \ \n --summary-file $align_log  \ \n --time --threads 2 --seed 123 "
#fr-unstranded (default)
/usr/local/bin/hisat2 \
-x $genome \
--known-splicesite-infile /home/anaorsi/genomas/hisat2/splicesites84.txt \
-1 $fastq_in/${fq}.R1.clean.fq.gz -2 $fastq_in/${fq}.R2.clean.fq.gz \
-S $align_out \
--time \
--summary-file $align_log \
--threads $thr \
--seed 123 \
--rna-strandness RF --avoid-pseudogene

#-u 2000: align the first <int> reads or read pairs from the input (after the -s/--skip reads or pairs have been skipped), then stop. Default: no limit.
#mv -t "/media/jozuza/Backup Plus/Documents/GSE109555/HRA000128/Processing/01.Clean_data/" ${clean_fq_out}/${fq}

#Processing aligment files-----------------------------------------
# Create BAM from SAM
echo " 2. Create BAM from SAM"
~/tools/sambamba-0.8.1-linux-amd64-static view -h -S -t $thr -f bam $align_out > $align_bam
#/usr/local/bin/samtools view -h -S -b -o $align_bam $align_out
#     -h: include header 
#     -S: input format is auto-detected
#     -b: output BAM
#     -o: output file name
 
rm $align_out
 
# Sort BAM file by genomic coordinates
echo " 3. Sort BAM file by genomic coordinates"
~/tools/sambamba-0.8.1-linux-amd64-static sort -t $thr -o $align_sorted $align_bam
#    -q: quiet mode 
#    -t: number of threads / cores
#    -o: output file name

####rm $align_bam
# Filter out multi-mappers and duplicates
echo " 4. Keep only uniquely mapping reads"

#    -q: quiet mode 
#    -h: print SAM header before reads
#    -t: number of threads / cores
#    -f: format of output file (default is SAM)
#    -F: set custom filter - we will be using the filter to remove duplicates, multimappers and unmapped reads.
~/tools/sambamba-0.8.1-linux-amd64-static view -h -t $thr -f bam -F "not unmapped and not duplicate" $align_sorted > $align_filtered
#rm $align_sorted
#rm *sorted.bam.bai
# Create index for visualization
echo " 5. Create index"
/usr/local/bin/samtools index $align_filtered

# Remove intermediate bam files 
rm $aln_results/${fq}*sorted* 

#Read count--------------------------------------------------------
echo " 6. Read count"
date
#echo "--> command line:\n htseq-count --quiet --stranded no --format bam \ \n --samout $htseq_bam_out --samout-format BAM \ \n $align_filtered \ \n $gtf \ \n > $htseq_out"
/usr/local/bin/htseq-count --quiet --order pos --stranded reverse -m intersection-nonempty --format bam $align_filtered $gtf > $htseq_out
#-n: Number of parallel CPU processes to use 
#--minaqual: default 10
#o meu nao vai ter nada disso pq estou usando ensemble ID
#--type: default, suitable for Ensembl GTF files: exon
#--idattr : default, suitable for Ensembl GTF files: gene_id
#grep -v -P '^ERCC-|^RGC-|MIR|SNORD|Mir|Snord' $htseq_out > $htseq_clean

date

###/usr/local/bin/samtools view -S -b $htseq_sam_out > $htseq_bam_out

# Create index for visualization
###echo " 7. Create index"
###/usr/local/bin/samtools index $htseq_bam_out

# Move back to external hardrive........................
#echo "Copying file to external hard drive"
#mv -t "/media/jozuza/Backup Plus/Documents/GSE109555/HRA000128/Processing/Log_Quality_files/" ${log_dir}
#mv -t "/media/jozuza/Backup Plus/Documents/GSE109555/HRA000128/Processing/02.HISAT2/" ${aln_results}
#mv -t "/media/jozuza/Backup Plus/Documents/GSE109555/HRA000128/Processing/03.HTSeq/" ${read_count}

echo "--> process completed on:"
date


