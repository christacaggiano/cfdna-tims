
#!/bin/bash
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=32G
#$ -l h_rt=6:00:00
#$ -t 1-96

R1="NcfT_C"$SGE_TASK_ID"_S"$SGE_TASK_ID"_L004_R1_001.fastq.gz"
R2="NcfT_C"$SGE_TASK_ID"_S"$SGE_TASK_ID"_L004_R3_001.fastq.gz"
I1="NcfT_C"$SGE_TASK_ID"_S"$SGE_TASK_ID"_L004_R2_001.fastq.gz"
OUTPUT="mapped_bsbolt/"$SGE_TASK_ID
METH="bsbolt_meth/"$SGE_TASK_ID
TEMP="als_v1/temp/"$SGE_TASK_ID

##########################################################################
### QC on raw fastq files  

fastqc $R1 $R2 --outdir $OUTPUT


##########################################################################
#### add UMIs to fastq headers 

conda activate umi 
umi_tools extract --bc-pattern=NNNNNNNNN --stdin $I1 --read2-in $R1 --stdout $TEMP"_R1_umi.fastq.gz" --read2-stdout --log $TEMP"_umi.log"
umi_tools extract --bc-pattern=NNNNNNNNN --stdin $I1 --read2-in $R2 --stdout $TEMP"_R2_umi.fastq.gz" --read2-stdout --log $TEMP"_umi.log"


##########################################################################
#### perform trimming 
conda activate methyldackel

trim_galore --paired $TEMP"_R1_umi.fastq.gz" $TEMP"_R2_umi.fastq.gz" -o $TEMP"_trimmed" 

##########################################################################
#### align to bisulfite converted genome 

python -m bsbolt Align -t 4 -OT 4 -DB bsbolt_db -F1 $TEMP"_trimmed/"$SGE_TASK_ID"_R1_umi_val_1.fq.gz" \
	-F2 $TEMP"_trimmed/"$SGE_TASK_ID"_R2_umi_val_2.fq.gz"  -O $OUTPUT"_trimmed" > $OUTPUT"_log.txt"

##########################################################################
#### prepare for duplicate removal 

samtools fixmate -p -m $OUTPUT".bam" $OUTPUT".fixmates.bam" 
samtools sort -@ 4 -o $OUTPUT".sorted.bam" $OUTPUT".fixmates.bam"
samtools index $OUTPUT".sorted.bam"

##########################################################################
#### remove duplicates using UMIs

conda activate umi 
umi_tools dedup -I $OUTPUT".sorted.bam" --paired -S $OUTPUT".umi_dedup.bam"

samtools flagstat -@ 4 $OUTPUT".umi_dedup.bam"
samtools index -@ 4 $OUTPUT".umi_dedup.bam"


##########################################################################

#### call methylation 
conda activate methyldackel

python -m bsbolt CallMethylation -BG -remove-ccgg -I $OUTPUT".umi_dedup.bam" -DB bsbolt_db -O $METH"_test" 






