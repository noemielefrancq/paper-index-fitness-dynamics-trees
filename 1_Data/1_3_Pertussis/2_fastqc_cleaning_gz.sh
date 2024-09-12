#########################################################
## QC and cleaning of fastq files
#########################################################

#!/bin/sh
##file must be fastq (paired)

## Unzip files
pigz -d -p 2 --keep fastq_raw/$( basename ${file} _1.fastq.gz)_1.fastq.gz ;
pigz -d -p 2 --keep fastq_raw/$( basename ${file} _1.fastq.gz)_2.fastq.gz ;
mv fastq_raw/$( basename ${file} _1.fastq.gz)_1.fastq fastq_raw/$( basename ${file} _1.fastq.gz)_2.fastq ./

## Check quality of fastq
~/softwares/FastQC/fastqc -t 2 $( basename ${file} _1.fastq.gz)_1.fastq $( basename ${file} _1.fastq.gz)_2.fastq;

## Deduplicate
for f in $( basename ${file} _1.fastq.gz)_1.fastq $( basename ${file} _1.fastq.gz)_2.fastq; do echo $f; done > $( basename ${file} _1.fastq.gz)_input.list
fastuniq -i $( basename ${file} _1.fastq.gz)_input.list -tq -o $( basename ${file} _1.fastq.gz)_1_dedup.fastq -p $( basename ${file} _1.fastq.gz)_2_dedup.fastq

## TRIM READS
## Quality
~/.local/bin/cutadapt -q 30,30 --pair-filter=any -o $( basename ${file} _1.fastq.gz)_1_trim_dedup.fastq -p $( basename ${file} _1.fastq.gz)_2_trim_dedup.fastq $( basename ${file} _1.fastq.gz)_1_dedup.fastq $( basename ${file} _1.fastq.gz)_2_dedup.fastq ;
## Cut 20 FIRST bases
~/.local/bin/cutadapt -u 20 -U 20 --pair-filter=any -o $( basename ${file} _1.fastq.gz)_1_trim_dedup_tmp.fastq -p $( basename ${file} _1.fastq.gz)_2_trim_dedup_tmp.fastq $( basename ${file} _1.fastq.gz)_1_trim_dedup.fastq $( basename ${file} _1.fastq.gz)_2_trim_dedup.fastq ;
mv $( basename ${file} _1.fastq.gz)_1_trim_dedup_tmp.fastq $( basename ${file} _1.fastq.gz)_1_trim_dedup.fastq;
mv $( basename ${file} _1.fastq.gz)_2_trim_dedup_tmp.fastq $( basename ${file} _1.fastq.gz)_2_trim_dedup.fastq;
## Quality again
~/.local/bin/cutadapt -q 30,30 --pair-filter=any -o $( basename ${file} _1.fastq.gz)_1_trim_dedup_tmp.fastq -p $( basename ${file} _1.fastq.gz)_2_trim_dedup_tmp.fastq $( basename ${file} _1.fastq.gz)_1_trim_dedup.fastq $( basename ${file} _1.fastq.gz)_2_trim_dedup.fastq ;
mv $( basename ${file} _1.fastq.gz)_1_trim_dedup_tmp.fastq $( basename ${file} _1.fastq.gz)_1_trim_dedup.fastq;
mv $( basename ${file} _1.fastq.gz)_2_trim_dedup_tmp.fastq $( basename ${file} _1.fastq.gz)_2_trim_dedup.fastq;
## Look for Nextera adapters
~/.local/bin/cutadapt -a CTGTCTCTTATA -g CTGTCTCTTATA --pair-filter=any -o $( basename ${file} _1.fastq.gz)_1_trim_dedup_tmp.fastq -p $( basename ${file} _1.fastq.gz)_2_trim_dedup_tmp.fastq $( basename ${file} _1.fastq.gz)_1_trim_dedup.fastq $( basename ${file} _1.fastq.gz)_2_trim_dedup.fastq ;
mv $( basename ${file} _1.fastq.gz)_1_trim_dedup_tmp.fastq $( basename ${file} _1.fastq.gz)_1_trim_dedup.fastq;
mv $( basename ${file} _1.fastq.gz)_2_trim_dedup_tmp.fastq $( basename ${file} _1.fastq.gz)_2_trim_dedup.fastq;
## Minimum length
~/.local/bin/cutadapt --minimum-length=20 --pair-filter=any -o $( basename ${file} _1.fastq.gz)_1_trim_dedup_tmp.fastq -p $( basename ${file} _1.fastq.gz)_2_trim_dedup_tmp.fastq $( basename ${file} _1.fastq.gz)_1_trim_dedup.fastq $( basename ${file} _1.fastq.gz)_2_trim_dedup.fastq ;
mv $( basename ${file} _1.fastq.gz)_1_trim_dedup_tmp.fastq $( basename ${file} _1.fastq.gz)_1_trim_dedup.fastq;
mv $( basename ${file} _1.fastq.gz)_2_trim_dedup_tmp.fastq $( basename ${file} _1.fastq.gz)_2_trim_dedup.fastq;

rm $( basename ${file} _1.fastq.gz)_1.fastq $( basename ${file} _1.fastq.gz)_2.fastq

## CHECK QUALITY FASTQC TRIMMED	
~/softwares/FastQC/fastqc -t 2 $( basename ${file} _1.fastq.gz)_1_trim_dedup.fastq $( basename ${file} _1.fastq.gz)_2_trim_dedup.fastq ;

## Clean directory and zip raw fastq files again
rm $( basename ${file} _1.fastq.gz)_1_dedup.fastq ;
rm $( basename ${file} _1.fastq.gz)_2_dedup.fastq ;
pigz -p 2 $( basename ${file} _1.fastq.gz)_1_trim_dedup.fastq ;
pigz -p 2 $( basename ${file} _1.fastq.gz)_2_trim_dedup.fastq ;

rm $( basename ${file} _1.fastq.gz)_input.list

mv $( basename ${file} _1.fastq.gz)_1_trim_dedup.fastq.gz $( basename ${file} _1.fastq.gz)_2_trim_dedup.fastq.gz fastq_trimmed_dedup/

mv $( basename ${file} _1.fastq.gz)*dedup_fastqc.* reports_trimmed/
mv $( basename ${file} _1.fastq.gz)*1_fastqc.* $( basename ${file} _1.fastq.gz)*2_fastqc.* reports_not_trimmed/

