#########################################################
## BWA mapping + Variant calling with GATK
#########################################################

#!/bin/sh

#################################################### Mapping with BWA
#Create pe.sam files out of trimmed fastq
## Decompress trimmed fastqs
pigz -p 4 -d --keep fastq_trimmed_dedup/$( basename ${file} _1_trim_dedup.fastq.gz)_1_trim_dedup.fastq.gz ;
pigz -p 4 -d --keep fastq_trimmed_dedup/$( basename ${file} _1_trim_dedup.fastq.gz)_2_trim_dedup.fastq.gz ;

## Move files to current directory
mv fastq_trimmed_dedup/$( basename ${file} _1_trim_dedup.fastq.gz)_1_trim_dedup.fastq fastq_trimmed_dedup/$( basename ${file} _1_trim_dedup.fastq.gz)_2_trim_dedup.fastq ./

## Run bwa
R1=$( basename ${file} _1_trim_dedup.fastq.gz)""_1_trim_dedup.fastq
R2=$( basename ${file} _1_trim_dedup.fastq.gz)""_2_trim_dedup.fastq
/home/ncmjl2/softwares/bwa/bwa mem -t 4 /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta $R1 $R2 > $( basename ${file} _1_trim_dedup.fastq.gz)_trim.pe.sam

## Clean directory
rm $R1 ;
rm $R2 ;

##################################################### Variant calling
## Sort sam file
/home/ncmjl2/softwares/gatk-4.2.0.0/gatk SortSam -I $( basename ${file} _1_trim_dedup.fastq.gz)_trim.pe.sam -O $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_sorted.bam -SO coordinate
rm $( basename ${file} _1_trim_dedup.fastq.gz)_trim.pe.sam

## Mark duplicates, add read group and create a bam index
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk MarkDuplicates -I $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_sorted.bam -O $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup.bam -M $( basename ${file} _1_trim_dedup.fastq.gz)""metrics.txt 
java -jar /home/ncmjl2/softwares/picard/build/libs/picard.jar AddOrReplaceReadGroups I=$( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup.bam O=$( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$( basename ${file} _1_trim_dedup.fastq.gz)
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk BuildBamIndex -I $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bam

## Get first variant calling done
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk HaplotypeCaller -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta -I $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bam --base-quality-score-threshold 20 --minimum-mapping-quality 60 --mapping-quality-threshold-for-genotyping 60 -ploidy 1 -ERC GVCF --annotation StrandBiasBySample --annotation AlleleFraction -G AS_StandardAnnotation -O $( basename ${file} _1_trim_dedup.fastq.gz)_GATK.g.vcf --all-site-pls

## Compute coverage
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk CollectWgsMetrics --USE_FAST_ALGORITHM -CAP 500 --MINIMUM_BASE_QUALITY 20 --MINIMUM_MAPPING_QUALITY 60 --INCLUDE_BQ_HISTOGRAM -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta -I $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bam -O $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage.txt 

## Clean
rm $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup.bam 
rm $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_sorted.bam
rm $( basename ${file} _1_trim_dedup.fastq.gz)""metrics.txt 
rm $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage_sd.txt
rm $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage_mean.txt
rm $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage_line.txt
rm -r tmp_$( basename ${file} _1_trim_dedup.fastq.gz)

## Move GVCF files - directories must have been created beforehand
mv $( basename ${file} _1_trim_dedup.fastq.gz)""_GATK.g.vcf g_vcf_IS_masked/
mv $( basename ${file} _1_trim_dedup.fastq.gz)""_GATK.g.vcf.idx g_vcf_IS_masked/
mv $( basename ${file} _1_trim_dedup.fastq.gz)""_coverage.txt g_vcf_IS_masked/

##################################################### Filter variants, to take out IS
## Min coverage to call a SNP
threshold=5

## GATK Genotype
/home/ngmjl2/softwares/gatk-4.2.0.0/./gatk --java-options "-Xmx10g -Xms10g" GenotypeGVCFs \
   -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta \
   -V g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_GATK.g.vcf \
   --sample-ploidy 1 \
   -G AS_StandardAnnotation \
   -A AlleleFraction \
   -O g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_genotyped.vcf.gz
      
## Filter IS and SNPs - including mutations in IS elements
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk --java-options "-Xmx10g -Xms10g" VariantFiltration \
   -V g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_genotyped.vcf.gz \
   -filter "QD < 2.0" --filter-name "QD2" \
   -filter "MQ < 20.0" --filter-name "MQ20" \
   -filter "DP < $threshold" --filter-name "DPthre" \
   -filter "SOR > 5.0"  --filter-name "SORthre" \
   --genotype-filter-expression "AF < 0.75" --genotype-filter-name "AFthre" \
   --exclude-intervals /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/IS_all_Tohama_pertussis_full_positions_transposase_inc.intervals \
   -O g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.vcf.gz

## Excluse filtered IS and SNPs
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk SelectVariants \
    -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta \
    -V g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.vcf.gz \
    -O g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.vcf.gz \
    --exclude-filtered \
    --set-filtered-gt-to-nocall \
    --select-type-to-include SNP
    
bgzip -d g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.vcf.gz
cat g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.vcf | grep -v AFthre > g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS_tmp.vcf
mv g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS_tmp.vcf g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.vcf 
bgzip g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.vcf

##################################################### Compute mask, to put N is the reference
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk HaplotypeCaller -base-quality-score-threshold 20 --minimum-mapping-quality 60 --mapping-quality-threshold-for-genotyping 60 -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta -I $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bam -ploidy 1 -ERC BP_RESOLUTION  -G AS_StandardAnnotation  --annotation StrandBiasBySample -O g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_GATK_allsite.g.vcf

/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk GenotypeGVCFs \
    -V g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_GATK_allsite.g.vcf \
    -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta \
    -G AS_StandardAnnotation \
    --all-sites \
    -O g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_allsite_raw.vcf
    
/home/ncmjl2/softwares/gatk-4.2.0.0/./gatk SelectVariants \
    -R /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta \
    -V g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_allsite_raw.vcf \
    -O g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked1.vcf \
    -select "DP < $threshold || SOR > 5.0 || QD < 2.0 || MQ < 20.0" 
    
grep :0:0:0 g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_allsite_raw.vcf > g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked2.vcf
    
cat g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked1.vcf g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked2.vcf > g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked.vcf

rm g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked1.vcf 
rm g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked2.vcf

cat g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked.vcf |tail -n +36 |awk '{FS="\t";OFS="\t";print $1,$2-1,$2,$3, etc}' > g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked.bed

##################################################### Create pseudo consensus fasta file
/home/ncmjl2/softwares/bcftools-1.12/bcftools consensus -f /home/ncmjl2/rds/rds-hs743-arbodynamic/Noemie/Pertussis_fastq/Ref_seq/TohamaNC_0029292.fasta --mask g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked.bed g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.vcf.gz > g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_consensus.fasta

/home/ncmjl2/softwares/bcftools-1.12/bcftools annotate -x INFO,^FORMAT/GT -O z g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.vcf.gz > g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.GT.vcf.gz


##################################################### Clean direcotory 
rm g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_genotyped.vcf.gz* 
rm g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.vcf.gz*
mv g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.vcf.gz* filtered_PASS_vcf/
mv g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_masked.* filtered_PASS_vcf/
rm g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_GATK_allsite*
rm g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)*.idx
rm g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_allsite*
mv g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_consensus* consensus_fasta/
mv g_vcf_IS_masked/$( basename ${file} _1_trim_dedup.fastq.gz)_filtered.PASS.GT.vcf.gz filtered_PASS_no_GT_vcf/
mv $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bam bam_dedup/
mv $( basename ${file} _1_trim_dedup.fastq.gz)""_trim_dedup2.bai bam_dedup/
