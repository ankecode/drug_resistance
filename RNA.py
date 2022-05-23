#!/usr/bin/env python
# -*- coding: utf-8 -*-
__metaclass__ = type
# Author: An Ke

import re
import os
wd = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/rawdata/RNA_seq"
fqc_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/RNA_seq_new/1_fastqc"
fqc_path1 = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/RNA_seq_new/1_fastqc/before"
fqc_path2 = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/RNA_seq_new/1_fastqc/after"
cuta_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/RNA_seq_new/2_cutadapt"
trim_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/RNA_seq_new/3_trimmomatic"
hisat_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/RNA_seq_new/4_hisat2"
htseq_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/RNA_seq_new/5_htseq"
fpkm_tpm = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/RNA_seq_new/6_count_fpkm_tmp"


#os.mkdir(fqc_path)
#os.mkdir(fqc_path1)
#os.mkdir(fqc_path2)
#os.mkdir(cuta_path)
#os.mkdir(trim_path)
#os.mkdir(hisat_path)
#os.mkdir(htseq_path)
#os.mkdir(fpkm_tpm)
a = os.walk(wd)
next(a)
for i in a:
	outputdir1 = os.path.join(fqc_path1,os.path.basename(i[0]))
	outputdir2=os.path.join(fqc_path2,os.path.basename(i[0]))
	outputdir3=os.path.join(cuta_path,os.path.basename(i[0]))
	outputdir4=os.path.join(trim_path,os.path.basename(i[0]))
	outputdir5=os.path.join(hisat_path,os.path.basename(i[0]))
	outputdir6=os.path.join(htseq_path,os.path.basename(i[0]))
	outputdir7=os.path.join(fpkm_tpm,os.path.basename(i[0]))
	
#	os.mkdir(outputdir1)
#	os.mkdir(outputdir2)
#	os.mkdir(outputdir3)
#	os.mkdir(outputdir4)
#	os.mkdir(outputdir5)
#	os.mkdir(outputdir6)
#	os.mkdir(outputdir7)
	
	for fq in i[2]:
		if re.search(r"1.fq.gz",fq):
			fq1=os.path.join(i[0],fq)
		if re.search(r"2.fq.gz",fq):
			fq2=os.path.join(i[0],fq)
			sample_name=os.path.basename(i[0])
			cutadapt1=os.path.join(outputdir3,os.path.basename(i[0])+"_cuta_1.fq.gz")
			cutadapt2=os.path.join(outputdir3,os.path.basename(i[0])+"_cuta_2.fq.gz")
			trim1=os.path.join(outputdir4,os.path.basename(i[0])+"_trim_1.fq.gz")
			trim1_u=os.path.join(outputdir4,os.path.basename(i[0])+"_unpair_1.fq.gz")
			trim2=os.path.join(outputdir4,os.path.basename(i[0])+"_trim_2.fq.gz")
			trim2_u=os.path.join(outputdir4,os.path.basename(i[0])+"_unpair_2.fq.gz")
	script='''#!/bin/bash
#SBATCH -J %(sample_name)s
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH -p normal
#SBATCH -o %(sample_name)s_3.out
#SBATCH -e %(sample_name)s_3.err
#SBATCH -t 10-00:00:00


cd %(outputdir1)s
fastqc -o %(outputdir1)s %(fq1)s %(fq2)s
#
cd %(outputdir3)s
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -o %(cutadapt1)s -p %(cutadapt2)s %(fq1)s %(fq2)s
#
cd %(outputdir4)s
trimmomatic PE -phred33 %(cutadapt1)s %(cutadapt2)s %(trim1)s %(trim1_u)s %(trim2)s %(trim2_u)s LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:35
#
cd %(outputdir2)s
fastqc -o %(outputdir2)s %(trim1)s %(trim2)s
#
cd %(outputdir5)s
hisat2  -p 8 -N 1 --dta -x /public/home/zdyfy03/reference/human/genecode_rel38_hg38/Index/hisat2/GRCh38.primary_assembly.genome.fa --rna-strandness RF -1 %(trim1)s -2 %(trim2)s -S %(sample_name)s.sam --un-conc ./
#
gzip un-conc-mate.1
gzip un-conc-mate.2
#
echo "reads with no -q:" > %(sample_name)s.readsCount
cat %(sample_name)s.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> %(sample_name)s.readsCount
samtools view -h -Sb %(sample_name)s.sam > %(sample_name)s.bam
samtools view %(sample_name)s.bam -q 20 -h > %(sample_name)s_20.sam
#
echo "reads with -q20:" >> %(sample_name)s.readsCount
cat %(sample_name)s_20.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> %(sample_name)s.readsCount
#rm %(sample_name)s_20.sam
#
samtools view -S %(sample_name)s_20.sam -b -o %(sample_name)s.uniqmap.bam
samtools sort %(sample_name)s.uniqmap.bam -@ 8 -o %(sample_name)s.uniqmap_sort.bam
samtools index %(sample_name)s.uniqmap_sort.bam
bamToBed -i %(sample_name)s.uniqmap_sort.bam -split > %(sample_name)s.uniqmap_sort.bed
#

htseq-count -m union -f bam -s no %(outputdir5)s/%(sample_name)s.uniqmap_sort.bam /public/home/zdyfy03/reference/human/genecode_rel38_hg38/gencode.v38.primary_assembly.annotation.gtf > %(sample_name)s_count_union_no.out

grep -v "__" %(sample_name)s_count_union_no.out > %(sample_name)s_count_union_no_del_tail5.out
grep -v "PAR_Y" %(sample_name)s_count_union_no_del_tail5.out|awk 'BEGIN{OFS="\\t"} {gsub(/[:.:].+/,"",$1); print $0}' > %(sample_name)s_count_union_no_del_tail5_final.out

cd %(outputdir7)s
Rscript /public/home/zdyfy03/script/R/RPKM_and_TPM_human_genecode.R  %(outputdir6)s/%(sample_name)s_count_union_no_del_tail5.out %(sample_name)s

grep -v 'Count' %(sample_name)s_count_rpkm_tpm_tmp.txt | sed '1i\Gene_id\\t%(sample_name)s.Count\\t%(sample_name)s.rpkm\\t%(sample_name)s.TPM'  > %(sample_name)s_count_rpkm_tpm

cut -f 1,2 %(sample_name)s_count_rpkm_tpm | grep -v 'Gene_id' | sed '1i\Gene_id\\t%(sample_name)s' > %(sample_name)s_count
cut -f 1,3 %(sample_name)s_count_rpkm_tpm | grep -v 'Gene_id' | sed '1i\Gene_id\\t%(sample_name)s' > %(sample_name)s_rpkm
cut -f 1,4 %(sample_name)s_count_rpkm_tpm | grep -v 'Gene_id' | sed '1i\Gene_id\\t%(sample_name)s' > %(sample_name)s_tpm
rm %(sample_name)s_count_rpkm_tpm_tmp.txt

	''' %(locals())
	open(sample_name+"_RNA.sh","w").write(script)
	#os.system("qsub "+sample_name+"_RNA.sh")
