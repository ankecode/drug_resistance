#!/usr/bin/env python
# -*- coding: utf-8 -*-
__metaclass__ = type
# Author: An Ke

import re
import os
wd = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/rawdata/Ribo_seq"
fqc_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/Ribo_seq_new/1_fastqc"
fqc_path1 = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/Ribo_seq_new/1_fastqc/before"
fqc_path2 = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/Ribo_seq_new/1_fastqc/after"
cuta_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/Ribo_seq_new/2_cutadapt"
trim_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/Ribo_seq_new/3_trimmomatic"
hisat_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/Ribo_seq_new/4_bowtie"
htseq_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/Ribo_seq_new/5_htseq"
fpkm_tpm = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/Ribo_seq_new/6_cufflink_fpkm"


os.mkdir(fqc_path)
os.mkdir(fqc_path1)
os.mkdir(fqc_path2)
os.mkdir(cuta_path)
os.mkdir(trim_path)
os.mkdir(hisat_path)
os.mkdir(htseq_path)
os.mkdir(fpkm_tpm)
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
	
	os.mkdir(outputdir1)
	os.mkdir(outputdir2)
	os.mkdir(outputdir3)
	os.mkdir(outputdir4)
	os.mkdir(outputdir5)
	os.mkdir(outputdir6)
	os.mkdir(outputdir7)
	
	for fq in i[2]:
		if re.search(r"1.fq.gz",fq):
			fq1=os.path.join(i[0],fq)
		if re.search(r"2.fq.gz",fq):
			fq2=os.path.join(i[0],fq)
			sample_name=os.path.basename(i[0])
			cutadapt1=os.path.join(outputdir3,os.path.basename(i[0])+"_cuta_1.fq.gz")
			cutadapt1_1=os.path.join(outputdir3,os.path.basename(i[0])+"_cuta_1_1.fq.gz")
			cutadapt2=os.path.join(outputdir3,os.path.basename(i[0])+"_cuta_2.fq.gz")
			trim1=os.path.join(outputdir4,os.path.basename(i[0])+"_trim_1.fq.gz")
			trim1_u=os.path.join(outputdir4,os.path.basename(i[0])+"_unpair_1.fq.gz")
			trim2=os.path.join(outputdir4,os.path.basename(i[0])+"_trim_2.fq.gz")
			trim2_u=os.path.join(outputdir4,os.path.basename(i[0])+"_unpair_2.fq.gz")
	script='''#!/bin/bash
#SBATCH -J %(sample_name)s
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH -p normal
#SBATCH -o %(sample_name)s_3.out
#SBATCH -e %(sample_name)s_3.err
#SBATCH -t 20-00:00:00


cd %(outputdir1)s
fastqc -o %(outputdir1)s %(fq1)s %(fq2)s
#
cd %(outputdir3)s
cutadapt -a AGATCGGAAGAG -o %(cutadapt1)s %(fq1)s
cutadapt -a "A{10}" -u 3 -e 0.1 -m 18 -M 60 %(cutadapt1)s|gzip - > %(cutadapt1_1)s

##
cd %(outputdir4)s
trimmomatic SE -phred33 %(cutadapt1_1)s %(trim1)s LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:18
##
cd %(outputdir2)s
fastqc -o %(outputdir2)s %(cutadapt1_1)s
fastqc -o %(outputdir2)s %(trim1)s %(trim2)s
#
cd %(outputdir5)s
bowtie -p 32 --chunkmbs 1000 /public/home/zdyfy03/reference/human/rRNA/bowtie/human_all_rRNA.fa -q %(trim1)s -5 3 --un %(sample_name)s_un-rRNA.fq -S %(sample_name)s.rRNA.sam
samtools view -Sb -h -@ 32 %(sample_name)s.rRNA.sam > %(sample_name)s.rRNA.bam
rm %(sample_name)s.rRNA.sam
#
bowtie -p 32 -x /public/home/zdyfy03/reference/human/genecode_rel38_hg38/Index/bowtie/GRCh38.primary_assembly.genome.fa -q %(sample_name)s_un-rRNA.fq -S %(sample_name)s_un-rRNA.sam
#
samtools view -Sb -h -@ 32 %(sample_name)s_un-rRNA.sam > %(sample_name)s_un-rRNA.bam
#
samtools view %(sample_name)s_un-rRNA.bam -q 20 -h > %(sample_name)s_un-rRNA_20.sam
samtools view -Sb -h -@ 32 %(sample_name)s_un-rRNA_20.sam > %(sample_name)s_un-rRNA_20.bam
samtools sort -n %(sample_name)s_un-rRNA_20.bam -o %(sample_name)s_un-rRNA_20_sorted.bam
bedtools bamtobed -split -i %(sample_name)s_un-rRNA_20_sorted.bam > %(sample_name)s_un-rRNA_20_sorted.bed
#由于cufflink的原因，进行重新的bam获取
samtools view -H %(sample_name)s_un-rRNA_20.bam > header.sam
samtools reheader header.sam %(sample_name)s_un-rRNA_20.bam > %(sample_name)s_un-rRNA_20.reheader.bam
samtools sort %(sample_name)s_un-rRNA_20.reheader.bam > %(sample_name)s_un-rRNA_20.reheader.sort.bam

#
cd %(outputdir7)s
/public/home/zdyfy03/software/cufflinks-2.2.1.Linux_x86_64/cufflinks -p 32 -G /public/home/zdyfy03/reference/human/genecode_rel38_hg38/gencode.v38.primary_assembly.annotation.gtf -u --library-type fr-unstranded -o ./ %(outputdir5)s/%(sample_name)s_un-rRNA_20.reheader.sort.sam


	''' %(locals())
	open(sample_name+"_RNA_new_1.sh","w").write(script)
	#os.system("qsub "+sample_name+"_RNA.sh")
