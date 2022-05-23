#!/usr/bin/env python
# -*- coding: utf-8 -*-
__metaclass__ = type
# Author: An Ke

import re
import os
wd = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/rawdata/20210610_paris"
fqc_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/20210610_paris_new/1_fastqc"
fqc_path1 = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/20210610_paris_new/1_fastqc/before"
fqc_path2 = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/20210610_paris_new/1_fastqc/after"
cuta_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/20210610_paris_new/2_cutadapt"
trim_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/20210610_paris_new/3_trimmomatic"
hisat_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/20210610_paris_new/4_star"
gap_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/20210610_paris_new/5_gap_sam"
dg_path = "/public/home/zdyfy03/Data/shiby/20210515_naiyao_structure/analysis/20210610_paris_new/6_dg_cluster"


os.mkdir(fqc_path)
os.mkdir(fqc_path1)
os.mkdir(fqc_path2)
os.mkdir(cuta_path)
os.mkdir(trim_path)
os.mkdir(hisat_path)
os.mkdir(gap_path)
os.mkdir(dg_path)
a = os.walk(wd)
next(a)
for i in a:
	outputdir1 = os.path.join(fqc_path1,os.path.basename(i[0]))
	outputdir2=os.path.join(fqc_path2,os.path.basename(i[0]))
	outputdir3=os.path.join(cuta_path,os.path.basename(i[0]))
	outputdir4=os.path.join(trim_path,os.path.basename(i[0]))
	outputdir5=os.path.join(hisat_path,os.path.basename(i[0]))
	outputdir6=os.path.join(gap_path,os.path.basename(i[0]))
	outputdir7=os.path.join(dg_path,os.path.basename(i[0]))
	
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
			cutadapt1_1=os.path.join(outputdir3,os.path.basename(i[0])+"_cuta_1_1.fq")
			cutadapt_dup=os.path.join(outputdir3,os.path.basename(i[0])+"_cuta_dup.fq")
			cutadapt2=os.path.join(outputdir3,os.path.basename(i[0])+"_cuta_2.fq.gz")
			trim1=os.path.join(outputdir4,os.path.basename(i[0])+"_trim_1.fq.gz")
			trim1_1=os.path.join(outputdir4,os.path.basename(i[0])+"_trim_1.fq")
			trim1_u=os.path.join(outputdir4,os.path.basename(i[0])+"_unpair_1.fq.gz")
			trim2=os.path.join(outputdir4,os.path.basename(i[0])+"_trim_2.fq.gz")
			trim2_u=os.path.join(outputdir4,os.path.basename(i[0])+"_unpair_2.fq.gz")
	script='''#!/bin/bash
#SBATCH -J %(sample_name)s
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH -p normal
#SBATCH --mem 220G
#SBATCH -o %(sample_name)s_2.out
#SBATCH -e %(sample_name)s_2.err
#SBATCH -t 10-00:00:00


cd %(outputdir1)s
fastqc -o %(outputdir1)s %(fq1)s %(fq2)s

cd %(outputdir3)s
cutadapt -a AGATCGGAAGAG -o %(cutadapt1)s %(fq1)s
cutadapt -a "A{10}" -u 3 -e 0.1 -m 18 %(cutadapt1)s > %(cutadapt1_1)s
#cutadapt -a "A{10}" -u 3 -e 0.1 -m 18 -M 60 %(cutadapt2)s|perl -e '$/="@ST-E"; chomp; while(<>){@array=split /\\n/; $array[0]="\@ST-E"\.$array[0]; $array[1]=~s/^.{3}//;  $array[3]=~s/^.{3}//;  if(length($array[1])){print $array[0],"\\n",$array[1],"\\n",$array[2],"\\n",$array[3],"\\n";} }' |gzip - > %(cutadapt1_1)s

/public/home/zdyfy03/software/icSHAPE-master/scripts/readCollapse.pl -U %(cutadapt1_1)s -o %(cutadapt_dup)s

cd %(outputdir4)s
trimmomatic SE -phred33 %(cutadapt1_1)s %(trim1)s LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:18
gunzip %(trim1)s
#
cd %(outputdir2)s
fastqc -o %(outputdir2)s %(trim1)s 

cd %(outputdir5)s
/public/home/zdyfy03/software/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode alignReads --genomeDir /public/home/zdyfy03/reference/human/genecode_rel38_hg38/Index/Star2.5.3a/genome_index_hang100 --readFilesIn %(trim1_1)s --outFileNamePrefix %(sample_name)s_ --outReadsUnmapped Fastx --outFilterMultimapNmax 10 --outSAMattributes All --alignIntronMin 1 --outSAMmultNmax 2 --chimOutType WithinBAM SoftClip --outSAMtype BAM Unsorted --scoreGapNoncan -4 --scoreGapATAC -4 --chimSegmentMin 15 --limitOutSJcollapsed 9000000 --limitIObufferSize 950000000 --chimJunctionOverhangMin 15 --runThreadN 32
samtools view -h %(sample_name)s_Aligned.out.bam  > %(sample_name)s_starAligned.out.sam
#
cd %(outputdir6)s
/public/home/zdyfy03/software/icSHAPE-pipe-master/PsBL/Bin_Src/sam_group_trim -in %(outputdir5)s/%(sample_name)s_starAligned.out.sam,%(outputdir5)s/%(sample_name)s_Chimeric.out.sam -out %(sample_name)s_Align_chimeric_gap.sam
samtools view -bS %(sample_name)s_Align_chimeric_gap.sam > %(sample_name)s_Align_chimeric_gap.bam
samtools sort %(sample_name)s_Align_chimeric_gap.bam > %(sample_name)s_Align_chimeric_gap.sort.bam
/public/home/zdyfy03/software/cufflinks-2.2.1.Linux_x86_64/cufflinks -p 32 -G /public/home/zdyfy03/reference/human/genecode_rel38_hg38/gencode.v38.primary_assembly.annotation.gtf --library-type fr-firststrand -o ./ %(sample_name)s_Align_chimeric_gap.sort.bam

cd %(outputdir7)s
/public/home/zdyfy03/software/icSHAPE-pipe-master/PsBL/Bin_Src/sam2dg -in %(outputdir6)s/%(sample_name)s_Align_chimeric_gap.sam -out 1_%(sample_name)s_Align_chimeric.dg -s left
/public/home/zdyfy03/software/icSHAPE-pipe-master/PsBL/Bin_Src/dg_cluster -in 1_%(sample_name)s_Align_chimeric.dg -out 2_%(sample_name)s_duplex -tag_sam %(outputdir6)s/%(sample_name)s_Align_chimeric_gap.sam,2_%(sample_name)s_Align_chimeric.tag.sam -min_support 2
mkdir 3_duplex_anno

cd 3_duplex_anno
python2 /public/home/zdyfy03/software/paris-master/PARIS-master-lipan/scripts/DG2Frame.py -i ../2_%(sample_name)s_duplex -g /public/home/zdyfy03/reference/human/genecode_rel38_hg38/gencode.v38.primary_assembly.annotation.gtf -o 3_%(sample_name)s.duplex -s gencode --genomeCoor /public/home/zdyfy03/reference/human/genecode_rel38_hg38/Paris/gencode.v38.primary.genomeCoor.bed
grep -v "NULL" 3_%(sample_name)s.duplex.gframe|awk '$9 >=2'|sed 's/:/\t/g'|sed 's/||/\t/g'|cut -f1-13,16,21|awk 'FNR==NR{a[$14]=$0;next}($6 in a){print a[$6]"\t"$0}' - /public/home/zdyfy03/reference/human/genecode_rel38_hg38/Paris/gencode.v38.primary.genomeCoor_add_genename.bed|cut -f1-13,15,16-20,21,23|sort|uniq|awk 'FNR==NR{a[$6]=$0;next}($14 in a){print $0"\t"a[$14]}' /public/home/zdyfy03/reference/human/genecode_rel38_hg38/Paris/gencode.v38.primary.genomeCoor_add_genename.bed -|cut -f1-13,19-21,26-27,29|sort|uniq > 4_%(sample_name)s.duplex.gframe.add_anno

echo "all duplex" >> readcount
echo "cut -f1 4_%(sample_name)s.duplex.gframe.add_anno|sort|uniq|wc -l" >> readcount
echo "intre group" >> readcount
echo "awk '$15!=$18' 4_%(sample_name)s.duplex.gframe.add_anno|cut -f1|sort|uniq|wc -l" >> readcount
echo "intra group" >> readcount
echo "awk '$15==$18' 4_%(sample_name)s.duplex.gframe.add_anno|cut -f1|sort|uniq|wc -l" >> readcount



	''' %(locals())
	open(sample_name+"_paris.sh","w").write(script)
	#os.system("qsub "+sample_name+"_RNA.sh")
