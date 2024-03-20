# Genome indexing ([STAR](https://github.com/alexdobin/STAR))

Fasta file: [hg19.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)

GTF file: [hg19.refGene.gtf](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz)

```
./STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir /home/ec2-user/UCSC_hg19_index \
--genomeFastaFiles /home/ec2-user/UCSC/hg19.fa \
--sjdbGTFfile /home/ec2-user/UCSC/hg19.refGene.gtf
````

# Genome mapping ([STAR](https://github.com/alexdobin/STAR))

```
$ for f in `ls *.fastq.gz | sed 's/_[12].fastq.gz//g' | sort -u`; do echo $f && \
STAR --runThreadN 16 \
--genomeDir /home/ec2-user/UCSC_hg19_index \
--readFilesCommand zcat \
--readFilesIn ${f}_1.fastq.gz ${f}_2.fastq.gz \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /home/ec2-user/STAR/${f}_; done
```

# Gene counting ([featureCounts](https://subread.sourceforge.net/featureCounts.html))

```
./featureCounts -T 16 \
-M --fraction -p \
-a /home/ec2-user/UCSC/hg19.refGene.gtf \
-o /home/ec2-user/featureCounts \
/home/ec2-user/STAR/*.bam
```

# Transcript counting ([cuffdiff](https://cole-trapnell-lab.github.io/cufflinks/cuffdiff/))

```
./cuffdiff /home/ec2-user/UCSC/hg19.refGene.gtf \
[group1 bam list] \
[group2 bam list] \
-o /home/ec2-user/cuffdiff/ \
-p 16 -L group1,group2 --library-type fr-firststrand
```

# Splicing analysis ([rMATS_turbo](https://github.com/Xinglab/rmats-turbo))

```
python3 rmats.py \
--b1 /home/ec2-user/STAR/GW.txt \
--b2 /home/ec2-user/STAR/WT.txt \
--gtf /home/ec2-user/UCSC/hg19.refGene.gtf \
-t paired --nthread 16 \
--od /home/ec2-user/rMATS \
--tmp /home/ec2-user/rMATS/tmp \
--cstat 0.0001 --readLength 101
```

# Peak calling ([MACS2](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html))

```
./macs2 callpeak \
-t /home/ec2-user/RIP/sample.bam \
-c /home/ec2-user/RIP/input_control.bam \
-f BAM -g hs \
--outdir /home/ec2-user/RIP/MACS2 \
-n sample -B -q 0.0001
```
