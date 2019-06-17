
## Table of contents

- [Software requirements](#software-requirements)
- [Libraries](#libraries)
- [Quality check](#quality-check)
- [Trim adapters and filter bases based on quality](#trim-adapters-and-filter-bases-based-on-quality)
- [Discard reads containing protecting group sequence](#discard-reads-containing-protecting-group-sequence)
- [Alignment](#alignment)
- [Spike-ins read count analysis](#spike-ins-read-count-analysis)
- [Filtering genomic libraries](#filtering-genomic-libraries)
- [Create tdf files](#create-tdf-files)
- [Insert size plots](#insert-size-plots)
- [Peak calling](#peak-calling)
- [Consensus peaks](#consensus-peaks)
- [Genomic Association Tester (GAT) analysis](#genomic-association-tester-gat-analysis)
- [Base composition tables](#base-composition-tables)



## Software requirements

- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [sambamba v0.6.5](https://academic.oup.com/bioinformatics/article/31/12/2032/214758)
- Standard Unix tools: awk, sort, uniq
- [igvtools v2.3.91](https://software.broadinstitute.org/software/igv/igvtools)
- [bedtools v2.27.0](http://bedtools.readthedocs.io/en/latest/)
- [tableCat.py](https://github.com/dariober/bioinformatics-cafe/blob/master/tableCat/tableCat.py)
- [macs2 v2.1.1.20160309](https://github.com/taoliu/MACS)
- [fastaRegexFinder.py v0.1.1](https://github.com/dariober/bioinformatics-cafe/tree/master/fastaRegexFinder)
- [gat-run.py](http://gat.readthedocs.io/en/latest/contents.html)
- [python v2.7.12](https://www.python.org/). Libraries:
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)
  - [GenomicFeatures v1.26.4](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
  - [ggplot2 v2.2.1](http://ggplot2.org/)
  - [edgeR v3.16.5](https://bioconductor.org/packages/release/bioc/html/edgeR.html)



## Libraries

ArrayExpress E-MTAB-7152 naming:

Library | Biological replicate | Sequencing type| script id
:------:|:--------------------:|:--------------:|:--------:
APE1-siRNA_snAP1 | rep1 | paired-end | Lib134
APE1-siRNA_snAP2 | rep2 | paired-end | Lib135
APE1-siRNA_snAP3 | rep3 | paired-end | Lib145
APE1-siRNA_snAP4 | rep4 | paired-end | Lib146
APE1-siRNA_input1 | rep1 | paired-end | Lib136
APE1-siRNA_input2 | rep2 | paired-end | Lib137
APE1-siRNA_input3 | rep3 | paired-end | Lib147
APE1-siRNA_input4 | rep4 | paired-end | Lib148
Cont-siRNA_snAP1 | rep1 | paired-end | Lib141
Cont-siRNA_snAP2 | rep2 | paired-end | Lib142
Cont-siRNA_snAP3 | rep3 | paired-end | Lib149
Cont-siRNA_snAP4 | rep4 | paired-end | Lib150
Cont-siRNA_input1 | rep1 | paired-end | Lib143
Cont-siRNA_input2 | rep2 | paired-end | Lib144
Cont-siRNA_input3 | rep3 | paired-end | Lib151
Cont-siRNA_input4 | rep4 | paired-end | Lib152

Bear the mapping above in mind for the code following from now:



## Renaming

`*_1.fastq.gz` -> `*_R1.fastq.gz` and `*_2.fastq.gz` -> `*_R2.fastq.gz`



## Quality check

```bash
cd ~/fastq # directory containing the fastq sequencing files
mkdir ../fastqc
for fq in *.fastq.gz
do
  fastqc --noextract -q -o ../fastqc $fq
done
```



## Trim adapters and filter bases based on quality

```bash
cd ~/fastq
mkdir ../fastq_trimmed

for fq1 in *R1*.fastq.gz
do
  fq2=${fq1/_R1/_R2}
  bname=${fq1%.fastq.gz}
  cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 15 -q 20 -o ../fastq_trimmed/$fq1 -p ../fastq_trimmed/$fq2 $fq1 $fq2 > ../fastq_trimmed/${bname}.txt
done
```



## Discard reads containing protecting group sequence

```bash
cd ~/fastq_trimmed
mkdir ../fastq_trimmed_protecting_sequence_optimal

for fq1 in *R1*.fastq.gz
do
  fq2=${fq1/_R1/_R2}
  bname=${fq1%_R1.fastq.gz}
  cutadapt -g GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -G GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -e 0.1 -O 20 --discard-trimmed --pair-filter=any -o ../fastq_trimmed_protecting_sequence_optimal/$fq1 -p ../fastq_trimmed_protecting_sequence_optimal/$fq2 $fq1 $fq2 > ../fastq_trimmed_protecting_sequence_optimal/$bname.txt
done
```



## Alignment

### Prepare and index references

```bash
cd ~/reference/

wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
tar -xvf Homo_sapiens_UCSC_hg38.tar.gz

cat genome.fa spikeins.fa > genome_spikeins.fa # see oligo.md for the contents of spikeins.fa
bwa index genome_spikeins.fa
samtools faidx genome_spikeins.fa
```


### Align and sort

```bash
cd ~/fastq_trimmed_protecting_sequence_optimal

mkdir ../bam

ref=../reference/genome_spikeins.fa

for fq1 in *R1*.fastq.gz
do
  bname=${fq1%_R1.fastq.gz}
  fq2=${fq1/_R1_/_R2_}
  bwa mem -t 20 -M $ref $fq1 $fq2 | \
  samtools view -@ 20 -b - | \
  samtools sort -@ 20 -T ~/tmp/$bname -o ../bam/$bname.tmp.bam -
done
```


### Merge lanes, mark duplicates, index and flagstat

```bash
cd ~/bam

mkdir ../flagstat

listOfIds="..." # add here the list of ids for the lanes that you would like to merge. If there is no need, then omit the samtools merge step in the loop below

for id in listOfIds
do
  bams=`echo ${id}_L00[1-4].tmp.bam`
  samtools merge -@ 20 -f ${id}.tmp.bam $bams && \
  sambamba markdup -t 20 ${id}.tmp.bam ${id}.bam 2> ${id}.markdup.txt && \
  samtools flagstat ${id}.bam > ../flagstat/${id}.txt
done
```



## Spike-ins read count analysis

### Overall

```bash
cd ~/bam

for bam in `ls *.bam | grep -v "tmp"`
do
bname=${bam%.bam}
AP1_fwd=`samtools view -c -F 16 $bam AP1`
AP1_rev=`samtools view -c -f 16 $bam AP1`
GCAT1_fwd=`samtools view -c -F 16 $bam GCAT1`
GCAT1_rev=`samtools view -c -f 16 $bam GCAT1`
fU1_fwd=`samtools view -c -F 16 $bam fU1`
fU1_rev=`samtools view -c -f 16 $bam fU1`
fC1_fwd=`samtools view -c -F 16 $bam fC1`
fC1_rev=`samtools view -c -f 16 $bam fC1`
AP2_fwd=`samtools view -c -F 16 $bam AP2`
AP2_rev=`samtools view -c -f 16 $bam AP2`
GCAT2_fwd=`samtools view -c -F 16 $bam GCAT2`
GCAT2_rev=`samtools view -c -f 16 $bam GCAT2`
fU2_fwd=`samtools view -c -F 16 $bam fU2`
fU2_rev=`samtools view -c -f 16 $bam fU2`
fC2_fwd=`samtools view -c -F 16 $bam fC2`
fC2_rev=`samtools view -c -f 16 $bam fC2`
echo -e "${bname}\t${AP1_fwd}\t${AP1_rev}\t${GCAT1_fwd}\t${GCAT1_rev}\t${fU1_fwd}\t${fU1_rev}\t${fC1_fwd}\t${fC1_rev}\t${AP2_fwd}\t${AP2_rev}\t${GCAT2_fwd}\t${GCAT2_rev}\t${fU2_fwd}\t${fU2_rev}\t${fC2_fwd}\t${fC2_rev}"
done | column -t
```


### All libraries

```bash
cd ~/bam

mkdir ../reference_distribution

for bam in `ls *.bam | grep -v "tmp"`
do
  samtools view $bam -F260 -f64 -f32 | cut -f3 | sort | uniq -c | grep -E 'AP1|AP2|fC1|fC2|fU1|fU2|GCAT1|GCAT2' > ../reference_distribution/${bam%.bam}_R1_fwd.txt && \
  samtools view $bam -F260 -f64 -f16 | cut -f3 | sort | uniq -c | grep -E 'AP1|AP2|fC1|fC2|fU1|fU2|GCAT1|GCAT2' > ../reference_distribution/${bam%.bam}_R1_rev.txt && \
  samtools view $bam -F260 -f128 -f32 | cut -f3 | sort | uniq -c | grep -E 'AP1|AP2|fC1|fC2|fU1|fU2|GCAT1|GCAT2' > ../reference_distribution/${bam%.bam}_R2_fwd.txt && \
  samtools view $bam -F260 -f128 -f16 | cut -f3 | sort | uniq -c | grep -E 'AP1|AP2|fC1|fC2|fU1|fU2|GCAT1|GCAT2' > ../reference_distribution/${bam%.bam}_R2_rev.txt
done
```

Generating a fragmentation pattern file for each library, R1/R2 and fwd/rev.

```bash
cd ~/bam

mkdir ../fragmentation_analysis

for bam in `ls *.bam | grep -v "tmp"`
do
  for ref in AP1 AP2 fC1 fC2 fU1 fU2 GCAT1 GCAT2
  do
      nohup samtools view $bam -F260 -f64 -f32 | awk -v ref=$ref '$3==ref {print $4}' | sort -k1,1n | uniq -c | awk '{t = $1; $1 = $2; $2 = t; print; }' > ../fragmentation_analysis/${bam%.bam}_R1_fwd_${ref}.txt &
      nohup samtools view $bam -F260 -f64 -f16 | awk -v ref=$ref '$3==ref {print $4}' | sort -k1,1n | uniq -c | awk '{ t = $1; $1 = $2; $2 = t; print; }' > ../fragmentation_analysis/${bam%.bam}_R1_rev_${ref}.txt &
      nohup samtools view $bam -F260 -f128 -f32 | awk -v ref=$ref '$3==ref {print $4}' | sort -k1,1n | uniq -c | awk '{ t = $1; $1 = $2; $2 = t; print; }' > ../fragmentation_analysis/${bam%.bam}_R2_fwd_${ref}.txt &
      nohup samtools view $bam -F260 -f128 -f16 | awk -v ref=$ref '$3==ref {print $4}' | sort -k1,1n | uniq -c | awk '{ t = $1; $1 = $2; $2 = t; print; }' > ../fragmentation_analysis/${bam%.bam}_R2_rev_${ref}.txt &
  done
done
```



## Filtering genomic libraries

Generate a whitelist from blacklisted regions:

```bash
cd ~/reference
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -A -e "select chrom, size from hg38.chromInfo" | awk -v OFS="\t" '{print $1, 0, $2}' | bedtools subtract -a - -b hg38.blacklist.bed.gz | sort -k1,1 -k2,2n > hg38.whitelist.bed
```

Filter:

```bash
cd ~/bam

whitelist=~/reference/hg38.whitelist.bed

for bam in `*.bam | grep -v "tmp"`
do
  bname=${bam%.bam}
  samtools idxstats $bam | cut -f1 | grep 'chr' | xargs samtools view -@ 20 -F 3844 -q 10 -L $whitelist -b $bam > $bname.clean.bam && samtools index $bname.clean.bam && \
  samtools view -@ 20 -f 64 -b $bname.clean.bam > $bname.clean.R1.bam && samtools index $bname.clean.R1.bam && \
  samtools view -@ 20 -f 64 -f 32 -b $bname.clean.bam > $bname.clean.R1.fwd.bam && samtools index $bname.clean.R1.fwd.bam && \
  samtools view -@ 20 -f 64 -f 16 -b $bname.clean.bam > $bname.clean.R1.rev.bam && samtools index $bname.clean.R1.rev.bam
done
```


### Filter R1.fwd and R1.rev by CIGAR

```bash
cd ~/bam

#R1.fwd
for bam in *.clean.R1.fwd.bam
do
  bname=${bam%.bam}
  samtools view -@ 20 -h $bam | awk '$6 !~ /^[0-9]+S/' | samtools view -@ 20 -b - > ${bname}.cigar.bam && samtools index ${bname}.cigar.bam
done

#R1.rev
for bam in *.clean.R1.rev.bam
do
  bname=${bam%.bam}
  nohup samtools view -@ 20 -h $bam | awk '$6 !~ /S$/' | samtools view -@ 20 -b - > ${bname}.cigar.bam && samtools index ${bname}.cigar.bam &
done
```



## Create tdf files

```bash
cd ~/bam

mkdir ../tdf

#paired-end
for bam in *.clean.bam
do
  bname=${bam%.bam}
  igvtools count $bam ../tdf/$bname.tdf ../reference/genome_spikeins.fa.fai --includeDuplicates --pairs -w 1 -e 0 --preExtFactor 0 --postExtFactor 0
done

#single-end
for bam in *.clean.R*.bam
do
  bname=${bam%.bam}
  igvtools count $bam ../tdf/$bname.tdf ../reference/genome_spikeins.fa.fai --includeDuplicates -w 1 -e 0 --preExtFactor 0 --postExtFactor 0
done
```



## Insert size plots

```bash
cd ~/bam

mkdir ../stats

for bam in *.clean.bam
do
  bname=${bam%.bam}
  mkdir ../stats/$bname
  samtools stats $bam > ../stats/$bname/$bname.txt && \
  plot-bamstats -p ../stats/$bname/ ../stats/$bname/$bname.txt
done
```



## Peak calling

### macs2

```bash
cd ~/bam

mkdir -p ../macs2/default

##############
# APE1-siRNA #
##############

# Lib134_HeLa_siRNA3_Z_AP1 vs. Lib136_HeLa_siRNA3_Z_Y1
bam_t=Lib134_HeLa_siRNA3_Z_AP1_S5.hg38.clean.bam
bam_c=Lib136_HeLa_siRNA3_Z_Y1_S4.hg38.clean.bam

bname=`basename ${bam_t%.bam}`

macs2 callpeak \
-t $bam_t \
-c $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel.p1e-5 \
--tempdir ~/tmp/ \
--nomodel \
-p 0.00001 \
--cutoff-analysis

bname=`basename ${bam_c%.bam}`

macs2 callpeak \
-t $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel \
--tempdir ~/tmp/ \
--nomodel \
--cutoff-analysis"

bed_a=../macs2_hg38/default/Lib136_HeLa_siRNA3_Z_Y1_S4.hg38.clean.nomodel_peaks.narrowPeak
bed_t=../macs2_hg38/default/Lib134_HeLa_siRNA3_Z_AP1_S5.hg38.clean.nomodel.p1e-5_peaks.narrowPeak
bname=${bed_t%_peaks.narrowPeak}
bedtools intersect -a $bed_t -b $bed_a -v > ${bname}.subtracted_peaks.narrowPeak



# Lib135_HeLa_siRNA3_Z_AP2 vs. Lib137_HeLa_siRNA3_Z_Y2
bam_t=Lib135_HeLa_siRNA3_Z_AP2_S1.hg38.clean.bam
bam_c=Lib137_HeLa_siRNA3_Z_Y2_S2.hg38.clean.bam

bname=`basename ${bam_t%.bam}`

macs2 callpeak \
-t $bam_t \
-c $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel.p1e-5 \
--tempdir ~/tmp/ \
--nomodel \
-p 0.00001 \
--cutoff-analysis

bname=`basename ${bam_c%.bam}`

macs2 callpeak \
-t $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel \
--tempdir ~/tmp/ \
--nomodel \
--cutoff-analysis"

bed_a=../macs2_hg38/default/Lib137_HeLa_siRNA3_Z_Y2_S2.hg38.clean.nomodel_peaks.narrowPeak
bed_t=../macs2_hg38/default/Lib135_HeLa_siRNA3_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak
bname=${bed_t%_peaks.narrowPeak}
bedtools intersect -a $bed_t -b $bed_a -v > ${bname}.subtracted_peaks.narrowPeak



# Lib145_HeLa_siRNA3_Z_AP3 vs. Lib147_HeLa_siRNA_Z_Y3
bam_t=Lib145_HeLa_siRNA3_Z_AP3_S1.hg38.clean.bam
bam_c=Lib147_HeLa_siRNA_Z_Y3_S2.hg38.clean.bam

bname=`basename ${bam_t%.bam}`

macs2 callpeak \
-t $bam_t \
-c $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel.p1e-5 \
--tempdir ~/tmp/ \
--nomodel \
-p 0.00001 \
--cutoff-analysis

bname=`basename ${bam_c%.bam}`

macs2 callpeak \
-t $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel \
--tempdir ~/tmp/ \
--nomodel \
--cutoff-analysis"

bed_a=../macs2_hg38/default/Lib147_HeLa_siRNA_Z_Y3_S2.hg38.clean.nomodel_peaks.narrowPeak
bed_t=../macs2_hg38/default/Lib145_HeLa_siRNA3_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak
bname=${bed_t%_peaks.narrowPeak}
bedtools intersect -a $bed_t -b $bed_a -v > ${bname}.subtracted_peaks.narrowPeak



# Lib146_HeLa_siRNA_Z_AP4 vs. Lib148_HeLa_siRNA_Z_Y4
bam_t=Lib146_HeLa_siRNA_Z_AP4_S1.hg38.clean.bam
bam_c=Lib148_HeLa_siRNA_Z_Y4_S2.hg38.clean.bam

bname=`basename ${bam_t%.bam}`

macs2 callpeak \
-t $bam_t \
-c $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel.p1e-5 \
--tempdir ~/tmp/ \
--nomodel \
-p 0.00001 \
--cutoff-analysis

bname=`basename ${bam_c%.bam}`

macs2 callpeak \
-t $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel \
--tempdir ~/tmp/ \
--nomodel \
--cutoff-analysis"

bed_a=../macs2_hg38/default/Lib148_HeLa_siRNA_Z_Y4_S2.hg38.clean.nomodel_peaks.narrowPeak
bed_t=../macs2_hg38/default/Lib146_HeLa_siRNA_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak
bname=${bed_t%_peaks.narrowPeak}
bedtools intersect -a $bed_t -b $bed_a -v > ${bname}.subtracted_peaks.narrowPeak



##############
# Cont-siRNA #
##############

# Lib141_HeLa_cp_Z_AP1 vs. Lib143_HeLa_cp_Z_Y1
bam_t=Lib141_HeLa_cp_Z_AP1_S1.hg38.clean.bam
bam_c=Lib143_HeLa_cp_Z_Y1_S2.hg38.clean.bam

bname=`basename ${bam_t%.bam}`

macs2 callpeak \
-t $bam_t \
-c $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel.p1e-5 \
--tempdir ~/tmp/ \
--nomodel \
-p 0.00001 \
--cutoff-analysis

bname=`basename ${bam_c%.bam}`

macs2 callpeak \
-t $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel \
--tempdir ~/tmp/ \
--nomodel \
--cutoff-analysis"

bed_a=../macs2_hg38/default/Lib143_HeLa_cp_Z_Y1_S2.hg38.clean.nomodel_peaks.narrowPeak
bed_t=../macs2_hg38/default/Lib141_HeLa_cp_Z_AP1_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak
bname=${bed_t%_peaks.narrowPeak}
bedtools intersect -a $bed_t -b $bed_a -v > ${bname}.subtracted_peaks.narrowPeak



# Lib142_HeLa_cp_Z_AP2 vs. Lib144_HeLa_cp_Z_Y2
bam_t=Lib142_HeLa_cp_Z_AP2_S1.hg38.clean.bam
bam_c=Lib144_HeLa_cp_Z_Y2_S2.hg38.clean.bam

bname=`basename ${bam_t%.bam}`

macs2 callpeak \
-t $bam_t \
-c $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel.p1e-5 \
--tempdir ~/tmp/ \
--nomodel \
-p 0.00001 \
--cutoff-analysis

bname=`basename ${bam_c%.bam}`

macs2 callpeak \
-t $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel \
--tempdir ~/tmp/ \
--nomodel \
--cutoff-analysis"

bed_a=../macs2_hg38/default/Lib144_HeLa_cp_Z_Y2_S2.hg38.clean.nomodel_peaks.narrowPeak
bed_t=../macs2_hg38/default/Lib142_HeLa_cp_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak
bname=${bed_t%_peaks.narrowPeak}
bedtools intersect -a $bed_t -b $bed_a -v > ${bname}.subtracted_peaks.narrowPeak



# Lib149_HeLa_cp_Z_AP3 vs. Lib151_HeLa_cp_Z_Y3
bam_t=Lib149_HeLa_cp_Z_AP3_S1.hg38.clean.bam
bam_c=Lib151_HeLa_cp_Z_Y3_S2.hg38.clean.bam

bname=`basename ${bam_t%.bam}`

macs2 callpeak \
-t $bam_t \
-c $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel.p1e-5 \
--tempdir ~/tmp/ \
--nomodel \
-p 0.00001 \
--cutoff-analysis

bname=`basename ${bam_c%.bam}`

macs2 callpeak \
-t $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel \
--tempdir ~/tmp/ \
--nomodel \
--cutoff-analysis"

bed_a=../macs2_hg38/default/Lib151_HeLa_cp_Z_Y3_S2.hg38.clean.nomodel_peaks.narrowPeak
bed_t=../macs2_hg38/default/Lib149_HeLa_cp_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak
bname=${bed_t%_peaks.narrowPeak}
bedtools intersect -a $bed_t -b $bed_a -v > ${bname}.subtracted_peaks.narrowPeak



# Lib150_HeLa_cp_Z_AP4 vs. Lib152_HeLa_cp_Z_Y4
bam_t=Lib150_HeLa_cp_Z_AP4_S1.hg38.clean.bam
bam_c=Lib152_HeLa_cp_Z_Y4_S2.hg38.clean.bam

bname=`basename ${bam_t%.bam}`

macs2 callpeak \
-t $bam_t \
-c $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel.p1e-5 \
--tempdir ~/tmp/ \
--nomodel \
-p 0.00001 \
--cutoff-analysis

bname=`basename ${bam_c%.bam}`

macs2 callpeak \
-t $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2_hg38/default/ \
-n $bname.nomodel \
--tempdir ~/tmp/ \
--nomodel \
--cutoff-analysis"

bed_a=../macs2_hg38/default/Lib152_HeLa_cp_Z_Y4_S2.hg38.clean.nomodel_peaks.narrowPeak
bed_t=../macs2_hg38/default/Lib150_HeLa_cp_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak
bname=${bed_t%_peaks.narrowPeak}
bedtools intersect -a $bed_t -b $bed_a -v > ${bname}.subtracted_peaks.narrowPeak
```



## Consensus peaks

```bash
cd ~/macs2/default

##############
# APE1-siRNA #
##############

tableCat.py -i \
180326_NS500222_0389_HMV3JBGX5_Lib134_HeLa_siRNA3_Z_AP1_S5.hg38.clean.nomodel.p1e-5.subtracted_peaks.narrowPeak \
180427_NS500222_0395_HTCGYBGX5_Lib135_HeLa_siRNA3_Z_AP2_S1.hg38.clean.nomodel.p1e-5.subtracted_peaks.narrowPeak \
180503_NS500222_0398_HTCJ3BGX5_Lib145_HeLa_siRNA3_Z_AP3_S1.hg38.clean.nomodel.p1e-5.subtracted_peaks.narrowPeak \
180504_NS500222_0399_HTGKFBGX5_Lib146_HeLa_siRNA_Z_AP4_S1.hg38.clean.nomodel.p1e-5.subtracted_peaks.narrowPeak | \
awk -v OFS="\t" '{print $1, $2, $3, $NF}' | \
sort -k1,1 -k2,2n | \
bedtools merge -c 4,4 -o distinct,count_distinct -i - | \
awk -v OFS="\t" '$5 > 2 {print $1, $2, $3}' > Lib134_Lib135_Lib145_Lib146_merge_p1e-5_subtracted.bed

wc -l Lib134_Lib135_Lib145_Lib146_merge_p1e-5_subtracted.bed # 25080


##############
# Cont-siRNA #
##############

tableCat.py -i \
180501_NS500222_0396_HTGMFBGX5_Lib141_HeLa_cp_Z_AP1_S1.hg38.clean.nomodel.p1e-5.subtracted_peaks.narrowPeak \
180502_NS500222_0397_HTGM5BGX5_Lib142_HeLa_cp_Z_AP2_S1.hg38.clean.nomodel.p1e-5.subtracted_peaks.narrowPeak \
180508_NS500222_0401_HMVKKBGX5_Lib149_HeLa_cp_Z_AP3_S1.hg38.clean.nomodel.p1e-5.subtracted_peaks.narrowPeak \
180509_NS500222_0402_HMJJHBGX5_Lib150_HeLa_cp_Z_AP4_S1.hg38.clean.nomodel.p1e-5.subtracted_peaks.narrowPeak | \
awk -v OFS="\t" '{print $1, $2, $3, $NF}' | \
sort -k1,1 -k2,2n | \
bedtools merge -c 4,4 -o distinct,count_distinct -i - | \
awk -v OFS="\t" '$5 > 2 {print $1, $2, $3}' > Lib141_Lib142_Lib149_Lib150_merge_p1e-5_subtracted.bed

wc -l Lib141_Lib142_Lib149_Lib150_merge_p1e-5_subtracted.bed # 14110
```

Intersections and union of `Lib134_Lib135_Lib145_Lib146_merge_p1e-5_subtracted.bed` and `Lib141_Lib142_Lib149_Lib150_merge_p1e-5_subtracted.bed`:

```bash
cd ~/macs2_hg38/default

wc -l Lib134_Lib135_Lib145_Lib146_merge_p1e-5_subtracted.bed Lib141_Lib142_Lib149_Lib150_merge_p1e-5_subtracted.bed
# 25080 Lib134_Lib135_Lib145_Lib146_merge_p1e-5_subtracted.bed
# 14110 Lib141_Lib142_Lib149_Lib150_merge_p1e-5_subtracted.bed

bedtools intersect \
-a Lib134_Lib135_Lib145_Lib146_merge_p1e-5_subtracted.bed \
-b Lib141_Lib142_Lib149_Lib150_merge_p1e-5_subtracted.bed \
-sorted \
-wa -u | wc -l # 10261 Lib134_Lib135_Lib145_Lib146_merge_p1e-5_subtracted regions intersect

bedtools intersect \
-a Lib141_Lib142_Lib149_Lib150_merge_p1e-5_subtracted.bed \
-b Lib134_Lib135_Lib145_Lib146_merge_p1e-5_subtracted.bed \
-sorted \
-wa -u | wc -l # 10387 (100*10387/14110 = 74%) Lib141_Lib142_Lib149_Lib150_merge_p1e-5_subtracted regions intersect
```



## Genomic Association Tester (GAT) analysis

### gene features

Using [GenomicFeatures](http://bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf):

```r
library(data.table)
library(GenomicFeatures)

# change width
options(width = 250)

# prepare coordinates table
txdb <- makeTxDbFromGFF("~/reference/genes.gtf", format="gtf")
seqlevels(txdb)
columns(txdb)
keytypes(txdb)

# promoters
promoter <- data.table(data.frame(promoters(genes(txdb), upstream=1000, downstream=0)))[!grepl("_", seqnames), c("seqnames", "start", "end", "gene_id", "strand")][order(seqnames, start)]
nrow(promoter) # 25043
promoter[, featuretype := "promoter_1kb"]
write.table(promoter, file = '~/annotation/hg38_promoter_gene.bed', row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

# 5'UTR
utr5 <- data.table(data.frame(fiveUTRsByTranscript(txdb, use.names = TRUE)))[!grepl("_", seqnames), c("seqnames", "start", "end", "group_name", "strand")][order(seqnames, start)]
nrow(utr5) # 65487
utr5[, featuretype := "5'UTR"]
write.table(utr5, file = '~/annotation/hg38_utr5_transcript.bed', row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

# exons
exons_gene <- data.table(data.frame(exonsBy(txdb, by = "gene")))[!grepl("_", seqnames), c("seqnames", "start", "end", "group_name", "strand")][order(seqnames, start)]
nrow(exons_gene) # 250153
exons_gene[, featuretype := "exon"]
write.table(exons_gene, file = '~/annotation/hg38_exon_gene.bed', row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

exons_transcript <- data.table(data.frame(exonsBy(txdb, by = "tx", use.names = TRUE)))[!grepl("_", seqnames), c("seqnames", "start", "end", "group_name", "strand")][order(seqnames, start)]
nrow(exons_transcript) # 484310
exons_transcript[, featuretype := "exon"]
write.table(exons_transcript, file = '~/annotation/hg38_exon_transcript.bed', row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

# introns
introns <- data.table(data.frame(intronsByTranscript(txdb, use.names = TRUE)))[!grepl("_", seqnames), c("seqnames", "start", "end", "group_name", "strand")][order(seqnames, start)]
nrow(introns) # 432392
introns[, featuretype := "intron"]
write.table(introns, file = '~/annotation/hg38_intron_transcript.bed', row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

# 3'UTR
utr3 <- data.table(data.frame(threeUTRsByTranscript(txdb, use.names = TRUE)))[!grepl("_", seqnames), c("seqnames", "start", "end", "group_name", "strand")][order(seqnames, start)]
nrow(utr3) # 41533
utr3[, featuretype := "3'UTR"]
write.table(utr3, file = '~/annotation/hg38_utr3_transcript.bed', row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

# genes
gene <- data.table(data.frame(genes(txdb)))[!grepl("_", seqnames), c("seqnames", "start", "end", "gene_id", "strand")][order(seqnames, start)]
nrow(gene) # 25043
gene[, featuretype := "gene"]
write.table(gene, file = '~/annotation/hg38_gene.bed', row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
```


### intergenic

```bash
cd ~/annotation/

genome=~/reference/genome.fa.fai

complementBed -g <(sort -k1,1 -k2,2n $genome) -i <(cat hg38_gene.bed hg38_exon_transcript.bed hg38_intron_transcript.bed | bedtools sort -i | bedtools merge) | \
grep -v -E "_|chrEBV" | \
sed 's/$/&\tintergenic/' > hg38_intergenic.bed

wc -l hg38_intergenic.bed # 35958
```



### combine gene features and intergenic regions

```bash
cd ~/annotation

cat <(cut -f1-3,6 hg38_promoter_gene.bed hg38_utr5_transcript.bed hg38_exon_transcript.bed hg38_intron_transcript.bed hg38_utr3_transcript.bed) \
hg38_intergenic.bed | \
bedtools sort -i - > hg38_promoter_utr5_exon_intron_utr3_intergenic.bed
```


### define mappable genome

Need to define two mappable files: APE1-siRNA and Cont-siRNA:

- Lib136_HeLa_siRNA3_Z_Y1_S4.hg38.clean.bam
- Lib137_HeLa_siRNA3_Z_Y2_S2.hg38.clean.bam
- Lib147_HeLa_siRNA_Z_Y3_S2.hg38.clean.bam
- Lib148_HeLa_siRNA_Z_Y4_S2.hg38.clean.bam

- Lib143_HeLa_cp_Z_Y1_S2.hg38.clean.bam
- Lib144_HeLa_cp_Z_Y2_S2.hg38.clean.bam
- Lib151_HeLa_cp_Z_Y3_S2.hg38.clean.bam
- Lib152_HeLa_cp_Z_Y4_S2.hg38.clean.bam

```bash
cd ~/bam

g=~/reference/genome.fa.fai

# APE1-siRNA
samtools merge -@ 20 APE1.Y.tmp.bam *siRNA*Y*.clean.bam
bedtools genomecov -ibam APE1.Y.tmp.bam -bg -g $g | grep -v -E '_|chrEBV' | sort -k1,1 -k2,2n | mergeBed > APE1.Y.mappable.bed
rm APE1.Y.tmp.bam

# Cont-siRNA
samtools merge -@ 20 WT.Y.tmp.bam *cp*Y*.clean.bam &
bedtools genomecov -ibam WT.Y.tmp.bam -bg -g $g | grep -v -E '_|chrEBV' | sort -k1,1 -k2,2n | mergeBed > WT.Y.mappable.bed
rm WT.Y.tmp.bam
```


### run gat.py

```bash
cd ~/macs2_hg38/default

mkdir ~/gat

annotations=~/annotation/hg38_promoter_utr5_exon_intron_utr3_intergenic.bed

# APE1-siRNA
mappable=~/bam/APE1.Y.mappable.bed

bed=Lib134_Lib135_Lib145_Lib146_merge_p1e-5_subtracted.bed
bname=`basename ${bed%.bed}`
gat-run.py -a $annotations -s <(grep -v -E '_|chrEBV' $bed) -w $mappable --ignore-segment-tracks -n 10000 -t 20 -L ../../gat/$bname.log > ../../gat/$bname.txt


# Cont-siRNA
mappable=~/bam/WT.Y.mappable.bed

bed=Lib141_Lib142_Lib149_Lib150_merge_p1e-5_subtracted.bed
bname=`basename ${bed%.bed}`
gat-run.py -a $annotations -s <(grep -v -E '_|chrEBV' $bed) -w $mappable --ignore-segment-tracks -n 10000 -t 20 -L ../../gat/$bname.log > ../../gat/$bname.txt
```



## Base composition tables

Counts for position 5' - 1 (R1.fwd) and 5' + 1 (R1.rev)

```bash
cd ~/bam

mkdir ../base_composition

g=~/reference/genome_spikeins.fa

for bam in *.clean.R1.fwd.cigar.bam; do
  bname=${bam%.bam}
  bedtools bamtobed -i $bam | cut -f1-3 | grep -v -E '_|chrEBV' | bedtools sort -i - | bedtools flank -i - -g $g.fai -l 1 -r 0 | bedtools getfasta -fi $g -bed - -tab | cut -f2 | tr a-z A-Z | sort | uniq -c > ../base_composition/$bname.txt
done

for bam in *.clean.R1.rev.cigar.bam; do
  bname=${bam%.bam}
  bedtools bamtobed -i $bam | cut -f1-3 | grep -v -E '_|chrEBV' | bedtools sort -i - | bedtools flank -i - -g $g.fai -l 0 -r 1 | bedtools getfasta -fi $g -bed - -tab | cut -f2 | tr a-z A-Z | sort | uniq -c > ../base_composition/$bname.txt
done
```

What are A, C, G, T frequencies in hg38?

```python
ifasta=open("~/reference/genome_spikeins.fa", "r")
chromosomes=ifasta.read().split(">")
ifasta.close()

sequence = ""

for chr in chromosomes:
  lines=chr.split("\n")
  name=lines[0]
  if name.startswith("chr") and ("_" not in name) and (name != "chrEBV"):
    print name
    seq="".join(lines[1:]).upper()
    sequence=sequence+seq


sequence.count('A') # 863155826
sequence.count('C') # 596433284
sequence.count('G') # 598488389
sequence.count('T') # 865655055
sequence.count('N') # 164553847
```

Compile all counts for R1.fwd and R1.rev in different libraries into one single file:

```python

names = ["Lib134_HeLa_siRNA3_Z_AP1_S5",
"Lib135_HeLa_siRNA3_Z_AP2_S1",
"Lib145_HeLa_siRNA3_Z_AP3_S1",
"Lib146_HeLa_siRNA_Z_AP4_S1",
"Lib136_HeLa_siRNA3_Z_Y1_S4",
"Lib137_HeLa_siRNA3_Z_Y2_S2",
"Lib147_HeLa_siRNA_Z_Y3_S2",
"Lib148_HeLa_siRNA_Z_Y4_S2",
"Lib141_HeLa_cp_Z_AP1_S1",
"Lib142_HeLa_cp_Z_AP2_S1",
"Lib149_HeLa_cp_Z_AP3_S1",
"Lib150_HeLa_cp_Z_AP4_S1",
"Lib143_HeLa_cp_Z_Y1_S2",
"Lib144_HeLa_cp_Z_Y2_S2",
"Lib151_HeLa_cp_Z_Y3_S2",
"Lib152_HeLa_cp_Z_Y4_S2"]

d = {}

for n in names:
  ifwd = open("~/base_composition/%s.clean.R1.fwd.cigar.txt" % n, "r")
  fwd = ifwd.readlines()
  ifwd.close()
  fwd_d = {}
  for line in fwd:
    fields=line.split()
    fwd_d[fields[1]] = fields[0]
  irev = open("~/base_composition/%s.clean.R1.rev.cigar.txt" % n, "r")
  rev = irev.readlines()
  irev.close()
  rev_d = {}
  for line in rev:
    fields=line.split()
    rev_d[fields[1]] = fields[0]
  n_short="_".join(n.split("_")[4:9])
  d[(n_short, "A")] = int(fwd_d["A"]) + int(rev_d["T"])
  d[(n_short, "C")] = int(fwd_d["C"]) + int(rev_d["G"])
  d[(n_short, "G")] = int(fwd_d["G"]) + int(rev_d["C"])
  d[(n_short, "T")] = int(fwd_d["T"]) + int(rev_d["A"])


len(d) # 72, ok 18*4


# Counts
ofile=open("~/base_composition/counts.txt", "w")
ofile.write("\t" + "\t".join(["A", "C", "G", "T"]) + "\n")

for n in names:
  n_short="_".join(n.split("_")[4:9])
  ofile.write("%s\t" % n_short)
  ofile.write("\t".join([str(d[n_short, l]) for l in ["A", "C", "G", "T"]]) + "\n")

ofile.close()


# Normalised counts
freq = {"A": 863155826, "C": 596433284, "G": 598488389, "T": 865655055} # obtained above

ofile=open("~/base_composition/normalisedcounts.txt", "w")
ofile.write("\t" + "\t".join(["A", "C", "G", "T"]) + "\n")

for n in names:
  n_short="_".join(n.split("_")[4:9])
  ofile.write("%s\t" % n_short)
  ofile.write("\t".join([str(round(float(d[n_short, l])/freq[l], 5)) for l in ["A", "C", "G", "T"]]) + "\n")

ofile.close()
```
