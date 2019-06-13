
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
- [Calling AP sites](#calling-ap-sites)



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

# Lib134 and Lib135
bedtools intersect \
-a Lib134_HeLa_siRNA3_Z_AP1_S5.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib135_HeLa_siRNA3_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 56780, 100*56780/148097 = 38%

bedtools intersect \
-a Lib135_HeLa_siRNA3_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib134_HeLa_siRNA3_Z_AP1_S5.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 56908, 100*56908/156856 = 36%


# Lib145 and Lib146
bedtools intersect \
-a Lib145_HeLa_siRNA3_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib146_HeLa_siRNA_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 34241, 100*34241/91364 = 37%

bedtools intersect \
-a Lib146_HeLa_siRNA_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib145_HeLa_siRNA3_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 35282, 100*35282/87183 = 40%


# Lib134 and Lib145
bedtools intersect \
-a Lib134_HeLa_siRNA3_Z_AP1_S5.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib145_HeLa_siRNA3_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 24709, 100*24709/148097 = 17%

bedtools intersect \
-a Lib145_HeLa_siRNA3_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib134_HeLa_siRNA3_Z_AP1_S5.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 23001, 100*23001/91364 = 25%


# Lib134 and Lib146
bedtools intersect \
-a Lib134_HeLa_siRNA3_Z_AP1_S5.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib146_HeLa_siRNA_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 35468, 100*35468/148097 = 24%

bedtools intersect \
-a Lib146_HeLa_siRNA_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib134_HeLa_siRNA3_Z_AP1_S5.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 32872, 100*32872/87183 = 38%


# Lib135 and Lib145
bedtools intersect \
-a Lib135_HeLa_siRNA3_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib145_HeLa_siRNA3_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 25512, 100*25512/156856 = 16%

bedtools intersect \
-a 180503_NS500222_0398_HTCJ3BGX5_Lib145_HeLa_siRNA3_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b 180427_NS500222_0395_HTCGYBGX5_Lib135_HeLa_siRNA3_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 23427, 100*23427/91364 = 26%


# Lib135 and Lib146
bedtools intersect \
-a Lib135_HeLa_siRNA3_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib146_HeLa_siRNA_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 36259, 100*36259/156856 = 23%

bedtools intersect \
-a Lib146_HeLa_siRNA_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib135_HeLa_siRNA3_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 33317, 100*33317/87183 = 38%


# multiinter
bedtools multiinter -i \
Lib134_HeLa_siRNA3_Z_AP1_S5.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib135_HeLa_siRNA3_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib145_HeLa_siRNA3_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib146_HeLa_siRNA_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-names Lib134 Lib135 Lib145 Lib146 | \
cut -f5 | sort | uniq -c | sort -k1,1nr
# 135975 Lib135
# 127583 Lib134
#  89546 Lib145
#  71165 Lib146
#  47875 Lib134,Lib135
#  37087 Lib145,Lib146
#  21113 Lib135,Lib146
#  20717 Lib134,Lib135,Lib146
#  20665 Lib134,Lib146
#  16428 Lib134,Lib135,Lib145,Lib146
#  14626 Lib135,Lib145,Lib146
#  13933 Lib134,Lib145,Lib146
#  10663 Lib135,Lib145
#  10106 Lib134,Lib145
#   8278 Lib134,Lib135,Lib145


# merge
tableCat.py -i \
Lib134_HeLa_siRNA3_Z_AP1_S5.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib135_HeLa_siRNA3_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib145_HeLa_siRNA3_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib146_HeLa_siRNA_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak | \
awk -v OFS="\t" '{print $1, $2, $3, $NF}' | \
sort -k1,1 -k2,2n | \
bedtools merge -c 4,4 -o distinct,count_distinct -i - | \
awk -v OFS="\t" '$5 > 2 {print $1, $2, $3}' > Lib134_Lib135_Lib145_Lib146_merge_p1e-5.bed

wc -l Lib134_Lib135_Lib145_Lib146_merge_p1e-5.bed # 27516



##############
# Cont-siRNA #
##############

# Lib141 and Lib142
bedtools intersect \
-a Lib141_HeLa_cp_Z_AP1_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib142_HeLa_cp_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 37664, 100*37664/59705 = 63%

bedtools intersect \
-a Lib142_HeLa_cp_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib141_HeLa_cp_Z_AP1_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 36459, 100*36459/93948 = 39%


# Lib149 and Lib150
bedtools intersect \
-a Lib149_HeLa_cp_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib150_HeLa_cp_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 7370, 100*7370/18639 = 40%

bedtools intersect \
-a Lib150_HeLa_cp_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib149_HeLa_cp_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 7676, 100*7676/92196 = 8%


# Lib141 and Lib149
bedtools intersect \
-a Lib141_HeLa_cp_Z_AP1_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib149_HeLa_cp_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 9508, 100*9508/59705 = 16%

bedtools intersect \
-a Lib149_HeLa_cp_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib141_HeLa_cp_Z_AP1_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 9486, 100*9486/18639 = 51%


# Lib141 and Lib150
bedtools intersect \
-a Lib141_HeLa_cp_Z_AP1_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib150_HeLa_cp_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 17559, 100*17559/59705 = 29%

bedtools intersect \
-a Lib150_HeLa_cp_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib141_HeLa_cp_Z_AP1_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 17508, 100*17508/92196 = 19%


# Lib142 and Lib149
bedtools intersect \
-a Lib142_HeLa_cp_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib149_HeLa_cp_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 11735, 100*11735/93948 = 12%

bedtools intersect \
-a Lib149_HeLa_cp_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib142_HeLa_cp_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 11678, 100*11678/18639 = 63%


# Lib142 and Lib150
bedtools intersect \
-a Lib142_HeLa_cp_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib150_HeLa_cp_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 22707, 100*22707/93948 = 24%

bedtools intersect \
-a Lib150_HeLa_cp_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib142_HeLa_cp_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 22434, 100*22434/92196 = 24%


# multiinter
bedtools multiinter -i \
Lib141_HeLa_cp_Z_AP1_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib142_HeLa_cp_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib149_HeLa_cp_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib150_HeLa_cp_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
-names Lib141 Lib142 Lib149 Lib150 | \
cut -f5 | sort | uniq -c | sort -k1,1nr
#  91801 Lib142
#  87714 Lib150
#  47181 Lib141
#  39526 Lib141,Lib142
#  19523 Lib142,Lib150
#  12831 Lib141,Lib142,Lib150
#  12386 Lib149
#  11474 Lib141,Lib150
#   7414 Lib142,Lib149
#   7239 Lib141,Lib142,Lib149
#   4244 Lib141,Lib142,Lib149,Lib150
#   4186 Lib141,Lib149
#   3929 Lib149,Lib150
#   2937 Lib142,Lib149,Lib150
#   1697 Lib141,Lib149,Lib150


# merge
tableCat.py -i \
Lib141_HeLa_cp_Z_AP1_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib142_HeLa_cp_Z_AP2_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib149_HeLa_cp_Z_AP3_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib150_HeLa_cp_Z_AP4_S1.hg38.clean.nomodel.p1e-5_peaks.narrowPeak | \
awk -v OFS="\t" '{print $1, $2, $3, $NF}' | \
sort -k1,1 -k2,2n | \
bedtools merge -c 4,4 -o distinct,count_distinct -i - | \
awk -v OFS="\t" '$5 > 2 {print $1, $2, $3}' > Lib141_Lib142_Lib149_Lib150_merge_p1e-5.bed

wc -l Lib141_Lib142_Lib149_Lib150_merge_p1e-5.bed # 16835
```

Intersections of `Lib134_Lib135_Lib145_Lib146_merge_p1e-5.bed` and `Lib141_Lib142_Lib149_Lib150_merge_p1e-5.bed`:

```bash
cd ~/macs2_hg38/default

wc -l Lib134_Lib135_Lib145_Lib146_merge_p1e-5.bed Lib141_Lib142_Lib149_Lib150_merge_p1e-5.bed
#  27516 Lib134_Lib135_Lib145_Lib146_merge_p1e-5.bed
#  16835 Lib141_Lib142_Lib149_Lib150_merge_p1e-5.bed

bedtools intersect \
-a Lib134_Lib135_Lib145_Lib146_merge_p1e-5.bed \
-b Lib141_Lib142_Lib149_Lib150_merge_p1e-5.bed \
-sorted \
-wa -u | wc -l # 12492 Lib134_Lib135_Lib145_Lib146_merge_p1e-5 regions intersect

bedtools intersect \
-a Lib141_Lib142_Lib149_Lib150_merge_p1e-5.bed \
-b Lib134_Lib135_Lib145_Lib146_merge_p1e-5.bed \
-sorted \
-wa -u | wc -l # 12606 (100*12606/16835 = 75%) Lib141_Lib142_Lib149_Lib150_merge_p1e-5 regions intersect
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


### PQS

```bash
cd ~/annotation/

ref=~/reference/genome.fa

nohup fastaRegexFinder.py -f $ref -q | \
bedtools sort -i | \
cut -f1-3 | \
grep -v -E "_|chrEBV" | \
sed 's/$/&\tPQS/' > hg38_pqs.bed &

wc -l hg38_pqs.bed # 362816
```


### combine gene features, intergenic and PQS

```bash
cd ~/annotation

cat <(cut -f1-3,6 hg38_promoter_gene.bed hg38_utr5_transcript.bed hg38_exon_transcript.bed hg38_intron_transcript.bed hg38_utr3_transcript.bed) \
hg38_intergenic.bed \
hg38_pqs.bed | \
bedtools sort -i - > hg38_promoter_utr5_exon_intron_utr3_intergenic_pqs.bed
```


### define mappable genome

Need to define two mappable files: APE1 k.d. and WT:

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

annotations=~/annotation/hg38_promoter_utr5_exon_intron_utr3_intergenic_pqs.bed

# APE1-siRNA
mappable=~/bam/APE1.Y.mappable.bed

bed=Lib134_Lib135_Lib145_Lib146_merge_p1e-5.bed
bname=`basename ${bed%.bed}`
gat-run.py -a $annotations -s <(grep -v -E '_|chrEBV' $bed) -w $mappable --ignore-segment-tracks -n 10000 -t 20 -L ../../gat/$bname.log > ../../gat/$bname.txt


# Cont-siRNA
mappable=~/bam/WT.Y.mappable.bed

bed=Lib141_Lib142_Lib149_Lib150_merge_p1e-5.bed
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
"Lib169_HeLa_cp_Z_AP5_S3",
"Lib143_HeLa_cp_Z_Y1_S2",
"Lib144_HeLa_cp_Z_Y2_S2",
"Lib151_HeLa_cp_Z_Y3_S2",
"Lib152_HeLa_cp_Z_Y4_S2",
"Lib171_HeLa_cp_Z_Y5_S2"]

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



## Calling AP sites

### R1.fwd APE1-siRNA

#### Coverage of alignment start sites

```bash
cd ~/bam

mkdir -p ../ap/cov


### Coverage by chromosome
g=~/reference/genome_spikeins.fa.fai

for bam in Lib134_HeLa_siRNA3_Z_AP1_S5.hg38.clean.R1.fwd.cigar.bam \
Lib135_HeLa_siRNA3_Z_AP2_S1.hg38.clean.R1.fwd.cigar.bam \
Lib145_HeLa_siRNA3_Z_AP3_S1.hg38.clean.R1.fwd.cigar.bam \
Lib146_HeLa_siRNA_Z_AP4_S1.hg38.clean.R1.fwd.cigar.bam \
Lib136_HeLa_siRNA3_Z_Y1_S4.hg38.clean.R1.fwd.cigar.bam \
Lib137_HeLa_siRNA3_Z_Y2_S2.hg38.clean.R1.fwd.cigar.bam \
Lib147_HeLa_siRNA_Z_Y3_S2.hg38.clean.R1.fwd.cigar.bam \
Lib148_HeLa_siRNA_Z_Y4_S2.hg38.clean.R1.fwd.cigar.bam
do
  bname=`basename ${bam%.bam}`
  for chr in `cut -f1 $g | grep 'chr*' | grep -v -E "_|EBV"`
  do
    bedtools genomecov -d -5 -ibam $bam -g $g | grep -P '^$chr\t' | pigz -p 20 > ../ap/cov/$bname.$chr.cov.gz
  done
done


### Concatenate by chromosome
cd ../ap/cov

for chr in `cut -f1 $g | grep 'chr*' | grep -v -E "_|EBV"`
do
  paste \
  <(zcat Lib134_HeLa_siRNA3_Z_AP1_S5.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib135_HeLa_siRNA3_Z_AP2_S1.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib145_HeLa_siRNA3_Z_AP3_S1.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib146_HeLa_siRNA_Z_AP4_S1.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib136_HeLa_siRNA3_Z_Y1_S4.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib137_HeLa_siRNA3_Z_Y2_S2.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib147_HeLa_siRNA_Z_Y3_S2.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib148_HeLa_siRNA_Z_Y4_S2.hg38.clean.R1.fwd.cigar.$chr.cov.gz) | \
  cut -f1-3,6,9,12,15,18,21,24 | \
  pigz -p 20 > APE1.R1.fwd.cigar.$chr.cov.gz
done
```


#### Obtain sites

```bash
cd ~/ap/cov/

mkdir ../tables

for cov in APE1.R1.fwd.cigar.*.cov.gz; do
  bname=${cov%.gz}
  Rscript --vanilla ../../scripts/20180725_apcaller_human.R $cov
done
```

chrM separately:

```r
library(data.table)
library(edgeR)

# Load data
print("Load data")
data <- fread("zcat APE1.R1.fwd.cigar.chrM.cov.gz")
setnames(data, c("chr", "pos", "Lib134", "Lib135", "Lib145", "Lib146", "Lib136", "Lib137", "Lib147", "Lib148"))

# Define group
group <- factor(c('ap', 'ap', 'ap', 'ap', 'input', 'input', 'input', 'input'))

# Define DGEList object
print("Define DGEList object")
y <- DGEList(counts = data[,-c(1,2)], group = group, genes = data[,c(1,2)])

# Filter and get the top 100k sites according to cpm
print("Filter and get the top 100k sites according to cpm")
if(grepl("chrM", "APE1.R1.fwd.cigar.chrM.cov.gz")){
y_f <- y
} else {
y_f <- y[names(sort(rowMeans(cpm(y)), decreasing=T)[1:1e5]),]
}

# Define design matrix
des <- model.matrix(~ 0 + group, data = y_f$samples)
colnames(des) <- levels(factor(y_f$samples$group))

# Calculate normalization factors
y_f <- calcNormFactors(y_f, method = "TMM")

# Estimate dispersion
print("Estimate dispersion")
y_f <- estimateDisp(y_f, des, min.row.sum=7)

# Fit linear model
fit <- glmFit(y_f, des)

# Define contrasts
my.contrasts <- makeContrasts(apVSinput = ap - input, levels = des)

# Obtain likelihoods
lrt <- glmLRT(fit, contrast=my.contrasts[,"apVSinput"])

# Table
detable <- topTags(lrt, n = Inf)$table
detable_e <- data.table(cbind(detable, data.frame(data[,-c(1,2)][as.numeric(row.names(detable)),])))
detable_e[, start := as.integer(pos - 2)]
detable_e[, end := as.integer(start + 1)]
detable_e <- detable_e[,.(chr, start, end, Lib134, Lib135, Lib145, Lib146, Lib136, Lib137, Lib147, Lib148, logFC, logCPM, LR, PValue, FDR)]
print("Writing table")
write.table(detable_e, "../tables/APE1.R1.fwd.cigar.chrM.cov.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
```

Plotting:

```r
library(data.table)
library(ggplot2)

# Enlarge the view width when printing tables
options(width = 250)

# Concatenate tables
data <- fread('tableCat.py -i ~/ap/tables/APE1.R1.fwd.cigar.*.cov.txt -H')

# Explore
nrow(data[FDR < 1e-1 & logFC > 0]) # 417
nrow(data[FDR < 1e-1 & logFC < 0]) # 2293
nrow(data[FDR < 0.05 & logFC > 0]) # 229
nrow(data[FDR < 0.05 & logFC < 0]) # 1512
nrow(data[FDR < 0.01 & logFC > 0]) # 86
nrow(data[FDR < 0.01 & logFC < 0]) # 692
data[logFC > 0][1,]$FDR #
data[logFC < 0][1,]$FDR # 4.021293e-16

# Volcano plot
df <- data.table(x = data$logFC, y = -log10(data$FDR), d = densCols(data$logFC, -log10(data$FDR), colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

gg <- ggplot(data = df, aes(x, y, col = d)) +
geom_point(size = 0.1) +
scale_color_identity() +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
geom_vline(xintercept = 0, linetype = 'dotted') +
theme_bw() +
coord_cartesian(xlim = c(-7, 7))

ggsave('~/figures/APE1.R1.fwd.cigar.all.png', width = 12, height = 12, units = 'cm')
```


### R1.fwd Cont-siRNA

#### Coverage of alignment start sites

```bash
cd ~/bam

mkdir -p ../ap/cov


### Coverage by chromosome
g=~/reference/genome_spikeins.fa.fai

for bam in Lib141_HeLa_cp_Z_AP1_S1.hg38.clean.R1.fwd.cigar.bam \
Lib142_HeLa_cp_Z_AP2_S1.hg38.clean.R1.fwd.cigar.bam \
Lib149_HeLa_cp_Z_AP3_S1.hg38.clean.R1.fwd.cigar.bam \
Lib150_HeLa_cp_Z_AP4_S1.hg38.clean.R1.fwd.cigar.bam \
Lib143_HeLa_cp_Z_Y1_S2.hg38.clean.R1.fwd.cigar.bam \
Lib144_HeLa_cp_Z_Y2_S2.hg38.clean.R1.fwd.cigar.bam \
Lib151_HeLa_cp_Z_Y3_S2.hg38.clean.R1.fwd.cigar.bam \
Lib152_HeLa_cp_Z_Y4_S2.hg38.clean.R1.fwd.cigar.bam
do
  bname=`basename ${bam%.bam}`
  for chr in `cut -f1 $g | grep 'chr*' | grep -v -E "_|EBV"`
  do
    bedtools genomecov -d -5 -ibam $bam -g $g | grep -P '^$chr\t' | pigz -p 20 > ../ap/cov/$bname.$chr.cov.gz
  done
done


### Concatenate by chromosome
cd ../ap/cov

for chr in `cut -f1 $g | grep 'chr*' | grep -v -E "_|EBV"`
do
  nohup paste \
  <(zcat Lib141_HeLa_cp_Z_AP1_S1.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib142_HeLa_cp_Z_AP2_S1.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib149_HeLa_cp_Z_AP3_S1.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib150_HeLa_cp_Z_AP4_S1.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib143_HeLa_cp_Z_Y1_S2.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib144_HeLa_cp_Z_Y2_S2.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib151_HeLa_cp_Z_Y3_S2.hg38.clean.R1.fwd.cigar.$chr.cov.gz) \
  <(zcat Lib152_HeLa_cp_Z_Y4_S2.hg38.clean.R1.fwd.cigar.$chr.cov.gz) | \
  cut -f1-3,6,9,12,15,18,21,24 | \
  pigz -p 20 > WT.R1.fwd.cigar.$chr.cov.gz &
done
```


#### Obtain sites

```bash
cd ~/ap/cov/

for cov in WT.R1.fwd.cigar.*.cov.gz; do
  bname=${cov%.gz}
  Rscript --vanilla ../../scripts/20180727_apcaller_human.R $cov
done
```

chr1, chr2 and chr3 separately:

```bash
cd ~/ap/cov/

for cov in WT.R1.fwd.cigar.chr{1..3}.cov.gz; do
  bname=${cov%.gz}
  Rscript --vanilla ../../scripts/20180727_apcaller_human.R $cov
done

tail ../tables/*.log # all OK, just chrM left
```

chrM:

```r
library(data.table)
library(edgeR)

# Load data
print("Load data")
data <- fread("zcat WT.R1.fwd.cigar.chrM.cov.gz")
setnames(data, c("chr", "pos", "Lib141", "Lib142", "Lib149", "Lib150", "Lib143", "Lib144", "Lib151", "Lib152"))

# Define group
group <- factor(c('ap', 'ap', 'ap', 'ap', 'input', 'input', 'input', 'input'))

# Define DGEList object
print("Define DGEList object")
y <- DGEList(counts = data[,-c(1,2)], group = group, genes = data[,c(1,2)])

# Filter and get the top 100k sites according to cpm
print("Filter and get the top 100k sites according to cpm")
if(grepl("chrM", "WT.R1.fwd.cigar.chrM.cov.gz")){
y_f <- y
} else {
y_f <- y[names(sort(rowMeans(cpm(y)), decreasing=T)[1:1e5]),]
}

# Define design matrix
des <- model.matrix(~ 0 + group, data = y_f$samples)
colnames(des) <- levels(factor(y_f$samples$group))

# Calculate normalization factors
y_f <- calcNormFactors(y_f, method = "TMM")

# Estimate dispersion
print("Estimate dispersion")
y_f <- estimateDisp(y_f, des, min.row.sum=7)

# Fit linear model
fit <- glmFit(y_f, des)

# Define contrasts
my.contrasts <- makeContrasts(apVSinput = ap - input, levels = des)

# Obtain likelihoods
lrt <- glmLRT(fit, contrast=my.contrasts[,"apVSinput"])

# Table
detable <- topTags(lrt, n = Inf)$table
detable_e <- data.table(cbind(detable, data.frame(data[,-c(1,2)][as.numeric(row.names(detable)),])))
detable_e[, start := as.integer(pos - 2)]
detable_e[, end := as.integer(start + 1)]
detable_e <- detable_e[,.(chr, start, end, Lib141, Lib142, Lib149, Lib150, Lib143, Lib144, Lib151, Lib152, logFC, logCPM, LR, PValue, FDR)]
print("Writing table")
write.table(detable_e, "../tables/WT.R1.fwd.cigar.chrM.cov.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
```

Plotting:

```r
library(data.table)
library(ggplot2)

# Enlarge the view width when printing tables
options(width = 250)

# Concatenate tables
data <- fread('tableCat.py -i ~/ap/tables/WT.R1.fwd.cigar.*.cov.txt -H')

# Explore
nrow(data[FDR < 1e-1 & logFC > 0]) # 511
nrow(data[FDR < 1e-1 & logFC < 0]) # 2376
nrow(data[FDR < 0.05 & logFC > 0]) # 335
nrow(data[FDR < 0.05 & logFC < 0]) # 1689
nrow(data[FDR < 0.01 & logFC > 0]) # 138
nrow(data[FDR < 0.01 & logFC < 0]) # 925
data[order(FDR)][logFC > 0][1,]$FDR # 2.559161e-16
data[order(FDR)][logFC < 0][1,]$FDR # 1.173508e-21

# Volcano plot
df <- data.table(x = data$logFC, y = -log10(data$FDR), d = densCols(data$logFC, -log10(data$FDR), colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

gg <- ggplot(data = df, aes(x, y, col = d)) +
geom_point(size = 0.1) +
scale_color_identity() +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
geom_vline(xintercept = 0, linetype = 'dotted') +
theme_bw() +
coord_cartesian(xlim = c(-7, 7))

ggsave('~/figures/WT.R1.fwd.cigar.all.png', width = 12, height = 12, units = 'cm')
```
