
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
- [Calling AP sites](#calling-ap-sites)
- [Comparison of SMUG1-snAP-seq sites with peaks obtained in Kawasaki et al., Genome biology, 2017](#comparison-of-smug1-snap-seq-sites-with-peaks-obtained-in-kawasaki-et-al-genome-biology-2017)
- [Merge technical replicates and calculate coverage](#merge-technical-replicates-and-calculate-coverage)
- [Base composition profiles](#base-composition-profiles)
- [Coverage profiles centered around the sites](#coverage-profiles-centered-around-the-sites)
- [Sequence logos](#sequence-logos)
- [Exploring TG enrichment using frequencies of dinucleotides counts, motifs and tests](#exploring-tg-enrichment-using-frequencies-of-dinucleotides-counts-motifs-and-tests)
- [Peak calling](#peak-calling)
- [Consensus peaks](#consensus-peaks)
- [Overlap with 5hmU and baseJ peaks](#overlap-with-5hmu-and-basej-peaks)



## Software requirements

- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [sambamba v0.6.5](https://academic.oup.com/bioinformatics/article/31/12/2032/214758)
- Standard Unix tools: awk, sort, uniq
- [igvtools v2.3.91](https://software.broadinstitute.org/software/igv/igvtools)
- [deeptools v2.4.2-5-f439d22](https://deeptools.readthedocs.io/en/develop/index.html)
- [bedtools v2.27.0](http://bedtools.readthedocs.io/en/latest/)
- [tableCat.py](https://github.com/dariober/bioinformatics-cafe/blob/master/tableCat/tableCat.py)
- [EMBOSS v6.6.0.0](http://emboss.sourceforge.net/)
- [fastaRegexFinder.py v0.1.1](https://github.com/dariober/bioinformatics-cafe/tree/master/fastaRegexFinder)
- [meme v4.11.2](http://meme-suite.org/)
- [macs2 v2.1.1.20160309](https://github.com/taoliu/MACS)
- [python v2.7.12](https://www.python.org/). Libraries:
  - [os](https://docs.python.org/2/library/os.html)
  - [collections](https://docs.python.org/2/library/collections.html)
  - [random](https://docs.python.org/2/library/random.html)
  - [scipy](https://www.scipy.org/)
  - [numpy](http://www.numpy.org/)
  - [gzip](https://docs.python.org/2.7/library/gzip.html)
  - [re](https://docs.python.org/2/library/re.html)
  - [pandas](https://pandas.pydata.org/)
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [ggplot2 v2.2.1](http://ggplot2.org/)
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)
  - [edgeR v3.16.5](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
  - [ggseqlogo v0.1](https://cran.rstudio.com/web/packages/ggseqlogo/index.html)
  - [Biostrings v2.42.1](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
  - [GenomicFeatures v1.26.4](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)



## Libraries

SMUG1 batch:

Library | Biological replicate
:------:|:-------------------:
SMUG1-snAP-seq | rep1
SMUG1-snAP-seq | rep2
snAP-seq | rep1
snAP-seq | rep2
Y-input | rep1
Y-input | rep2

UNG batch:

Library | Biological replicate
:------:|:-------------------:
UNG-snAP-seq | rep1
UNG-snAP-seq | rep2
Y-input | rep1
Y-input | rep2



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
  fq2=${fq1/_R1_/_R2_}
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
  fq2=${fq1/_R1_/_R2_}
  bname=${fq1%_R1_001.fastq.gz}
  cutadapt -g GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -G GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -e 0.1 -O 20 --discard-trimmed --pair-filter=any -o ../fastq_trimmed_protecting_sequence_optimal/$fq1 -p ../fastq_trimmed_protecting_sequence_optimal/$fq2 $fq1 $fq2 > ../fastq_trimmed_protecting_sequence_optimal/$bname.txt
done
```



## Alignment

### Prepare and index references

```bash
cd ~/reference/

wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/Leishmania/major/current_gff3/Lmajor.genome.fasta.gz
wget ftp://ftp.sanger.ac.uk/pub/project/pathogens/Leishmania/major/current_gff3/Lmajor.gff3.gz
gzip -d Lmajor.genome.fasta.gz

cat Lmajor.genome.fasta spikeins.fa > Lmajor.genome_spikeins.fasta

bwa index Lmajor.genome_spikeins.fasta
samtools faidx Lmajor.genome_spikeins.fasta
```


### Align and sort

```bash
cd ~/fastq_trimmed_protecting_sequence_optimal

mkdir ../bam

ref=../reference/Lmajor.genome_spikeins.fasta

for fq1 in *R1*.fastq.gz
do
  bname=${fq1%_R1_001.fastq.gz}
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

listOfIds="..." # add here the list of ids for the libraries you want to merge

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

```bash
cd ~/bam

for bam in `*.bam | grep -v "tmp"`
do
  bname=${bam%.bam}
  samtools idxstats $bam | cut -f1 | grep -v -E 'AP1|AP2|fC1|fC2|fU1|fU2|GCAT1|GCAT2' | xargs samtools view -@ 20 -F 3844 -q 10 -b $bam > $bname.clean.bam && samtools index $bname.clean.bam && \
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
  igvtools count $bam ../tdf/$bname.tdf ../reference/Lmajor.genome_spikeins.fasta.fai --includeDuplicates --pairs -w 1 -e 0 --preExtFactor 0 --postExtFactor 0
done

#single-end
for bam in *.clean.R*.bam
do
  bname=${bam%.bam}
  igvtools count $bam ../tdf/$bname.tdf ../reference/Lmajor.genome_spikeins.fasta.fai --includeDuplicates -w 1 -e 0 --preExtFactor 0 --postExtFactor 0
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



## Calling AP sites

### R1.fwd

#### Coverage of alignment start sites

```bash
cd ~/bam

mkdir ../ap

g=../reference/Lmajor.genome_spikeins.fasta.fai

for bam in *.clean.R1.fwd.cigar.bam
do
  bname=`basename ${bam%.bam}`
  bedtools genomecov -d -5 -ibam $bam -g $g | grep 'LmjF*' > ../ap/$bname.cov
done
```


#### Obtain sites

```r
library(data.table)
library(edgeR)
library(ggplot2)

# Enlarge the view width when printing tables
options(width = 250)


##########################################
# Lib96 and Lib100 vs. Lib99b and Lib103 #
##########################################

# Load data
data <- fread("tableCat.py -i ~/ap/*Lib{96,100,99b,103}*fwd*.cov -r .clean.R1.fwd.cigar.cov")
setnames(data, c("chr", "pos", "count", "lib"))
table(data$lib)

# Cast
data_cast <- dcast.data.table(data = data, chr + pos ~ lib, value.var = "count")

# Define group
group <- factor(c('ap', 'input', 'ap', 'input'))

# Define DGEList object
y <- DGEList(counts = data_cast[,-c(1,2)], group = group, genes = data_cast[,c(1,2)])
y$samples

# Filter and get the top 100k sites according to cpm
y_f <- y[names(sort(rowMeans(cpm(y)), decreasing=T)[1:1e5]),]

# Define design matrix
des <- model.matrix(~ 0 + group, data = y_f$samples)
colnames(des) <- levels(factor(y_f$samples$group))

# Calculate normalization factors
y_f <- calcNormFactors(y_f, method = "TMM")

# Estimate dispersion
y_f <- estimateDisp(y_f, des)

# Fit linear model
fit <- glmFit(y_f, des, prior.count=1)

# Define contrasts
my.contrasts <- makeContrasts(apVSinput = ap - input, levels = des)

# Obtain likelihoods
lrt <- glmLRT(fit, contrast=my.contrasts[,"apVSinput"])

# Tables, bed files and volcano plots
## tables
detable <- topTags(lrt, n = Inf)$table
detable_e <- data.table(cbind(detable, data.frame(data_cast[,-c(1,2)][as.numeric(row.names(detable)),])))
detable_e[, start := as.integer(pos - 2)]
detable_e[, end := as.integer(start + 1)]
detable_e <- detable_e[,.(chr, start, end, Lib96, Lib100, Lib99b, Lib103, logFC, logCPM, LR, PValue, FDR)]

nrow(detable_e[FDR < 1e-15 & logFC > 0]) # 492 - very high confidence sites
nrow(detable_e[FDR < 1e-15 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-10 & logFC > 0]) # 1560 - high confidence sites
nrow(detable_e[FDR < 1e-10 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-5 & logFC > 0]) # 5635 - sites of confidence
nrow(detable_e[FDR < 1e-5 & logFC < 0]) # 6
nrow(detable_e[FDR < 1e-1 & logFC > 0]) # 38575 - more sites
nrow(detable_e[FDR < 1e-1 & logFC < 0]) # 3354
detable_e[logFC < 0][1,]$FDR # 2.086796e-09 - the best FDR in the opposite direction logFC<0 could be a good threshold for our logFC>0

write.table(detable_e[FDR < 1e-10 & logFC > 0], "~/ap/Lib96_Leish_SMUG1-AP_S4.Lib100_Leish_SMUG1-AP_2_S7.clean.R1.fwd.cigar_FDR1e-10.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

## bed
write.table(detable_e[FDR < 1e-10 & logFC > 0][, c("chr", "start", "end", "FDR")][start>-1], "~/ap/Lib96_Leish_SMUG1-AP_S4.Lib100_Leish_SMUG1-AP_2_S7.clean.R1.fwd.cigar_FDR1e-10.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

## volcano
df <- data.table(x = detable_e$logFC, y = -log10(detable_e$FDR), d = densCols(detable_e$logFC, -log10(detable_e$FDR), colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

gg <- ggplot(data = df, aes(x, y, col = d)) +
geom_point(size = 0.1) +
scale_color_identity() +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
geom_vline(xintercept = 0, linetype = 'dotted') +
theme_bw() +
coord_cartesian(xlim = c(-10, 10))

ggsave('~/figures/Lib96_Leish_SMUG1-AP.Lib100_Leish_SMUG1-AP_2.clean.R1.fwd.cigar.png', width = 12, height = 12, units = 'cm')


##########################################
# Lib97 and Lib101 vs. Lib99b and Lib103 #
##########################################

# Load data
data <- fread("tableCat.py -i ~/ap/*Lib{97,101,99b,103}*fwd*.cov -r .clean.R1.fwd.cigar.cov")
setnames(data, c("chr", "pos", "count", "lib"))
table(data$lib)

# Cast
data_cast <- dcast.data.table(data = data, chr + pos ~ lib, value.var = "count")

# Define group
group <- factor(c('ap', 'input', 'ap', 'input'))

# Define DGEList object
y <- DGEList(counts = data_cast[,-c(1,2)], group = group, genes = data_cast[,c(1,2)])
y$samples

# Filter and get the top 100k sites according to cpm
y_f <- y[names(sort(rowMeans(cpm(y)), decreasing=T)[1:1e5]),]

# Define design matrix
des <- model.matrix(~ 0 + group, data = y_f$samples)
colnames(des) <- levels(factor(y_f$samples$group))

# Calculate normalization factors
y_f <- calcNormFactors(y_f, method = "TMM")

# Estimate dispersion
y_f <- estimateDisp(y_f, des)

# Fit linear model
fit <- glmFit(y_f, des, prior.count=1)

# Define contrasts
my.contrasts <- makeContrasts(apVSinput = ap - input, levels = des)

# Obtain likelihoods
lrt <- glmLRT(fit, contrast=my.contrasts[,"apVSinput"])

# Tables, bed files and volcano plots
## tables
detable <- topTags(lrt, n = Inf)$table
detable_e <- data.table(cbind(detable, data.frame(data_cast[,-c(1,2)][as.numeric(row.names(detable)),])))
detable_e[, start := as.integer(pos - 2)]
detable_e[, end := as.integer(start + 1)]
detable_e <- detable_e[,.(chr, start, end, Lib97, Lib101, Lib99b, Lib103, logFC, logCPM, LR, PValue, FDR)]

nrow(detable_e[FDR < 1e-15 & logFC > 0]) # 0 - very high confidence sites
nrow(detable_e[FDR < 1e-15 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-10 & logFC > 0]) # 0 - high confidence sites
nrow(detable_e[FDR < 1e-10 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-5 & logFC > 0]) # 0 - sites of confidence
nrow(detable_e[FDR < 1e-5 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-1 & logFC > 0]) # 1 - more sites
nrow(detable_e[FDR < 1e-1 & logFC < 0]) # 3
detable_e[logFC < 0][1,]$FDR # 0.01005186

## volcano
df <- data.table(x = detable_e$logFC, y = -log10(detable_e$FDR), d = densCols(detable_e$logFC, -log10(detable_e$FDR), colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

gg <- ggplot(data = df, aes(x, y, col = d)) +
geom_point(size = 0.1) +
scale_color_identity() +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
geom_vline(xintercept = 0, linetype = 'dotted') +
theme_bw() +
coord_cartesian(xlim = c(-10, 10))

ggsave('~/figures/Lib97_Leish_AP.Lib101_Leish_AP_2.clean.R1.fwd.cigar.png', width = 12, height = 12, units = 'cm')


############################
# UNG-snAP-seq vs. Y-input #
############################

# Load data
data <- fread("tableCat.py -i ~/ap/*Lib13{0..3}*fwd*cov -r .clean.R1.fwd.cigar.cov")
setnames(data, c("chr", "pos", "count", "lib"))
table(data$lib)

# Cast
data_cast <- dcast.data.table(data = data, chr + pos ~ lib, value.var = "count")

# Define group
group <- factor(c('ap', 'ap', 'input', 'input'))

# Define DGEList object
y <- DGEList(counts = data_cast[,-c(1,2)], group = group, genes = data_cast[,c(1,2)])
y$samples

# Filter and get the top 100k sites according to cpm
y_f <- y[names(sort(rowMeans(cpm(y)), decreasing=T)[1:1e5]),]

# Define design matrix
des <- model.matrix(~ 0 + group, data = y_f$samples)
colnames(des) <- levels(factor(y_f$samples$group))

# Calculate normalization factors
y_f <- calcNormFactors(y_f, method = "TMM")

# Estimate dispersion
y_f <- estimateDisp(y_f, des)

# Fit linear model
fit <- glmFit(y_f, des, prior.count=1)

# Define contrasts
my.contrasts <- makeContrasts(apVSinput = ap - input, levels = des)

# Obtain likelihoods
lrt <- glmLRT(fit, contrast=my.contrasts[,"apVSinput"])

# Tables, bed files and volcano plots
## tables
detable <- topTags(lrt, n = Inf)$table
detable_e <- data.table(cbind(detable, data.frame(data_cast[,-c(1,2)][as.numeric(row.names(detable)),])))
detable_e[, start := as.integer(pos - 2)]
detable_e[, end := as.integer(start + 1)]
detable_e <- detable_e[,.(chr, start, end, Lib130, Lib131, Lib132, Lib133, logFC, logCPM, LR, PValue, FDR)]

nrow(detable_e[FDR < 1e-15 & logFC > 0]) # 0 - very high confidence sites
nrow(detable_e[FDR < 1e-15 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-10 & logFC > 0]) # 0 - high confidence sites
nrow(detable_e[FDR < 1e-10 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-5 & logFC > 0]) # 0 - sites of confidence
nrow(detable_e[FDR < 1e-5 & logFC < 0]) # 7
nrow(detable_e[FDR < 1e-1 & logFC > 0]) # 21 - more sites
nrow(detable_e[FDR < 1e-1 & logFC < 0]) # 136
detable_e[logFC < 0][1,]$FDR # 2.588202e-09 - the best FDR in the opposite direction logFC<0 could be a good threshold for our logFC>0

## volcano
df <- data.table(x = detable_e$logFC, y = -log10(detable_e$FDR), d = densCols(detable_e$logFC, -log10(detable_e$FDR), colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

gg <- ggplot(data = df, aes(x, y, col = d)) +
geom_point(size = 0.1) +
scale_color_identity() +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
geom_vline(xintercept = 0, linetype = 'dotted') +
theme_bw() +
coord_cartesian(xlim = c(-10, 10))

ggsave('~/figures/Lib130_Leish_UNG-AP1.Lib131_Leish_UNG-AP2.clean.R1.fwd.cigar.png', width = 12, height = 12, units = 'cm')
```


### R1.rev

#### Coverage of alignment start sites

```bash
cd ~/bam

g=../reference/Lmajor.genome_spikeins.fasta.fai

for bam in *Leish*.clean.R1.rev.cigar.bam
do
  bname=`basename ${bam%.bam}`
  bedtools genomecov -d -5 -ibam $bam -g $g | grep 'LmjF*' > ../ap/$bname.cov
done
```


#### Obtain sites

```r
library(data.table)
library(edgeR)
library(ggplot2)

# Enlarge the view width when printing tables
options(width = 250)


##############################
# SMUG1-snAP-seq vs. Y-input #
##############################

# Load data
data <- fread("tableCat.py -i ~/ap/*Lib{96,100,99b,103}*rev*.cov -r .clean.R1.rev.cigar.cov")
setnames(data, c("chr", "pos", "count", "lib"))
table(data$lib)

# Cast
data_cast <- dcast.data.table(data = data, chr + pos ~ lib, value.var = "count")

# Define group
group <- factor(c('ap', 'input', 'ap', 'input'))

# Define DGEList object
y <- DGEList(counts = data_cast[,-c(1,2)], group = group, genes = data_cast[,c(1,2)])
y$samples

# Filter and get the top 100k sites according to cpm
y_f <- y[names(sort(rowMeans(cpm(y)), decreasing=T)[1:1e5]),]

# Define design matrix
des <- model.matrix(~ 0 + group, data = y_f$samples)
colnames(des) <- levels(factor(y_f$samples$group))

# Calculate normalization factors
y_f <- calcNormFactors(y_f, method = "TMM")

# Estimate dispersion
y_f <- estimateDisp(y_f, des)

# Fit linear model
fit <- glmFit(y_f, des, prior.count=1)

# Define contrasts
my.contrasts <- makeContrasts(apVSinput = ap - input, levels = des)

# Obtain likelihoods
lrt <- glmLRT(fit, contrast=my.contrasts[,"apVSinput"])

# Tables, bed files and volcano plots
## tables
detable <- topTags(lrt, n = Inf)$table
detable_e <- data.table(cbind(detable, data.frame(data_cast[,-c(1,2)][as.numeric(row.names(detable)),])))
detable_e[, start := as.integer(pos)]
detable_e[, end := as.integer(start + 1)]
detable_e <- detable_e[,.(chr, start, end, Lib96, Lib100, Lib99b, Lib103, logFC, logCPM, LR, PValue, FDR)]

nrow(detable_e[FDR < 1e-15 & logFC > 0]) # 545 - very high confidence sites
nrow(detable_e[FDR < 1e-15 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-10 & logFC > 0]) # 1640 - high confidence sites
nrow(detable_e[FDR < 1e-10 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-5 & logFC > 0]) # 6081 - sites of confidence
nrow(detable_e[FDR < 1e-5 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-1 & logFC > 0]) # 38994 - more sites
nrow(detable_e[FDR < 1e-1 & logFC < 0]) # 3272
detable_e[logFC < 0][1,]$FDR # 1.365118e-05 - the best FDR in the opposite direction logFC<0 could be a good threshold for our logFC>0

write.table(detable_e[FDR < 1e-10 & logFC > 0], "~/ap/Lib96_Leish_SMUG1-AP_S4.Lib100_Leish_SMUG1-AP_2_S7.clean.R1.rev.cigar_FDR1e-10.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

## bed
write.table(detable_e[FDR < 1e-10 & logFC > 0][, c("chr", "start", "end", "FDR")][start>-1], "~/ap/Lib96_Leish_SMUG1-AP_S4.Lib100_Leish_SMUG1-AP_2_S7.clean.R1.rev.cigar_FDR1e-10.bed", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

## volcano
df <- data.table(x = detable_e$logFC, y = -log10(detable_e$FDR), d = densCols(detable_e$logFC, -log10(detable_e$FDR), colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

gg <- ggplot(data = df, aes(x, y, col = d)) +
geom_point(size = 0.1) +
scale_color_identity() +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
geom_vline(xintercept = 0, linetype = 'dotted') +
theme_bw() +
coord_cartesian(xlim = c(-10, 10))

ggsave('~/figures/Lib96_Leish_SMUG1-AP.Lib100_Leish_SMUG1-AP_2.clean.R1.rev.cigar.png', width = 12, height = 12, units = 'cm')


########################
# snAP-seq vs. Y-input #
########################

# Load data
data <- fread("tableCat.py -i ~/ap/*Lib{97,101,99b,103}*rev*.cov -r .clean.R1.rev.cigar.cov")
setnames(data, c("chr", "pos", "count", "lib"))
table(data$lib)

# Cast
data_cast <- dcast.data.table(data = data, chr + pos ~ lib, value.var = "count")

# Define group
group <- factor(c('ap', 'input', 'ap', 'input'))

# Define DGEList object
y <- DGEList(counts = data_cast[,-c(1,2)], group = group, genes = data_cast[,c(1,2)])
y$samples

# Filter and get the top 100k sites according to cpm
y_f <- y[names(sort(rowMeans(cpm(y)), decreasing=T)[1:1e5]),]

# Define design matrix
des <- model.matrix(~ 0 + group, data = y_f$samples)
colnames(des) <- levels(factor(y_f$samples$group))

# Calculate normalization factors
y_f <- calcNormFactors(y_f, method = "TMM")

# Estimate dispersion
y_f <- estimateDisp(y_f, des)

# Fit linear model
fit <- glmFit(y_f, des, prior.count=1)

# Define contrasts
my.contrasts <- makeContrasts(apVSinput = ap - input, levels = des)

# Obtain likelihoods
lrt <- glmLRT(fit, contrast=my.contrasts[,"apVSinput"])

# Tables, bed files and volcano plots
## tables
detable <- topTags(lrt, n = Inf)$table
detable_e <- data.table(cbind(detable, data.frame(data_cast[,-c(1,2)][as.numeric(row.names(detable)),])))
detable_e[, start := as.integer(pos)]
detable_e[, end := as.integer(start + 1)]
detable_e <- detable_e[,.(chr, start, end, Lib97, Lib101, Lib99b, Lib103, logFC, logCPM, LR, PValue, FDR)]

nrow(detable_e[FDR < 1e-15 & logFC > 0]) # 0 - very high confidence sites
nrow(detable_e[FDR < 1e-15 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-10 & logFC > 0]) # 0 - high confidence sites
nrow(detable_e[FDR < 1e-10 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-5 & logFC > 0]) # 0 - sites of confidence
nrow(detable_e[FDR < 1e-5 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-1 & logFC > 0]) # 0 - more sites
nrow(detable_e[FDR < 1e-1 & logFC < 0]) # 0
detable_e[logFC < 0][1,]$FDR # 0.2735397

## volcano
df <- data.table(x = detable_e$logFC, y = -log10(detable_e$FDR), d = densCols(detable_e$logFC, -log10(detable_e$FDR), colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

gg <- ggplot(data = df, aes(x, y, col = d)) +
geom_point(size = 0.1) +
scale_color_identity() +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
geom_vline(xintercept = 0, linetype = 'dotted') +
theme_bw() +
coord_cartesian(xlim = c(-10, 10))

ggsave('~/figures/Lib97_Leish_AP.Lib101_Leish_AP_2.clean.R1.rev.cigar.png', width = 12, height = 12, units = 'cm')


############################
# UNG-snAP-seq vs. Y-input #
############################

# Load data
data <- fread("tableCat.py -i ~/ap/*Lib13{0..3}*rev*.cov -r .clean.R1.rev.cigar.cov")
setnames(data, c("chr", "pos", "count", "lib"))
table(data$lib)

# Cast
data_cast <- dcast.data.table(data = data, chr + pos ~ lib, value.var = "count")

# Define group
group <- factor(c('ap', 'ap', 'input', 'input'))

# Define DGEList object
y <- DGEList(counts = data_cast[,-c(1,2)], group = group, genes = data_cast[,c(1,2)])
y$samples

# Filter and get the top 100k sites according to cpm
y_f <- y[names(sort(rowMeans(cpm(y)), decreasing=T)[1:1e5]),]

# Define design matrix
des <- model.matrix(~ 0 + group, data = y_f$samples)
colnames(des) <- levels(factor(y_f$samples$group))

# Calculate normalization factors
y_f <- calcNormFactors(y_f, method = "TMM")

# Estimate dispersion
y_f <- estimateDisp(y_f, des)

# Fit linear model
fit <- glmFit(y_f, des, prior.count=1)

# Define contrasts
my.contrasts <- makeContrasts(apVSinput = ap - input, levels = des)

# Obtain likelihoods
lrt <- glmLRT(fit, contrast=my.contrasts[,"apVSinput"])

# Tables, bed files and volcano plots
## tables
detable <- topTags(lrt, n = Inf)$table
detable_e <- data.table(cbind(detable, data.frame(data_cast[,-c(1,2)][as.numeric(row.names(detable)),])))
detable_e[, start := as.integer(pos)]
detable_e[, end := as.integer(start + 1)]
detable_e <- detable_e[,.(chr, start, end, Lib130, Lib131, Lib132, Lib133, logFC, logCPM, LR, PValue, FDR)]

nrow(detable_e[FDR < 1e-15 & logFC > 0]) # 0 - very high confidence sites
nrow(detable_e[FDR < 1e-15 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-10 & logFC > 0]) # 0 - high confidence sites
nrow(detable_e[FDR < 1e-10 & logFC < 0]) # 0
nrow(detable_e[FDR < 1e-5 & logFC > 0]) # 0 - sites of confidence
nrow(detable_e[FDR < 1e-5 & logFC < 0]) # 4
nrow(detable_e[FDR < 1e-1 & logFC > 0]) # 15 - more sites
nrow(detable_e[FDR < 1e-1 & logFC < 0]) # 66
detable_e[logFC < 0][1,]$FDR # 3.998957e-09 - the best FDR in the opposite direction logFC<0 could be a good threshold for our logFC>0

## volcano
df <- data.table(x = detable_e$logFC, y = -log10(detable_e$FDR), d = densCols(detable_e$logFC, -log10(detable_e$FDR), colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

gg <- ggplot(data = df, aes(x, y, col = d)) +
geom_point(size = 0.1) +
scale_color_identity() +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
geom_vline(xintercept = 0, linetype = 'dotted') +
theme_bw() +
coord_cartesian(xlim = c(-10, 10))

ggsave('~/figures/Lib130_Leish_UNG-AP1.Lib131_Leish_UNG-AP2.clean.R1.rev.cigar.png', width = 12, height = 12, units = 'cm')
```



## Comparison of SMUG1-snAP-seq sites with peaks obtained in [Kawasaki et al., Genome biology, 2017](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1150-1)

Processed files:

- [fk_hmU.bed](files/fk_hmU.bed)
- [fk_hmU_baseJ.bed](files/fk_hmU_baseJ.bed)

Obtain chemseq 5hmU peaks in [GSE83384](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83384):

- fk050_Lmaj_chem1 [GSM2200481](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2200481)
- fk051_Lmaj_chem2 [GSM2200482](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2200482)
- fk054_Lmaj_chem3 [GSM2200484](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2200484)

Create bed files from the `_peaks.txt.gz` files

```bash
cd ~/bed

zcat GSM2200481_fk050_Lmaj_chem1_peaks.txt.gz | grep "^LmjF.*" | cut -f1-3 > GSM2200481_fk050_Lmaj_chem1_peaks.bed
zcat GSM2200482_fk051_Lmaj_chem2_peaks.txt.gz | grep "^LmjF.*" | cut -f1-3 > GSM2200482_fk051_Lmaj_chem2_peaks.bed
zcat GSM2200484_fk054_Lmaj_chem3_peaks.txt.gz | grep "^LmjF.*" | cut -f1-3 > GSM2200484_fk054_Lmaj_chem3_peaks.bed

wc -l *.bed
#206 GSM2200481_fk050_Lmaj_chem1_peaks.bed
#188 GSM2200482_fk051_Lmaj_chem2_peaks.bed
#145 GSM2200484_fk054_Lmaj_chem3_peaks.bed
```

Intersection:

```bash
cd ~/bed

for bed1 in GSM2200481*.bed fk_*.bed
do
  echo `wc -l $bed1`
  echo `wc -l ../ap/*Lib96*R1.fwd*FDR1e-10.bed`
  echo `wc -l ../ap/*Lib96*R1.rev*FDR1e-10.bed`
  bedtools intersect \
  -sorted \
  -a <(bedtools sort -i $bed1) \
  -b <(cat ../ap/*$id*FDR1e-10.bed | bedtools sort -i -) \
  -wa -u | wc -l
  bedtools intersect \
  -sorted \
  -a <(cat ../ap/*$id*FDR1e-10.bed | bedtools sort -i -) \
  -b <(bedtools sort -i $bed1) \
  -wa -u | wc -l
done
```



## Merge technical replicates and calculate coverage

```bash
cd ~/bam

mkdir ../bam_merge
mkdir ../bw_merge

##################
# SMUG1-snAP-seq #
##################

# fwd
samtools merge -@ 20 ../bam_merge/Leish_SMUG1-snAP-seq.clean.R1.fwd.cigar.bam Leish_SMUG1-AP_S4.clean.R1.fwd.cigar.bam LLeish_SMUG1-AP_2_S7.clean.R1.fwd.cigar.bam && \
samtools index ../bam_merge/Leish_SMUG1-snAP-seq.clean.R1.fwd.cigar.bam && \
bamCoverage -b ../bam_merge/Leish_SMUG1-snAP-seq.clean.R1.fwd.cigar.bam -o ../bw_merge/Leish_SMUG1-snAP-seq.clean.R1.fwd.cigar.bw -of bigwig --binSize 1 -p 20 --normalizeUsingRPKM

# rev
samtools merge -@ 20 ../bam_merge/Leish_SMUG1-snAP-seq.clean.R1.rev.cigar.bam Leish_SMUG1-AP_S4.clean.R1.rev.cigar.bam LLeish_SMUG1-AP_2_S7.clean.R1.rev.cigar.bam && \
samtools index ../bam_merge/Leish_SMUG1-snAP-seq.clean.R1.rev.cigar.bam && \
bamCoverage -b ../bam_merge/Leish_SMUG1-snAP-seq.clean.R1.rev.cigar.bam -o ../bw_merge/Leish_SMUG1-snAP-seq.clean.R1.rev.cigar.bw -of bigwig --binSize 1 -p 20 --normalizeUsingRPKM
```

Same for snAP-seq and input-Y libraries.



## Base composition profiles

### R1.fwd

#### Generate `.bam` and `.fasta` files

```bash
cd ~/ap

mkdir ../base_composition_profiles
mkdir ../base_composition_profiles/bam
mkdir ../base_composition_profiles/fasta
mkdir ../base_composition_profiles/tables

g=../reference/Lmajor.genome_spikeins.fasta

bed=Lib96_Leish_SMUG1-AP_S4.Lib100_Leish_SMUG1-AP_2_S7.clean.R1.fwd.cigar_FDR1e-10.bed

bam=${bed%.clean.R1.fwd.cigar_FDR1e-10.bed}.clean.R1.fwd.cigar.bam
bname=${bed%.bed}

cat $bed | cut -f1-3 | bedtools sort -i - | bedtools slop -i - -g $g.fai -b 10 | bedtools getfasta -fi $g -bed - > ../base_composition_profiles/fasta/$bname.fasta && \
echo "1" && \
samtools view ../bam_merge/$bam -L <(cat $bed | cut -f1-3 | bedtools sort -i - | bedtools slop -i - -g $g.fai -b 10) -b > ../base_composition_profiles/bam/$bname.inpeaks.bam && \
bedtools bamtobed -i ../base_composition_profiles/bam/$bname.inpeaks.bam | cut -f1-3 | bedtools sort -i - | bedtools flank -i - -g $g.fai -l 11 -r 0 | bedtools slop -i - -g $g.fai -l 0 -r 10 | awk '($3-$2) == 21' | bedtools getfasta -fi $g -bed - > ../base_composition_profiles/fasta/$bname.inpeaks.fasta && \
echo "2" && \
samtools view ../bam_merge/Leish_input_Y.clean.R1.fwd.cigar.bam -L <(cat $bed | cut -f1-3 | bedtools sort -i - | bedtools slop -i - -g $g.fai -b 10) -b > ../base_composition_profiles/bam/$bname.input_Y.bam && \
bedtools bamtobed -i ../base_composition_profiles/bam/$bname.input_Y.bam | cut -f1-3 | bedtools sort -i - | bedtools flank -i - -g $g.fai -l 11 -r 0 | bedtools slop -i - -g $g.fai -l 0 -r 10 | awk '($3-$2) == 21' | bedtools getfasta -fi $g -bed - > ../base_composition_profiles/fasta/$bname.input_Y.fasta && \
echo "3" && \
samtools view ../bam_merge/$bam -L <(cat $bed | cut -f1-3 | bedtools sort -i - | bedtools slop -i - -g $g.fai -b 10 | bedtools shuffle -seed 123 -noOverlapping -i - -g $g.fai) -b > ../base_composition_profiles/bam/$bname.shuffle.bam && \
bedtools bamtobed -i ../base_composition_profiles/bam/$bname.shuffle.bam | cut -f1-3 | bedtools sort -i - | bedtools flank -i - -g $g.fai -l 11 -r 0 | bedtools slop -i - -g $g.fai -l 0 -r 10 | awk '($3-$2) == 21' | bedtools getfasta -fi $g -bed - > ../base_composition_profiles/fasta/$bname.shuffle.fasta && \
echo "4" && \
bedtools bamtobed -i ../bam_merge/$bam | cut -f1-3 | bedtools sort -i - | bedtools flank -i - -g $g.fai -l 11 -r 0 | bedtools slop -i - -g $g.fai -l 0 -r 10 | awk '($3-$2) == 21' | bedtools getfasta -fi $g -bed - > ../base_composition_profiles/fasta/${bam%.bam}.fasta

# All reads from Leish_input_Y
bedtools bamtobed -i ../bam_merge/Leish_input_Y.clean.R1.fwd.cigar.bam | cut -f1-3 | bedtools sort -i - | bedtools flank -i - -g $g.fai -l 11 -r 0 | bedtools slop -i - -g $g.fai -l 0 -r 10 | awk '($3-$2) == 21' | bedtools getfasta -fi $g -bed - > ../base_composition_profiles/fasta/Leish_input_Y.clean.R1.fwd.cigar.fasta
```


#### Convert fasta file into table

Create tables from the `.fasta` files:

```python
import os

fasta_files = os.listdir("~/base_composition_profiles/fasta")

for f in fasta_files:
  print f
  ifasta = open("~/base_composition_profiles/fasta/%s" % f, "r")
  ilines = ifasta.read()
  ifasta.close()
  gcat_dict = {}
  for entry in ilines.split(">")[1:]:
    seq = list(entry.split()[1].upper())
    for i in range(0, len(seq)):
      base = seq[i]
      i_m = i - 10
      if base in "ACGT":
        if (i_m, base) in gcat_dict:
          gcat_dict[(i_m, base)]+=1
        else:
          gcat_dict[(i_m, base)]=1
  otable = open("~/base_composition_profiles/tables/%s" % f.replace(".fasta", "")+".txt", "w")
  otable.write("pos\tA\tC\tG\tT\n")
  for pos in range(-10, 11):
    otable.write("%s" % str(pos))
    for base in ["A", "C", "G", "T"]:
      if (pos, base) in gcat_dict:
        otable.write("\t%s" % str(gcat_dict[(pos, base)]))
      else:
        otable.write("\t0")
    otable.write("\n")
  otable.close()

```


#### Plotting profiles

```r
library(data.table)
library(ggplot2)

setwd("~/base_composition_profiles/tables")

for (f in list.files(".")){
  print(f)
  # Input data
  data <- fread(f)
  # Calculate percentages
  data <- cbind(data[, "pos"], 100*data[, c("A", "C", "G", "T")]/rowSums(data[, c("A", "C", "G", "T")]))
  # Melt data
  data_melt <- melt(data, id.vars = "pos", variable.name = "base", value.name = "pct")
  # Plot
  gg <- ggplot(data = data_melt, aes(x = pos, y = pct, colour = base)) +
  geom_line(size = 1) +
  #stat_smooth(aes(x = pos, y = pct), method = "lm", formula = y ~ poly(x, 10), se = FALSE) +
  #geom_point() +
  #geom_smooth(span = 0.1, se = FALSE) +
  theme_bw() +
  #theme_classic() +
  ylab("Percentage") +
  xlab("") +
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_continuous(breaks = seq(-10, 10), minor_breaks = NULL) +
  scale_y_continuous(breaks = seq(0, 100, 10), minor_breaks = NULL) +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=24), legend.text=element_text(size=20), legend.title=element_blank())
  # Save
  ggsave(gsub(".txt", ".v2.png", f), width = 20, units = "cm")
}
```


### R1.rev

#### Generate `.bam` and `.fasta` files

```bash
cd ~/ap

g=../reference/Lmajor.genome_spikeins.fasta

bed=Lib96_Leish_SMUG1-AP_S4.Lib100_Leish_SMUG1-AP_2_S7.clean.R1.rev.cigar_FDR1e-10.bed

bam=${bed%.clean.R1.rev.cigar_FDR1e-10.bed}.clean.R1.rev.cigar.bam
bname=${bed%.bed}

cat $bed | cut -f1-3 | bedtools sort -i - | bedtools slop -i - -g $g.fai -b 10 | bedtools getfasta -fi $g -bed - > ../base_composition_profiles/fasta/$bname.fasta && \
echo "1" && \
samtools view ../bam_merge/$bam -L <(cat $bed | cut -f1-3 | bedtools sort -i - | bedtools slop -i - -g $g.fai -b 10) -b > ../base_composition_profiles/bam/$bname.inpeaks.bam && \
bedtools bamtobed -i ../base_composition_profiles/bam/$bname.inpeaks.bam | cut -f1-3 | bedtools sort -i - | bedtools flank -i - -g $g.fai -l 0 -r 11 | bedtools slop -i - -g $g.fai -l 10 -r 0 | awk '($3-$2) == 21' | bedtools getfasta -fi $g -bed - > ../base_composition_profiles/fasta/$bname.inpeaks.fasta && \
echo "2" && \
samtools view ../bam_merge/Leish_input_Y.clean.R1.rev.cigar.bam -L <(cat $bed | cut -f1-3 | bedtools sort -i - | bedtools slop -i - -g $g.fai -b 10) -b > ../base_composition_profiles/bam/$bname.input_Y.bam && \
bedtools bamtobed -i ../base_composition_profiles/bam/$bname.input_Y.bam | cut -f1-3 | bedtools sort -i - | bedtools flank -i - -g $g.fai -l 0 -r 11 | bedtools slop -i - -g $g.fai -l 10 -r 0 | awk '($3-$2) == 21' | bedtools getfasta -fi $g -bed - > ../base_composition_profiles/fasta/$bname.input_Y.fasta && \
echo "3" && \
samtools view ../bam_merge/$bam -L <(cat $bed | cut -f1-3 | bedtools sort -i - | bedtools slop -i - -g $g.fai -b 10 | bedtools shuffle -seed 123 -noOverlapping -i - -g $g.fai) -b > ../base_composition_profiles/bam/$bname.shuffle.bam && \
bedtools bamtobed -i ../base_composition_profiles/bam/$bname.shuffle.bam | cut -f1-3 | bedtools sort -i - | bedtools flank -i - -g $g.fai -l 0 -r 11 | bedtools slop -i - -g $g.fai -l 10 -r 0 | awk '($3-$2) == 21' | bedtools getfasta -fi $g -bed - > ../base_composition_profiles/fasta/$bname.shuffle.fasta && \
echo "4" && \
bedtools bamtobed -i ../bam_merge/$bam | cut -f1-3 | bedtools sort -i - | bedtools flank -i - -g $g.fai -l 0 -r 11 | bedtools slop -i - -g $g.fai -l 10 -r 0 | awk '($3-$2) == 21' | bedtools getfasta -fi $g -bed - > ../base_composition_profiles/fasta/${bam%.bam}.fasta

# All reads from Leish_input_Y
bedtools bamtobed -i ../bam_merge/Leish_input_Y.clean.R1.rev.cigar.bam | cut -f1-3 | bedtools sort -i - | bedtools flank -i - -g $g.fai -l 0 -r 11 | bedtools slop -i - -g $g.fai -l 10 -r 0 | awk '($3-$2) == 21' | bedtools getfasta -fi $g -bed - > ../base_composition_profiles/fasta/Leish_input_Y.clean.R1.rev.cigar.fasta
```


#### Convert fasta file into table

Create tables from the `.fasta` files:

```python
import os

fasta_files = os.listdir("~/base_composition_profiles/fasta")

for f in fasta_files:
  if "R1.rev" in f:
    print f
    ifasta = open("~/base_composition_profiles/fasta/%s" % f, "r")
    ilines = ifasta.read()
    ifasta.close()
    gcat_dict = {}
    for entry in ilines.split(">")[1:]:
      seq = list(entry.split()[1].upper())
      for i in range(0, len(seq)):
        base = seq[i]
        i_m = i - 10
        if base in "ACGT":
          if (i_m, base) in gcat_dict:
            gcat_dict[(i_m, base)]+=1
          else:
            gcat_dict[(i_m, base)]=1
    otable = open("~/base_composition_profiles/tables/%s" % f.replace(".fasta", "")+".txt", "w")
    otable.write("pos\tA\tC\tG\tT\n")
    for pos in range(-10, 11):
      otable.write("%s" % str(pos))
      for base in ["A", "C", "G", "T"]:
        if (pos, base) in gcat_dict:
          otable.write("\t%s" % str(gcat_dict[(pos, base)]))
        else:
          otable.write("\t0")
      otable.write("\n")
    otable.close()

```


#### Plotting profiles

```r
library(data.table)
library(ggplot2)

setwd("~/base_composition_profiles/tables")

for (f in list.files(".")){
  if (grepl("R1.rev", f)){
    print(f)
    # Input data
    data <- fread(f)
    # Calculate percentages
    data <- cbind(data[, "pos"], 100*data[, c("A", "C", "G", "T")]/rowSums(data[, c("A", "C", "G", "T")]))
    # Melt data
    data_melt <- melt(data, id.vars = "pos", variable.name = "base", value.name = "pct")
    # Plot
    gg <- ggplot(data = data_melt, aes(x = pos, y = pct, colour = base)) +
    geom_line(size = 1) +
    #stat_smooth(aes(x = pos, y = pct), method = "lm", formula = y ~ poly(x, 10), se = FALSE) +
    #geom_point() +
    #geom_smooth(span = 0.1, se = FALSE) +
    theme_bw() +
    #theme_classic() +
    ylab("Percentage") +
    xlab("") +
    coord_cartesian(ylim = c(0, 100)) +
    scale_x_continuous(breaks = seq(-10, 10), minor_breaks = NULL) +
    scale_y_continuous(breaks = seq(0, 100, 10), minor_breaks = NULL) +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=24), legend.text=element_text(size=20), legend.title=element_blank())
    # Save
    ggsave(gsub(".txt", ".v2.png", f), width = 20, units = "cm")
  }
}
```



## Coverage profiles centered around the sites

computeMatrix and plotProfile:

```bash
cd ~/ap

mkdir ../deeptools

bed_SMUG1_snAP_seq_fwd='Lib96_Leish_SMUG1-AP_S4.Lib100_Leish_SMUG1-AP_2_S7.clean.R1.fwd.cigar_FDR1e-10.bed'
bed_SMUG1-snAP-seq_rev='Lib96_Leish_SMUG1-AP_S4.Lib100_Leish_SMUG1-AP_2_S7.clean.R1.rev.cigar_FDR1e-10.bed'

bw_SMUG1_snAP_seq_fwd="../bw_merge/Leish_SMUG1-snAP-seq.clean.R1.fwd.cigar.bw"
bw_SMUG1_snAP_seq_rev="../bw_merge/Leish_SMUG1-snAP-seq.clean.R1.rev.cigar.bw"
bw_snAP_seq_fwd="../bw_merge/Leish_snAP-seq.clean.R1.fwd.cigar.bw"
bw_snAP_seq_rev="../bw_merge/Leish_snAP-seq.clean.R1.rev.cigar.bw"
bw_input_Y_fwd="../bw_merge/Leish_input_Y.clean.R1.fwd.cigar.bw"
bw_input_Y_rev="../bw_merge/Leish_input_Y.clean.R1.rev.cigar.bw"


#######################
# reference-point fwd # 1000bp
#######################

## computeMatrix
computeMatrix reference-point \
-R $bed_SMUG1_snAP_seq_fwd \
-S $bw_SMUG1_snAP_seq_fwd $bw_snAP_seq_fwd $bw_input_Y_fwd \
-out ../deeptools/SMUG1_snAP_seq_fwd_referencepoint_bs1_rpkm_1000.mat.gz \
--referencePoint "center" \
-b 1000 \
-a 1000 \
-bs 1 \
--sortRegions no \
--skipZeros \
-p "max"

## plotProfile --plotType "se"
plotProfile \
--matrixFile ../deeptools/SMUG1_snAP_seq_fwd_referencepoint_bs1_rpkm_1000.mat.gz \
-out ../deeptools/SMUG1_snAP_seq_fwd_referencepoint_bs1_rpkm_1000_plotprofile_se.png \
--dpi 300 \
--plotHeight 9 \
--plotWidth 10 \
--plotType "se" \
--colors blue red green \
--refPointLabel "sites" \
--regionsLabel "" \
--samplesLabel "SMUG1-snAP" "snAP" "input" \
-y "RPKM" \
--perGroup


#######################
# reference-point rev # 1000bp
#######################

## computeMatrix
computeMatrix reference-point \
-R $bed_SMUG1_snAP_seq_rev \
-S $bw_SMUG1_snAP_seq_rev $bw_snAP_seq_rev $bw_input_Y_rev \
-out ../deeptools/SMUG1_snAP_seq_fwd_referencepoint_bs1_rpkm_1000.mat.gz \
--referencePoint "center" \
-b 1000 \
-a 1000 \
-bs 1 \
--sortRegions no \
--skipZeros \
-p "max"

## plotProfile --plotType "se"
plotProfile \
--matrixFile ../deeptools/SMUG1_snAP_seq_fwd_referencepoint_bs1_rpkm_1000.mat.gz \
-out ../deeptools/SMUG1_snAP_seq_fwd_referencepoint_bs1_rpkm_1000_plotprofile_se.png \
--dpi 300 \
--plotHeight 9 \
--plotWidth 10 \
--plotType "se" \
--colors blue red green \
--refPointLabel "sites" \
--regionsLabel "" \
--samplesLabel "SMUG1-snAP" "snAP" "input" \
-y "RPKM" \
--perGroup
```



## Sequence logos

### Add 5bp flanks and create fasta files

```bash
cd ~/ap

mkdir -p ../logos/fasta
mkdir -p ../logos/ggseqlogo

g=../reference/Lmajor.genome_spikeins.fasta

for bed in *.bed
do
  bname=${bed%.bed}
  cat $bed | cut -f1-3 | bedtools sort -i - | bedtools slop -i - -g $g.fai -b 5 | bedtools getfasta -fi $g -bed - > ../logos/fasta/$bname.fasta
done
```


### ggseqlogo

```r
# Load the required packages
require(ggplot2)
require(ggseqlogo)
require(Biostrings)

# Set working directory
setwd("~/logos/fasta")


##########
# R1.fwd #
##########

# Create list with all sequences at once
r1.fwd <- list()
for (f in list.files(pattern = "*R1.fwd*")){
  if (grepl("SMUG1-snAP-seq", f) & grepl("FDR1e-1\\.", f)){r1.fwd[[paste("SMUG1-snAP-seq FDR < 0.1 (n = ",length(as.data.frame(readDNAStringSet(f))$x),")", sep="")]] <- as.data.frame(readDNAStringSet(f))$x[which(sapply(as.data.frame(readDNAStringSet(f))$x, function(x) nchar(x)==11))]}
  if (grepl("SMUG1-snAP-seq", f) & grepl("FDR1e-5\\.", f)){r1.fwd[[paste("SMUG1-snAP-seq FDR < 10^(-5) (n = ",length(as.data.frame(readDNAStringSet(f))$x),")", sep="")]] <- as.data.frame(readDNAStringSet(f))$x[which(sapply(as.data.frame(readDNAStringSet(f))$x, function(x) nchar(x)==11))]}
  if (grepl("SMUG1-snAP-seq", f) & grepl("FDR1e-10\\.", f)){r1.fwd[[paste("SMUG1-snAP-seq FDR < 10^(-10) (n = ",length(as.data.frame(readDNAStringSet(f))$x),")", sep="")]] <- as.data.frame(readDNAStringSet(f))$x[which(sapply(as.data.frame(readDNAStringSet(f))$x, function(x) nchar(x)==11))]}
  if (grepl("SMUG1-snAP-seq", f) & grepl("FDR1e-15\\.", f)){r1.fwd[[paste("SMUG1-snAP-seq FDR < 10^(-15) (n = ",length(as.data.frame(readDNAStringSet(f))$x),")", sep="")]] <- as.data.frame(readDNAStringSet(f))$x[which(sapply(as.data.frame(readDNAStringSet(f))$x, function(x) nchar(x)==11))]}
}

# Create custom colour scheme
cs <- make_col_scheme(chars=c('A', 'C', 'G', 'T'), cols=hue_pal()(4))

# Plot sequence logo bits and prob
gg <- ggplot() + geom_logo(r1.fwd[c(1,4,2,3,5,8,6,7)], method = 'bits', col_scheme = cs) + theme_logo() + facet_wrap(~seq_group, ncol = 4, scales='free_x') + theme(axis.line.y = element_line(color="black"))
ggsave('../figures/SMUG1-snAP-seq.R1.fwd_bits.png', width = 36, height = 12, units = 'cm')

gg <- ggplot() + geom_logo(r1.fwd[c(1,4,2,3,5,8,6,7)], method = 'prob', col_scheme = cs) + theme_logo() + facet_wrap(~seq_group, ncol = 4, scales='free_x') + theme(axis.line.y = element_line(color="black"))
ggsave('../figures/SMUG1-snAP-seq.R1.fwd_prob.png', width = 36, height = 12, units = 'cm')

# e.g. just for SMUG1-snAP-seq FDR < 10^(-10) (n = 1560)
ggseqlogo(r1.fwd[2], method = 'bits', col_scheme = cs)
ggplot() + geom_logo(r1.fwd[2], method = 'bits', col_scheme = cs) + theme_logo()
ggplot() + geom_logo(r1.fwd[2], method = 'bits', col_scheme = cs) + theme_logo() + theme(axis.line.y = element_line(color="black"))
ggplot() + geom_logo(r1.fwd[2], method = 'bits', col_scheme = cs) + theme_logo() + theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))


##########
# R1.rev #
##########

# Create list with all sequences at once
r1.rev <- list()
for (f in list.files(pattern = "*R1.rev*")){
  if (grepl("Lib96", f) & grepl("FDR1e-1\\.", f)){r1.rev[[paste("Lib96 and Lib100 FDR < 0.1 (n = ",length(as.data.frame(readDNAStringSet(f))$x),")", sep="")]] <- as.data.frame(readDNAStringSet(f))$x[which(sapply(as.data.frame(readDNAStringSet(f))$x, function(x) nchar(x)==11))]}
  if (grepl("Lib96", f) & grepl("FDR1e-5\\.", f)){r1.rev[[paste("Lib96 and Lib100 FDR < 10^(-5) (n = ",length(as.data.frame(readDNAStringSet(f))$x),")", sep="")]] <- as.data.frame(readDNAStringSet(f))$x[which(sapply(as.data.frame(readDNAStringSet(f))$x, function(x) nchar(x)==11))]}
  if (grepl("Lib96", f) & grepl("FDR1e-10\\.", f)){r1.rev[[paste("Lib96 and Lib100 FDR < 10^(-10) (n = ",length(as.data.frame(readDNAStringSet(f))$x),")", sep="")]] <- as.data.frame(readDNAStringSet(f))$x[which(sapply(as.data.frame(readDNAStringSet(f))$x, function(x) nchar(x)==11))]}
  if (grepl("Lib96", f) & grepl("FDR1e-15\\.", f)){r1.rev[[paste("Lib96 and Lib100 FDR < 10^(-15) (n = ",length(as.data.frame(readDNAStringSet(f))$x),")", sep="")]] <- as.data.frame(readDNAStringSet(f))$x[which(sapply(as.data.frame(readDNAStringSet(f))$x, function(x) nchar(x)==11))]}
}

# Create custom colour scheme
cs <- make_col_scheme(chars=c('A', 'C', 'G', 'T'), cols=hue_pal()(4))

# Plot sequence logo bits and prob
gg <- ggplot() + geom_logo(r1.rev[c(1,4,2,3,5,8,6,7)], method = 'bits', col_scheme = cs) + theme_logo() + facet_wrap(~seq_group, ncol = 4, scales='free_x') + theme(axis.line.y = element_line(color="black"))
ggsave('../figures/SMUG1-snAP-seq.R1.rev_bits.png', width = 36, height = 12, units = 'cm')

gg <- ggplot() + geom_logo(r1.rev[c(1,4,2,3,5,8,6,7)], method = 'prob', col_scheme = cs) + theme_logo() + facet_wrap(~seq_group, ncol = 4, scales='free_x') + theme(axis.line.y = element_line(color="black"))
ggsave('../figures/SMUG1-snAP-seq.R1.rev_prob.png', width = 36, height = 12, units = 'cm')
```



## Exploring TG enrichment using frequencies of dinucleotides counts, motifs and tests

### frequencies

- Calculate dinucleotide frequencies in the entire Leishmania genome:
- Calculate dinucleotide frequencies in Fumi's 5hmU peaks using the Leishmania genome as expected frequencies

```bash
cd ~/reference/

compseq Lmajor.genome.fasta -word 2 Lmajor.genome.fasta.comp

#
# Output from 'compseq'
#
# The Expected frequencies are calculated on the (false) assumption that every
# word has equal frequency.
#
# The input sequences are:
#	LmjF.01
#	LmjF.02
#	LmjF.03
#	LmjF.04
#	LmjF.05
#	LmjF.06
#	LmjF.07
#	LmjF.08
#	LmjF.09
#	LmjF.10
# ... et al.
#
#Word size	2
#Total count	32855059
#
#
# Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
#
#AA	1380153		0.0420073	0.0625000	0.6721171
#AC	1990373		0.0605804	0.0625000	0.9692866
#AG	2073185		0.0631009	0.0625000	1.0096150
#AT	1116719		0.0339893	0.0625000	0.5438281
#CA	2387571		0.0726698	0.0625000	1.1627170
#CC	2363431		0.0719351	0.0625000	1.1509611
#CG	2974892		0.0905459	0.0625000	1.4487349
#CT	2126964		0.0647378	0.0625000	1.0358047
#GA	2053568		0.0625039	0.0625000	1.0000618
#GC	3390395		0.1031925	0.0625000	1.6510797
#GG	2318901		0.0705797	0.0625000	1.1292756
#GT	2005230		0.0610326	0.0625000	0.9765218
#TA	739139		0.0224970	0.0625000	0.3599514
#TC	2108650		0.0641804	0.0625000	1.0268860
#TG	2401122		0.0730823	0.0625000	1.1693162
#TT	1424746		0.0433646	0.0625000	0.6938334
#
#Other	20		0.0000006	0.0000000	10000000000.0000000


cd ~/bed
bedtools getfasta -fi ../reference/Lmajor.genome.fasta -bed fk_hmU.bed > fk_hmU.fasta
samtools faidx fk_hmU.fasta
compseq fk_hmU.fasta -word 2 fk_hmU.fasta.comp -in ~/reference/Lmajor.genome.fasta.comp

#
# Output from 'compseq'
#
# The Expected frequencies are taken from the file: ~/reference/Lmajor.genome.fasta.comp
#
# The input sequences are:
#	256134-259249
#	260820-267767
#	268567-268988
#	263559-269641
#	246789-253447
#	258512-259983
#	383220-384273
#	141-415
#	701-1604
#	125200-132538
# ... et al.
#
#
#Word size	2
#Total count	260140
#
#
# Word	Obs Count	Obs Frequency	Exp Frequency	Obs/Exp Frequency
#
#AA	18324		0.0704390	0.0420073	1.6768275
#AC	13471		0.0517837	0.0605804	0.8547922
#AG	13922		0.0535173	0.0631009	0.8481232
#AT	8277		0.0318175	0.0339893	0.9361029
#CA	16822		0.0646652	0.0726698	0.8898494
#CC	24727		0.0950527	0.0719351	1.3213670
#CG	20519		0.0788768	0.0905459	0.8711246
#CT	13651		0.0524756	0.0647378	0.8105866
#GA	12819		0.0492773	0.0625039	0.7883878
#GC	24967		0.0959752	0.1031925	0.9300603
#GG	25643		0.0985738	0.0705797	1.3966317
#GT	13415		0.0515684	0.0610326	0.8449318
#TA	6010		0.0231029	0.0224970	1.0269345
#TC	12567		0.0483086	0.0641804	0.7527002
#TG	16762		0.0644345	0.0730823	0.8816709
#TT	18242		0.0701238	0.0433646	1.6170743
#
#Other	2		0.0000077	0.0000006	12.8136132
```


### dreme

1. Create fasta file containing +/- 2bp from FDR 10^(-10) R1.fwd SMUG1-snAP-seq sites

2. Create negative files:
  - +/- 2bp from all Ts genome-wide
  - +/- 2bp from all Ts within 5hmU peaks

```bash
cd ~
mkdir tg
cd tg

mkdir -p dreme/bed
mkdir -p dreme/fasta
mkdir -p dreme/run

g=~/reference//Lmajor.genome.fasta

cat ../ap/Lib96_Leish_SMUG1-AP_S4.Lib100_Leish_SMUG1-AP_2_S7.clean.R1.fwd.cigar_FDR1e-10.bed | \
cut -f1-3 | \
bedtools sort -i - | \
bedtools slop -i - -g $g.fai -b 2 | \
bedtools getfasta -fi $g -bed - > dreme/fasta/Lib96_Leish_SMUG1-AP_S4.Lib100_Leish_SMUG1-AP_2_S7.clean.R1.fwd.cigar_FDR1e-10.fasta

fastaRegexFinder.py -f $g -r T --noreverse | cut -f1-3 | bedtools slop -i - -g $g.fai -b 2 | awk -v OFS="\t" '$3 - $2 == 5' | bedtools getfasta -fi $g -bed - > dreme/fasta/Lmajor.genome.T.2bp.fasta
grep ">" dreme/fasta/Lmajor.genome.T.2bp.fasta | wc -l # 6673621

ref=~/bed/fk_hmU.fasta
fastaRegexFinder.py -f $ref -r T --noreverse | cut -f1-3 | bedtools slop -i - -g $ref.fai -b 2 | awk -v OFS="\t" '$3 - $2 == 5' | bedtools getfasta -fi $ref -bed - > dreme/fasta/fk_hmU.T.2bp.fasta
grep ">" dreme/fasta/fk_hmU.T.2bp.fasta | wc -l # 53505
```

3. Run dreme with -k 1-5 against negative files from second step

```bash
cd ~/tg

p=dreme/fasta/Lib96_Leish_SMUG1-AP_S4.Lib100_Leish_SMUG1-AP_2_S7.clean.R1.fwd.cigar_FDR1e-10.fasta
n=dreme/fasta/Lmajor.genome.T.2bp.fasta

for k in {1..5}
do
  dreme -k $k -norc -oc dreme/run/Lmajor.genome.T.2bp_k$k -p $p -n $n -png
done

p=dreme/fasta/Lib96_Leish_SMUG1-AP_S4.Lib100_Leish_SMUG1-AP_2_S7.clean.R1.fwd.cigar_FDR1e-10.fasta
n=dreme/fasta/fk_hmU.T.2bp.fasta

for k in {1..5}
do
  dreme -k $k -norc -oc dreme/run/fk_hmU.T.2bp_k$k -p $p -n $n -png
done
```


### Test for TG enrichment in SMUG1-snAP-seq

Following discussions with Jane:

- Calculate TG and TX frequencies in Lib96+Lib100 R1.fwd sites at FDR 10^(-10) (1560)
- Generate samples (e.g. 10000) of size 1560 from Ts in Fumi's hmU regions and Ts genome-wide and calculate TG and TX frequencies in each where X = A, C, T
- Generate samples (e.g. 10000) of size 1560 from synthetic 5-hmU N-oligo libraries and calculate TG and TX frequencies in each

```python
import collections
import random
from scipy import stats
import numpy
import gzip
import os
import re
import pandas as pd


## Calculate TG and TX frequencies in SMUG1-snAP-seq R1.fwd sites at FDR 10^(-10) (1560)
lib96_lib100_r1_fwd_ifile = open("~/tg/dreme/fasta/Lib96_Leish_SMUG1-AP_S4.Lib100_Leish_SMUG1-AP_2_S7.clean.R1.fwd.cigar_FDR1e-10.fasta", "r")
lib96_lib100_r1_fwd_fasta = lib96_lib100_r1_fwd_ifile.read().split(">")[1:]
lib96_lib100_r1_fwd_ifile.close()

len(lib96_lib100_r1_fwd_fasta) # 1560, good

lib96_lib100_r1_fwd_cnt = collections.Counter()

for i in lib96_lib100_r1_fwd_fasta:
  dn = i.split()[1][2:4]
  lib96_lib100_r1_fwd_cnt[dn] += 1

lib96_lib100_r1_fwd_cnt['tg'] # 1044
lib96_lib100_r1_fwd_cnt['ta'] + lib96_lib100_r1_fwd_cnt['tc'] + lib96_lib100_r1_fwd_cnt['tt'] # 483
lib96_lib100_r1_fwd_cnt['ta'] + lib96_lib100_r1_fwd_cnt['tc'] + lib96_lib100_r1_fwd_cnt['tg'] + lib96_lib100_r1_fwd_cnt['tt'] # 1527
100*float(lib96_lib100_r1_fwd_cnt['tg'])/(lib96_lib100_r1_fwd_cnt['ta'] + lib96_lib100_r1_fwd_cnt['tc'] + lib96_lib100_r1_fwd_cnt['tg'] + lib96_lib100_r1_fwd_cnt['tt']) # 68% (68.36935166994107)
100*float(lib96_lib100_r1_fwd_cnt['ta'] + lib96_lib100_r1_fwd_cnt['tc'] + lib96_lib100_r1_fwd_cnt['tt'])/(lib96_lib100_r1_fwd_cnt['ta'] + lib96_lib100_r1_fwd_cnt['tc'] + lib96_lib100_r1_fwd_cnt['tg'] + lib96_lib100_r1_fwd_cnt['tt']) # 32% (31.63064833005894)



## Generate samples of size 1560 from Ts in hmU regions and calculate TG and TX frequencies in each
fk_hmU_ifile = open("~/tg/dreme/fasta/fk_hmU.T.2bp.fasta", "r")
fk_hmU_fasta = fk_hmU_ifile.read().split(">")[1:]
fk_hmU_ifile.close()

len(fk_hmU_fasta) # 53505, good

fk_hmU_cnt = collections.Counter()

for i in fk_hmU_fasta:
  dn = i.split()[1][2:4]
  fk_hmU_cnt[dn] += 1

fk_hmU_cnt['tg'] # 16743
fk_hmU_cnt['ta'] + fk_hmU_cnt['tc'] + fk_hmU_cnt['tt'] # 36761
fk_hmU_cnt['ta'] + fk_hmU_cnt['tc'] + fk_hmU_cnt['tg'] + fk_hmU_cnt['tt'] # 53504
100*float(fk_hmU_cnt['tg'])/(fk_hmU_cnt['ta'] + fk_hmU_cnt['tc'] + fk_hmU_cnt['tg'] + fk_hmU_cnt['tt']) # 31% (31.29298744019139)
100*float(fk_hmU_cnt['ta'] + fk_hmU_cnt['tc'] + fk_hmU_cnt['tt'])/(fk_hmU_cnt['ta'] + fk_hmU_cnt['tc'] + fk_hmU_cnt['tg'] + fk_hmU_cnt['tt']) # 69% (68.70701255980862)

fk_hmU_s = []

for i in range(1, 10001):
  random.seed(i)
  print(i)
  fk_hmU_fasta_s = random.sample(fk_hmU_fasta, 1560)
  fk_hmU_cnt_s = collections.Counter()
  for i in fk_hmU_fasta_s:
    dn = i.split()[1][2:4]
    fk_hmU_cnt_s[dn] += 1
  tg_pct = 100*float(fk_hmU_cnt_s['tg'])/(fk_hmU_cnt_s['ta'] + fk_hmU_cnt_s['tc'] + fk_hmU_cnt_s['tg'] + fk_hmU_cnt_s['tt'])
  tx_pct = 100*float(fk_hmU_cnt_s['ta'] + fk_hmU_cnt_s['tc'] + fk_hmU_cnt_s['tt'])/(fk_hmU_cnt_s['ta'] + fk_hmU_cnt_s['tc'] + fk_hmU_cnt_s['tg'] + fk_hmU_cnt_s['tt'])
  fk_hmU_s.append((tg_pct, tx_pct))


stats.describe([s[0] for s in fk_hmU_s]) # DescribeResult(nobs=10000, minmax=(27.115384615384617, 35.96153846153846), mean=31.283412596009939, variance=1.3495955748538351, skewness=0.04975492476248911, kurtosis=-0.05909367536955523)
stats.describe([s[1] for s in fk_hmU_s]) # DescribeResult(nobs=10000, minmax=(64.038461538461533, 72.884615384615387), mean=68.716587403990061, variance=1.3495955748538355, skewness=-0.04975492476248899, kurtosis=-0.059093675369556564)
stats.ttest_1samp([s[0] for s in fk_hmU_s], 100*float(lib96_lib100_r1_fwd_cnt['tg'])/(lib96_lib100_r1_fwd_cnt['ta'] + lib96_lib100_r1_fwd_cnt['tc'] + lib96_lib100_r1_fwd_cnt['tg'] + lib96_lib100_r1_fwd_cnt['tt'])) # pvalue=0.0



## Generate samples of size 1560 from Ts genome-wide and calculate TG and TX frequencies in each
lmajor_genome_ifile = open("~/tg/dreme/fasta/Lmajor.genome.T.2bp.fasta", "r")
lmajor_genome_fasta = lmajor_genome_ifile.read().split(">")[1:]
lmajor_genome_ifile.close()

len(lmajor_genome_fasta) # 6673621, good

lmajor_genome_cnt = collections.Counter()

for i in lmajor_genome_fasta:
  dn = i.split()[1][2:4]
  lmajor_genome_cnt[dn] += 1

lmajor_genome_cnt['tg'] # 2401121
lmajor_genome_cnt['ta'] + lmajor_genome_cnt['tc'] + lmajor_genome_cnt['tt'] # 4272497
lmajor_genome_cnt['ta'] + lmajor_genome_cnt['tc'] + lmajor_genome_cnt['tg'] + lmajor_genome_cnt['tt'] # 6673618
100*float(lmajor_genome_cnt['tg'])/(lmajor_genome_cnt['ta'] + lmajor_genome_cnt['tc'] + lmajor_genome_cnt['tg'] + lmajor_genome_cnt['tt']) # 36% (35.97929938453175)
100*float(lmajor_genome_cnt['ta'] + lmajor_genome_cnt['tc'] + lmajor_genome_cnt['tt'])/(lmajor_genome_cnt['ta'] + lmajor_genome_cnt['tc'] + lmajor_genome_cnt['tg'] + lmajor_genome_cnt['tt']) # 64% (64.02070061546826)

lmajor_genome_s = []

for i in range(1, 10001):
  random.seed(i)
  print(i)
  lmajor_genome_fasta_s = random.sample(lmajor_genome_fasta, 1560)
  lmajor_genome_cnt_s = collections.Counter()
  for i in lmajor_genome_fasta_s:
    dn = i.split()[1][2:4]
    lmajor_genome_cnt_s[dn] += 1
  tg_pct = 100*float(lmajor_genome_cnt_s['tg'])/(lmajor_genome_cnt_s['ta'] + lmajor_genome_cnt_s['tc'] + lmajor_genome_cnt_s['tg'] + lmajor_genome_cnt_s['tt'])
  tx_pct = 100*float(lmajor_genome_cnt_s['ta'] + lmajor_genome_cnt_s['tc'] + lmajor_genome_cnt_s['tt'])/(lmajor_genome_cnt_s['ta'] + lmajor_genome_cnt_s['tc'] + lmajor_genome_cnt_s['tg'] + lmajor_genome_cnt_s['tt'])
  lmajor_genome_s.append((tg_pct, tx_pct))


stats.describe([s[0] for s in lmajor_genome_s]) # DescribeResult(nobs=10000, minmax=(30.192307692307693, 40.512820512820511), mean=35.98778068206115, variance=1.4912509477666627, skewness=0.018695718913297632, kurtosis=-0.01373407151315753)
stats.describe([s[1] for s in lmajor_genome_s]) # DescribeResult(nobs=10000, minmax=(59.487179487179489, 69.807692307692307), mean=64.012219317938872, variance=1.4912509477666627, skewness=-0.01869571891335036, kurtosis=-0.013734071513155754)
stats.ttest_1samp([s[0] for s in lmajor_genome_s], 100*float(lib96_lib100_r1_fwd_cnt['tg'])/(lib96_lib100_r1_fwd_cnt['ta'] + lib96_lib100_r1_fwd_cnt['tc'] + lib96_lib100_r1_fwd_cnt['tg'] + lib96_lib100_r1_fwd_cnt['tt'])) # pvalue=0.0



## Generate samples of size 1560 from Lib138 and calculate TG and TX frequencies in each

# Define variables
d = "~/fastq_trimmed"
files = os.listdir(d)
lib = "Lib138"

# Load reads into memory
reads = "".join([gzip.open(d + "/" + f, 'rt').read() for f in files if (lib in f) and ("_R1_001.fastq.gz" in f)]).split("@")[1:]
len(reads) # 594549, ok

# Create dictionary rname : right_barcode
fwd_seq = "GTAGTAGTCGACTAG"

rn_barcode = {}

for r in reads:
  fields = r.split('\n')
  rn = fields[0]
  s = fields[1]
  if fwd_seq in s:
    idx = re.search(fwd_seq, s).start()
    barcode = s[idx-10:idx]
    if len(barcode) == 10 and "N" not in barcode:
      rn_barcode[rn] = barcode


len(rn_barcode.values()) # 407781 barcodes before deduplication
len(set(rn_barcode.values())) # 217029 (53%) unique barcodes when deduplicating

# Deduplicate dictionary
barcode_rn = {v: k for k, v in rn_barcode.iteritems()}
rn_barcode_dedup = {v: k for k, v in barcode_rn.iteritems()}

# Counting like above
lib138_cnt = collections.Counter()

for i in rn_barcode_dedup.values():
  n = i[0]
  lib138_cnt[n] += 1

lib138_cnt['G'] # 86907
lib138_cnt['A'] + lib138_cnt['C'] + lib138_cnt['T'] # 130122
lib138_cnt['A'] + lib138_cnt['C'] + lib138_cnt['G'] + lib138_cnt['T'] # 217029
100*float(lib138_cnt['G'])/(lib138_cnt['A'] + lib138_cnt['C'] + lib138_cnt['G'] + lib138_cnt['T']) # 40% (40.04395725916813)
100*float(lib138_cnt['A'] + lib138_cnt['C'] + lib138_cnt['T'])/(lib138_cnt['A'] + lib138_cnt['C'] + lib138_cnt['G'] + lib138_cnt['T']) # 60% (59.95604274083187)

lib138_s = []

for i in range(1, 10001):
  random.seed(i)
  print(i)
  rn_barcode_dedup_s = random.sample(rn_barcode_dedup.values(), 1560)
  lib138_cnt_s = collections.Counter()
  for i in rn_barcode_dedup_s:
    n = i[0]
    lib138_cnt_s[n] += 1
  tg_pct = 100*float(lib138_cnt_s['G'])/(lib138_cnt_s['A'] + lib138_cnt_s['C'] + lib138_cnt_s['G'] + lib138_cnt_s['T'])
  tx_pct = 100*float(lib138_cnt_s['A'] + lib138_cnt_s['C'] + lib138_cnt_s['T'])/(lib138_cnt_s['A'] + lib138_cnt_s['C'] + lib138_cnt_s['G'] + lib138_cnt_s['T'])
  lib138_s.append((tg_pct, tx_pct))


stats.describe([s[0] for s in lib138_s]) # DescribeResult(nobs=10000, minmax=(35.064102564102562, 44.807692307692307), mean=40.057192307692304, variance=1.5296910162745398, skewness=0.00239954716476055, kurtosis=-0.052475425671872244)
stats.describe([s[1] for s in lib138_s]) # DescribeResult(nobs=10000, minmax=(55.192307692307693, 64.935897435897431), mean=59.942807692307703, variance=1.5296910162745398, skewness=-0.002399547164777818, kurtosis=-0.05247542567187313)



## Generate samples of size 1560 from Lib139 and calculate TG and TX frequencies in each

# Define variables
d = "~/fastq_trimmed"
files = os.listdir(d)
lib = "Lib139"

# Load reads into memory
reads = "".join([gzip.open(d + "/" + f, 'rt').read() for f in files if (lib in f) and ("_R2_001.fastq.gz" in f)]).split("@")[1:]
len(reads) # 594549, ok

# Create dictionaries N10 rname : right_barcode and N21 rname : left_barcode + modified_base + right_barcode
fwd_seq = "GTAGTAGTCGACTAG"

rn_barcode_N10 = {}
rn_barcode_N21_1 = {}
rn_barcode_N21_2 = {}

for r in reads:
  fields = r.split('\n')
  rn = fields[0]
  s = fields[1]
  if fwd_seq in s:
    idx = re.search(fwd_seq, s).start()
    barcode_N10 = s[idx-10:idx]
    barcode_N21 = s[idx-21:idx]
    if len(barcode_N10) == 10 and "N" not in barcode_N10:
      rn_barcode_N10[rn] = barcode_N10
      rn_barcode_N21_1[rn] = barcode_N21
    if len(barcode_N21) == 21 and "N" not in barcode_N21:
      rn_barcode_N21_2[rn] = barcode_N21


len(rn_barcode_N10.values()) # 574272 barcodes before deduplication
len(set(rn_barcode_N10.values())) # 281120 (49%) unique barcodes when deduplicating

len(rn_barcode_N21_1.values()) # 574272 barcodes before deduplication
len(set(rn_barcode_N21_1.values())) # 552203 (96%) unique barcodes when deduplicating

len(rn_barcode_N21_2.values()) # 563264 barcodes before deduplication
len(set(rn_barcode_N21_2.values())) # 552169 (98%) unique barcodes when deduplicating

# Deduplicate N10 dictionary
barcode_rn_N10 = {v: k for k, v in rn_barcode_N10.iteritems()}
rn_barcode_N10_dedup = {v: k for k, v in barcode_rn_N10.iteritems()}
len(set(rn_barcode_N10_dedup.values())) # 281120

# Reduce N21 dictionary based on the deduplication of N10 dictionary
rn_barcode_N21_1_dedup = {rn: rn_barcode_N21_1[rn] for rn in rn_barcode_N10_dedup.keys() if len(rn_barcode_N21_1[rn]) == 21 and "N" not in rn_barcode_N21_1[rn]}
len(rn_barcode_N21_1_dedup) # 276402
len(set(rn_barcode_N21_1_dedup.values())) # 276402

# Deduplicate N21 dictionary
barcode_rn_N21_2 = {v: k for k, v in rn_barcode_N21_2.iteritems()}
rn_barcode_N21_2_dedup = {v: k for k, v in barcode_rn_N21_2.iteritems()}
len(set(rn_barcode_N21_2_dedup.values())) # 552169

# Counting like above
lib139_cnt = collections.Counter()

for i in rn_barcode_N10_dedup.values():
  n = i[0]
  lib139_cnt[n] += 1

lib139_cnt['G'] # 101968
lib139_cnt['A'] + lib139_cnt['C'] + lib139_cnt['T'] # 179152
lib139_cnt['A'] + lib139_cnt['C'] + lib139_cnt['G'] + lib139_cnt['T'] # 281120
100*float(lib139_cnt['G'])/(lib139_cnt['A'] + lib139_cnt['C'] + lib139_cnt['G'] + lib139_cnt['T']) # 36% (36.2720546385885)
100*float(lib139_cnt['A'] + lib139_cnt['C'] + lib139_cnt['T'])/(lib139_cnt['A'] + lib139_cnt['C'] + lib139_cnt['G'] + lib139_cnt['T']) # 64% (63.7279453614115)

lib139_s = []

for i in range(1, 10001):
  random.seed(i)
  print(i)
  rn_barcode_N10_dedup_s = random.sample(rn_barcode_N10_dedup.values(), 1560)
  lib139_cnt_s = collections.Counter()
  for i in rn_barcode_N10_dedup_s:
    n = i[0]
    lib139_cnt_s[n] += 1
  tg_pct = 100*float(lib139_cnt_s['G'])/(lib139_cnt_s['A'] + lib139_cnt_s['C'] + lib139_cnt_s['G'] + lib139_cnt_s['T'])
  tx_pct = 100*float(lib139_cnt_s['A'] + lib139_cnt_s['C'] + lib139_cnt_s['T'])/(lib139_cnt_s['A'] + lib139_cnt_s['C'] + lib139_cnt_s['G'] + lib139_cnt_s['T'])
  lib139_s.append((tg_pct, tx_pct))


stats.describe([s[0] for s in lib139_s]) # DescribeResult(nobs=10000, minmax=(31.153846153846153, 41.53846153846154), mean=36.284782051282058, variance=1.4653650820374611, skewness=-0.010417139801069421, kurtosis=-0.0038983912495167417)
stats.describe([s[1] for s in lib139_s]) # DescribeResult(nobs=10000, minmax=(58.46153846153846, 68.84615384615384), mean=63.71521794871795, variance=1.4653650820374606, skewness=0.010417139801051465, kurtosis=-0.0038983912495158535)
stats.ttest_1samp([s[0] for s in lib139_s], stats.describe([s[0] for s in lib138_s]).mean) # pvalue=0.0
stats.ttest_ind([s[0] for s in lib139_s], [s[0] for s in lib138_s]) # pvalue=0.0



## Generate samples of size 1560 from Lib140 and calculate TG and TX frequencies in each

# Define variables
d = "~/fastq_trimmed"
files = os.listdir(d)
lib = "Lib140"

# Load reads into memory
reads = "".join([gzip.open(d + "/" + f, 'rt').read() for f in files if (lib in f) and ("_R2_001.fastq.gz" in f)]).split("@")[1:]
len(reads) # 839914, ok

# Create dictionaries N10 rname : right_barcode and N21 rname : left_barcode + modified_base + right_barcode
fwd_seq = "GTAGTAGTCGACTAG"

rn_barcode_N10 = {}
rn_barcode_N21_1 = {}
rn_barcode_N21_2 = {}

for r in reads:
  fields = r.split('\n')
  rn = fields[0]
  s = fields[1]
  if fwd_seq in s:
    idx = re.search(fwd_seq, s).start()
    barcode_N10 = s[idx-10:idx]
    barcode_N21 = s[idx-21:idx]
    if len(barcode_N10) == 10 and "N" not in barcode_N10:
      rn_barcode_N10[rn] = barcode_N10
      rn_barcode_N21_1[rn] = barcode_N21
    if len(barcode_N21) == 21 and "N" not in barcode_N21:
      rn_barcode_N21_2[rn] = barcode_N21


len(rn_barcode_N10.values()) # 557467 barcodes before deduplication
len(set(rn_barcode_N10.values())) # 276071 (50%) unique barcodes when deduplicating

len(rn_barcode_N21_1.values()) # 557467 barcodes before deduplication
len(set(rn_barcode_N21_1.values())) # 500449 (89%) unique barcodes when deduplicating

len(rn_barcode_N21_2.values()) # 507145 barcodes before deduplication
len(set(rn_barcode_N21_2.values())) # 500414 (99%) unique barcodes when deduplicating

# Deduplicate N10 dictionary
barcode_rn_N10 = {v: k for k, v in rn_barcode_N10.iteritems()}
rn_barcode_N10_dedup = {v: k for k, v in barcode_rn_N10.iteritems()}
len(set(rn_barcode_N10_dedup.values())) # 276071

# Reduce N21 dictionary based on the deduplication of N10 dictionary
rn_barcode_N21_1_dedup = {rn: rn_barcode_N21_1[rn] for rn in rn_barcode_N10_dedup.keys() if len(rn_barcode_N21_1[rn]) == 21 and "N" not in rn_barcode_N21_1[rn]}
len(rn_barcode_N21_1_dedup) # 252457
len(set(rn_barcode_N21_1_dedup.values())) # 252457

# Deduplicate N21 dictionary
barcode_rn_N21_2 = {v: k for k, v in rn_barcode_N21_2.iteritems()}
rn_barcode_N21_2_dedup = {v: k for k, v in barcode_rn_N21_2.iteritems()}
len(set(rn_barcode_N21_2_dedup.values())) # 500414

# Counting like above
lib140_cnt = collections.Counter()

for i in rn_barcode_N10_dedup.values():
  n = i[0]
  lib140_cnt[n] += 1

lib140_cnt['G'] # 96033
lib140_cnt['A'] + lib140_cnt['C'] + lib140_cnt['T'] # 180038
lib140_cnt['A'] + lib140_cnt['C'] + lib140_cnt['G'] + lib140_cnt['T'] # 276071
100*float(lib140_cnt['G'])/(lib140_cnt['A'] + lib140_cnt['C'] + lib140_cnt['G'] + lib140_cnt['T']) # 35% (34.78561674351888)
100*float(lib140_cnt['A'] + lib140_cnt['C'] + lib140_cnt['T'])/(lib140_cnt['A'] + lib140_cnt['C'] + lib140_cnt['G'] + lib140_cnt['T']) # 65% (65.21438325648113)

lib140_s = []

for i in range(1, 10001):
  random.seed(i)
  print(i)
  rn_barcode_N10_dedup_s = random.sample(rn_barcode_N10_dedup.values(), 1560)
  lib140_cnt_s = collections.Counter()
  for i in rn_barcode_N10_dedup_s:
    n = i[0]
    lib140_cnt_s[n] += 1
  tg_pct = 100*float(lib140_cnt_s['G'])/(lib140_cnt_s['A'] + lib140_cnt_s['C'] + lib140_cnt_s['G'] + lib140_cnt_s['T'])
  tx_pct = 100*float(lib140_cnt_s['A'] + lib140_cnt_s['C'] + lib140_cnt_s['T'])/(lib140_cnt_s['A'] + lib140_cnt_s['C'] + lib140_cnt_s['G'] + lib140_cnt_s['T'])
  lib140_s.append((tg_pct, tx_pct))


stats.describe([s[0] for s in lib140_s]) # DescribeResult(nobs=10000, minmax=(30.256410256410255, 40.0), mean=34.770615384615382, variance=1.4294577610292258, skewness=0.055306571905339115, kurtosis=0.02198226879298204)
stats.describe([s[1] for s in lib140_s]) # DescribeResult(nobs=10000, minmax=(60.0, 69.743589743589737), mean=65.229384615384618, variance=1.4294577610292258, skewness=-0.05530657190533867, kurtosis=0.02198226879298293)
stats.ttest_1samp([s[0] for s in lib140_s], stats.describe([s[0] for s in lib138_s]).mean) # pvalue=0.0
stats.ttest_ind([s[0] for s in lib140_s], [s[0] for s in lib138_s]) # pvalue=0.0
stats.ttest_1samp([s[0] for s in lib140_s], stats.describe([s[0] for s in lib139_s]).mean) # pvalue=0.0
stats.ttest_ind([s[0] for s in lib140_s], [s[0] for s in lib139_s]) # pvalue=0.0



# Is the difference of differences significant?
# Lib96+Lib100 (1)   -   Ts genome-wide (10000)
# Lib138 (mean of 10000)   -   Lib140 (10000)
# Lib139 (mean of 10000)   -   Lib140 (10000)

lib96_lib100_r1_fwd_tg = 100*float(lib96_lib100_r1_fwd_cnt['tg'])/(lib96_lib100_r1_fwd_cnt['ta'] + lib96_lib100_r1_fwd_cnt['tc'] + lib96_lib100_r1_fwd_cnt['tg'] + lib96_lib100_r1_fwd_cnt['tt'])
lib138_tg = stats.describe([s[0] for s in lib138_s]).mean
lib139_tg = stats.describe([s[0] for s in lib139_s]).mean
stats.ttest_ind(lib96_lib100_r1_fwd_tg - numpy.array([s[0] for s in lmajor_genome_s]), lib138_tg - numpy.array([s[0] for s in lib140_s])) # Ttest_indResult(statistic=1585.4222082646975, pvalue=0.0)
stats.ttest_ind(lib96_lib100_r1_fwd_tg / numpy.array([s[0] for s in lmajor_genome_s]), lib138_tg / numpy.array([s[0] for s in lib140_s])) # Ttest_indResult(statistic=985.71895480418505, pvalue=0.0)
stats.ttest_ind(lib96_lib100_r1_fwd_tg - numpy.array([s[0] for s in lmajor_genome_s]), lib139_tg - numpy.array([s[0] for s in lib140_s])) # Ttest_indResult(statistic=1806.159034591137, pvalue=0.0)
stats.ttest_ind(lib96_lib100_r1_fwd_tg / numpy.array([s[0] for s in lmajor_genome_s]), lib139_tg / numpy.array([s[0] for s in lib140_s])) # Ttest_indResult(statistic=1157.5360338049106, pvalue=0.0)
stats.ttest_ind(lib138_tg - numpy.array([s[0] for s in lib140_s]), lib139_tg - numpy.array([s[0] for s in lib140_s])) # Ttest_indResult(statistic=223.10959899938607, pvalue=0.0)
stats.ttest_ind(lib138_tg / numpy.array([s[0] for s in lib140_s]), lib139_tg / numpy.array([s[0] for s in lib140_s])) # Ttest_indResult(statistic=202.64554735085812, pvalue=0.0)
stats.ttest_1samp(lib139_tg / numpy.array([s[0] for s in lib140_s]), 1) # Ttest_1sampResult(statistic=124.44687391129081, pvalue=0.0)



# Output files for visualisation
df = pd.DataFrame({"lib96_lib100_diff":lib96_lib100_r1_fwd_tg - numpy.array([s[0] for s in lmajor_genome_s]), "lib96_lib100_fc":lib96_lib100_r1_fwd_tg / numpy.array([s[0] for s in lmajor_genome_s]), "lib138_diff": lib138_tg - numpy.array([s[0] for s in lib140_s]), "lib138_fc": lib138_tg / numpy.array([s[0] for s in lib140_s]), "lib139_diff": lib139_tg - numpy.array([s[0] for s in lib140_s]), "lib139_fc": lib139_tg / numpy.array([s[0] for s in lib140_s])})
df.to_csv("~/tg/tables/tg_diff_fc.txt", sep='\t', index=False)
```

Visualisation:

```r
library(ggplot2)
library(data.table)

# Load data
data <- fread("~/tg/tables/tg_diff_fc.txt")

# Split data in diff and fc
data_diff <- data[,grepl("diff", colnames(data)), with=FALSE]
data_fc <- data[,grepl("fc", colnames(data)), with=FALSE]

########
# diff #
########

data_diff_melt <- melt(data_diff, variable.name = "lib", value.name = "diff")

gg <- ggplot(data_diff_melt, aes(x = factor(lib, levels=c("lib96_lib100_diff", "lib138_diff", "lib139_diff")), y = diff)) +
geom_jitter(colour = 'gray', alpha = 0.35, width = 0.3, size = 0.2) +
geom_boxplot(outlier.shape=NA, alpha = 0) +
coord_cartesian(ylim = c(0, 40)) +
ylab("%TG difference") +
xlab("") +
theme_classic() +
theme(axis.title.y=element_text(size=20), axis.text.y=element_text(size=16), axis.text.x=element_text(size=16)) +
scale_x_discrete(labels=c("lib96_lib100_diff" = "Lib96+100\nvs\ngenome-wide", "lib138_diff" = "Lib138\nvs\nLib140", "lib139_diff" = "Lib139\nvs\nLib140"))

ggsave('~/figures/tg_diff.png')

######
# fc #
######

data_fc_melt <- melt(data_fc, variable.name = "lib", value.name = "fc")

gg <- ggplot(data_fc_melt, aes(x = factor(lib, levels=c("lib96_lib100_fc", "lib138_fc", "lib139_fc")), y = fc)) +
geom_hline(yintercept = 1, linetype = "dotted") +
geom_jitter(colour = 'gray', alpha = 0.35, width = 0.3, size = 0.2) +
geom_boxplot(outlier.shape=NA, alpha = 0) +
coord_cartesian(ylim = c(0.5, 2.5)) +
ylab("%TG fold-change") +
xlab("") +
theme_classic() +
theme(axis.title.y=element_text(size=20), axis.text.y=element_text(size=16), axis.text.x=element_text(size=16)) +
scale_x_discrete(labels=c("lib96_lib100_fc" = "Lib96+100\nvs\ngenome-wide", "lib138_fc" = "Lib138\nvs\nLib140", "lib139_fc" = "Lib139\nvs\nLib140"))

ggsave('~/figures/tg_fc.png')
```



## Peak calling

### macs2

```bash
cd ~/bam

mkdir -p ../macs2/default


# Lib96_Leish_SMUG1-AP vs. Lib99b_Leish_input_Y
bam_t=Lib96_Leish_SMUG1-AP_S4.clean.bam
bam_c=Lib99b_Leish_input_Y_S5.clean.bam

bname=`basename ${bam_t%.bam}`

macs2 callpeak \
-t $bam_t \
-c $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2/default/ \
-n $bname.nomodel.p1e-5 \
--tempdir ~/tmp/ \
--nomodel \
-p 0.00001 \
--cutoff-analysis


# Lib97_Leish_AP_S7 vs. Lib99b_Leish_input_Y
bam_t=Lib97_Leish_AP_S7.clean.bam
bam_c=Lib99b_Leish_input_Y_S5.clean.bam

bname=`basename ${bam_t%.bam}`

macs2 callpeak \
-t $bam_t \
-c $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2/default/ \
-n $bname.nomodel.p1e-5 \
--tempdir ~/tmp/ \
--nomodel \
-p 0.00001 \
--cutoff-analysis


# Lib100_Leish_SMUG1-AP_2 vs. Lib103_Leish_input-Y
bam_t=Lib100_Leish_SMUG1-AP_2_S7.clean.bam
bam_c=Lib103_Leish_input-Y_S8.clean.bam

bname=`basename ${bam_t%.bam}`

macs2 callpeak \
-t $bam_t \
-c $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2/default/ \
-n $bname.nomodel.p1e-5 \
--tempdir ~/tmp/ \
--nomodel \
-p 0.00001 \
--cutoff-analysis


# Lib101_Leish_AP_2 vs. Lib103_Leish_input-Y
bam_t=Lib101_Leish_AP_2_S5.clean.bam
bam_c=Lib103_Leish_input-Y_S8.clean.bam

bname=`basename ${bam_t%.bam}`

macs2 callpeak \
-t $bam_t \
-c $bam_c \
-f BAMPE \
--keep-dup all \
--outdir ../macs2/default/ \
-n $bname.nomodel.p1e-5 \
--tempdir ~/tmp/ \
--nomodel \
-p 0.00001 \
--cutoff-analysis
```



## Consensus peaks

```bash
cd ~/macs2/default

bedtools intersect \
-a Lib96_Leish_SMUG1-AP_S4.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib100_Leish_SMUG1-AP_2_S7.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 173, 100*173/214 = 80%

bedtools intersect \
-a Lib100_Leish_SMUG1-AP_2_S7.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib96_Leish_SMUG1-AP_S4.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 184, 100*184/196 = 94%


# multiinter
bedtools multiinter -i \
Lib96_Leish_SMUG1-AP_S4.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib100_Leish_SMUG1-AP_2_S7.clean.nomodel.p1e-5_peaks.narrowPeak \
-names Lib96 Lib100 | \
cut -f5 | sort | uniq -c | sort -k1,1nr
#    340 Lib96
#    185 Lib96,Lib100
#     52 Lib100


# merge
tableCat.py -i \
Lib96_Leish_SMUG1-AP_S4.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib100_Leish_SMUG1-AP_2_S7.clean.nomodel.p1e-5_peaks.narrowPeak | \
awk -v OFS="\t" '{print $1, $2, $3, $NF}' | \
sort -k1,1 -k2,2n | \
bedtools merge -c 4,4 -o distinct,count_distinct -i - | \
awk -v OFS="\t" '$5 > 1 {print $1, $2, $3}' > Lib96_Lib100_merge_p1e-5.bed

wc -l Lib96_Lib100_merge_p1e-5.bed # 172


bedtools intersect \
-a Lib97_Leish_AP_S7.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib101_Leish_AP_2_S5.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 4, 100*4/49 = 8%

bedtools intersect \
-a Lib101_Leish_AP_2_S5.clean.nomodel.p1e-5_peaks.narrowPeak \
-b Lib97_Leish_AP_S7.clean.nomodel.p1e-5_peaks.narrowPeak \
-sorted \
-wa -u | wc -l # 4, 100*4/102 = 4%


# multiinter
bedtools multiinter -i \
Lib97_Leish_AP_S7.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib101_Leish_AP_2_S5.clean.nomodel.p1e-5_peaks.narrowPeak \
-names Lib97 Lib101 | \
cut -f5 | sort | uniq -c | sort -k1,1nr
#    101 Lib101
#     50 Lib97
#      4 Lib97,Lib101


# merge
tableCat.py -i \
Lib97_Leish_AP_S7.clean.nomodel.p1e-5_peaks.narrowPeak \
Lib101_Leish_AP_2_S5.clean.nomodel.p1e-5_peaks.narrowPeak | \
awk -v OFS="\t" '{print $1, $2, $3, $NF}' | \
sort -k1,1 -k2,2n | \
bedtools merge -c 4,4 -o distinct,count_distinct -i - | \
awk -v OFS="\t" '$5 > 1 {print $1, $2, $3}' > Lib97_Lib101_merge_p1e-5.bed

wc -l Lib97_Lib101_merge_p1e-5.bed # 4
```



## Overlap with 5hmU and baseJ peaks

```bash
cd ~/macs2

wc -l ~/bed/fk_hmU.bed # 139
cat ~/bed/fk_hmU.bed | awk '{print $3-$2}' | awk '{s+=$1} END {print s}' # 260279
wc -l ~/bed/fk_hmU_baseJ.bed # 211
cat ~/bed/fk_hmU_baseJ.bed | awk '{print $3-$2}' | awk '{s+=$1} END {print s}' # 300735

for bed1 in default/Lib96_Lib100_merge_p1e-5.bed \
default/Lib97_Lib101_merge_p1e-5.bed
do
  echo $bed1, `cat $bed1 | wc -l`, `cat $bed1 | awk '{print $3-$2}' | awk '{s+=$1} END {print s}'`
  for bed2 in ~/bed/fk_hmU.bed ~/bed/fk_hmU_baseJ.bed
  do
    bedtools intersect \
    -sorted \
    -a <(bedtools sort -i $bed1) \
    -b <(cut -f1-3 $bed2 | bedtools sort -i -) \
    -wa -u | wc -l
    bedtools intersect \
    -sorted \
    -a <(bedtools sort -i $bed1) \
    -b <(cut -f1-3 $bed2 | bedtools sort -i -) \
    -wo | awk '{s+=$7} END {print s}'
  done
done
```
