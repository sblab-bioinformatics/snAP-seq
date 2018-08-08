
## Software requirements

- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [sambamba v0.6.5](https://academic.oup.com/bioinformatics/article/31/12/2032/214758)
- Standard Unix tools: awk, sort, uniq
- [igvtools v2.3.91](https://software.broadinstitute.org/software/igv/igvtools)
- [tableCat.py](https://github.com/dariober/bioinformatics-cafe/blob/master/tableCat/tableCat.py)
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [ggplot2 v2.2.1](http://ggplot2.org/)
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)
  - [edgeR v3.16.5](https://bioconductor.org/packages/release/bioc/html/edgeR.html)



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

First `clean.tdf`:

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



## Calling AP sites (SMUG1 and untreated)

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
data <- fread("tableCat.py -i ~/ap_caller_nosubtraction_genome/*Lib{97,101,99b,103}*fwd*.cov -r .clean.R1.fwd.cigar.cov")
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

##########################################
# Lib96 and Lib100 vs. Lib99b and Lib103 #
##########################################

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


##########################################
# Lib97 and Lib101 vs. Lib99b and Lib103 #
##########################################

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
```

https://github.com/sblab-bioinformatics/projects/blob/master/20161006_abasic_sites_jane/20180108_leish_human.md#compare-ap_caller_nosubtraction_genome_i1-results-with-5hmu-and-j-peaks-from-genome-biology-paper
