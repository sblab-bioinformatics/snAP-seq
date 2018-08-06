
## Software requirements

- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- Standard Unix tools: awk, sort, uniq
- [python v2.7.12](https://www.python.org/). Libraries:
  - [os](https://docs.python.org/2/library/os.html)
  - [re](https://docs.python.org/2/library/re.html)
  - [gzip](https://docs.python.org/2/library/gzip.html)
  - [collections](https://docs.python.org/2/library/collections.html)
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [ggplot2 v2.2.1](http://ggplot2.org/)
  - [ggseqlogo v0.1](https://cran.rstudio.com/web/packages/ggseqlogo/index.html)
  - [Biostrings v2.42.1](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
  - [scales v0.4.1](https://cran.r-project.org/web/packages/scales/index.html)




## Quality check

```bash
cd ~/fastq # directory containing the fastq sequencing files
mkdir ../fastqc
fastqc --noextract -q -o ../fastqc *.fastq.gz
```



## Trim adapters and filter bases based on quality

```bash
cd ~/fastq
mkdir ../fastq_trimmed

# Single-end (Lib16 and Lib26)
for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  cutadapt -a AGATCGGAAGAGC -m 15 -q 20 -o ../fastq_trimmed/$fq $fq > ../fastq_trimmed/$bname.txt
done


# Paired-end (Lib22 and Lib138-140)
for fq1 in *R1*.fastq.gz
do
  fq2=${fq1/_R1_/_R2_}
  bname=${fq1%.fastq.gz}
  cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 15 -q 20 -o ../fastq_trimmed/$fq1 -p ../fastq_trimmed/$fq2 $fq1 $fq2 > ../fastq_trimmed/${bname}.txt
done
```



## Alignment

### Prepare and index reference

```bash
cd ~/reference
cat spikeins.fa
#>AP1
#CACACCGCCAGCCACAGCAACGAACGTGCAGCGCCCCTCACGCCACAGAACATCGCATTTACGACGATTGATGTACTAAATAGTGGGTGGTCGGTTCGCG
#>GCAT1
#GGCCACCACCCGCACATACTCTGGTACGATTACGAACACAGCCCGACACCACCTCTAATGAACGTCGCTTATAGTGATTAACGCCCCGTAGACACCATGG
#>fU1
#GTGGCTGGTGTGCTGTCGGTCTAAGTCTAAGCAGTTCGGTTGGTCTGGGCTGTCTGGGGACATAAGCAAGGAAGCAAGGCAGACACGAGACCACGCAGCCAAGCA
#>fC1
#TTACGGGAGTGTGCTGGGCGTAGCGGTAGCTGACGTCGTGCGGGTCTTTGGGCTGTGTAGAGCGTATGGAGGAATGAGTGTGGAGGAATGGTAGATAGGTAGTTA
#>AP2
#GAACCTACGCATACCTAGCTGCCTAGAATCCTGCAGCGGCGTCTCGCTGGCTACCGAATCGTTAGGCTCATTGGAGTCTAGGAGGGAATTCGCACTCCCTGACGC
#>GCAT2
#GGCCACCACCTGCACATACTCCATCTGTAGAGACGGTGGACACGTCCGAGTACGCTGGCATACCCGTAGTACTCTGCCTAATGAGGAAGCTACTTCCACCCACGG
#>fU2
#TAGCCATACTGCCTCGTCCGGACACGTCAGGAGGAAAGCCAAGACACACGAACCAAGAGAACCAAGCAAGACAGAAGAGCACAAGCAGACCAGCGAACAG
#>fC2
#TCGCACACGCTCAGTCAGGTAGAGATCTAGGAGGGTGGAGAGGTGGTTGGAGAGGGTTAGGAGGAAGAGTGAGGTAGTGAGAGGGTGGAGGTGAGTGAGG

bwa index spikeins.fa
```


### Align, convert to bam, sort and index

```bash
cd ~/fastq_trimmed
mkdir ../bam

# Single-end (Lib16 and Lib26)
for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  bwa mem -t 4 -M ../reference/spikeins.fa $fq | samtools view -b - | samtools sort -T ../tmp/$bname -o ../bam/$bname.bam - &&
  samtools index ../bam/$bname.bam
done


# Paired-end (Lib22)
for fq1 in *R1*.fastq.gz
do
  fq2=${fq1/_R1_/_R2_}
  bname=${fq1%.fastq.gz}
  bwa mem -t 4 -M ../reference/spikeins.fa $fq1 $fq2 | samtools view -b - | samtools sort -T ../tmp/$bname -o ../bam/$bname.bam - &&
  samtools index ../bam/$bname.bam
done
```



## Distribution of reference sequences

```bash
cd ~/bam

# Single-end (Lib16 and Lib26)
for bam in *.bam
do
  bname=${bam%.bam}
  samtools view $bam -F 256 -F 4 -F 16 | cut -f3 | sort | uniq -c > ${bname}_fwd.txt
  samtools view $bam -F 256 -F 4 -f 16 | cut -f3 | sort | uniq -c > ${bname}_rev.txt
done


# Paired-end (Lib22)
for bam in *.bam
do
  bname=${bam%.bam}
  samtools view $bam -F260 -f64 -f32 | cut -f3 | sort | uniq -c > ${bname}_R1_fwd.txt
  samtools view $bam -F260 -f64 -f16 | cut -f3 | sort | uniq -c > ${bname}_R1_rev.txt
  samtools view $bam -F260 -f128 -f32 | cut -f3 | sort | uniq -c > ${bname}_R2_fwd.txt
  samtools view $bam -F260 -f128 -f16 | cut -f3 | sort | uniq -c > ${bname}_R2_rev.txt
done
```



## Patterns of alignment start sites

```bash
cd ~/bam

# Single-end (Lib16 and Lib26)
for bam in *.bam
do
  bname=${bam%.bam}
   for ref in AP1 AP2 fC1 fC2 fU1 fU2 GCAT1 GCAT2
   do
      samtools view $bam -F 256 -F 4 -F 16 | awk -v ref=$ref '$3==ref {print $4}' | sort -k1,1n | uniq -c | awk ' { t = $1; $1 = $2; $2 = t; print; } ' >  ${bname}_${ref}_fwd.txt
      samtools view $bam -F 256 -F 4 -f 16 | awk -v ref=$ref '$3==ref {print $4}' | sort -k1,1n | uniq -c | awk ' { t = $1; $1 = $2; $2 = t; print; } ' >  ${bname}_${ref}_rev.txt
   done
done


# Paired-end (Lib22)
for bam in *.bam
do
  for ref in AP1 AP2 fC1 fC2 fU1 fU2 GCAT1 GCAT2
  do
    bname=${bam%.bam}
    samtools view $bam -F260 -f64 -f32 | awk -v ref=$ref '$3==ref {print $4}' | sort -k1,1n | uniq -c | awk '{t = $1; $1 = $2; $2 = t; print; }' > ${bname}_R1_fwd_${ref}.txt &
    samtools view $bam -F260 -f64 -f16 | awk -v ref=$ref '$3==ref {print $4}' | sort -k1,1n | uniq -c | awk '{ t = $1; $1 = $2; $2 = t; print; }' > ${bname}_R1_rev_${ref}.txt &
    samtools view $bam -F260 -f128 -f32 | awk -v ref=$ref '$3==ref {print $4}' | sort -k1,1n | uniq -c | awk '{ t = $1; $1 = $2; $2 = t; print; }' > ${bname}_R2_fwd_${ref}.txt &
    samtools view $bam -F260 -f128 -f16 | awk -v ref=$ref '$3==ref {print $4}' | sort -k1,1n | uniq -c | awk '{ t = $1; $1 = $2; $2 = t; print; }' > ${bname}_R2_rev_${ref}.txt &
    done
done
```



## Synthetic 5-hmU N-oligo libraries

### Analysis

```bash
                                     *
fwd:   GTCTACCTGAACGCCGCTGTNNNNNNNNNNUNNNNNNNNNNGTAGTAGTCGACTAGACGTCCAACCAACGGAAGGGTATTCGGACGAGGCAGTATGGCTA
                                                                            *
rev:   TAGCCATACTGCCTCGTCCGAATACCCTTCCGTTGGTTGGACGTCTAGTCGACTACTACNNNNNNNNNNANNNNNNNNNNACAGCGGCGTTCAGGTAGAC
```

```python
import gzip
import os
import re
from collections import Counter, defaultdict

##########
# Lib138 #
##########

# Define variables
d = "fastq_trimmed"
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

# Common barcodes
cnt = Counter()
for s in rn_barcode.values():
  cnt[s] += 1

cnt.most_common(10) # [('GGGGAGGGAG', 59), ('GGGAGGGAGG', 59), ('GGAGGAGGGG', 56), ('GAAGGGAGGG', 56), ('GGAAGGGAGG', 56), ('GGGGAGGGGG', 55), ('GGAGGGGGAG', 55), ('GAGGGGAGGG', 55), ('GGGGAGGGGA', 54), ('GGAGGGGAGG', 54)]

# Deduplicate dictionary
barcode_rn = {v: k for k, v in rn_barcode.iteritems()}
rn_barcode_dedup = {v: k for k, v in barcode_rn.iteritems()}

# Create fasta file
fasta = ''
for rn in rn_barcode_dedup:
  fasta = fasta + ">%s\n%s\n" % (rn, rn_barcode_dedup[rn])

ofasta = open("~/fasta/Lib138.fasta", "w")
ofasta.write(fasta)
ofasta.close()

# Create a file containing ACGT counts
index_base = defaultdict(int)

for s in rn_barcode_dedup.values():
  for i in range(0, len(s)):
    index_base[(i+1, list(s)[i])] += 1


ofile = open("~/tables/Lib138.txt", "w")

ofile.write("\t" + "\t".join(["A", "C", "G", "T"]) + "\n")

for i in range(1, 11):
  ofile.write("%s\t" % str(i) + "\t".join([str(index_base[(i, b)]) for b in ["A", "C", "G", "T"]]) + "\n")


ofile.close()  



##########
# Lib139 #
##########

# Define variables
d = "~/fastq_trimmed"
files = os.listdir(d)
lib = "Lib139"

# Load reads into memory
reads = "".join([gzip.open(d + "/" + f, 'rt').read() for f in files if (lib in f) and ("_R2_001.fastq.gz" in f)]).split("@")[1:]
len(reads) # 658699, ok

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


# Common N10 barcodes
cnt = Counter()
for s in rn_barcode_N10.values():
  cnt[s] += 1

cnt.most_common(10) # [('GGGGGAGGGG', 120), ('GGGGGGGAGG', 105), ('GGGGGGAGGG', 104), ('GGAGGGGGGG', 101), ('GGGGAGGGGG', 100), ('GGGGGGGGAG', 89), ('GGGAGGGGAG', 82), ('GGAAGGGGGG', 81), ('GGGGAAGGGG', 78), ('GGGGAGAGGG', 78)]

# Deduplicate N10 dictionary
barcode_rn_N10 = {v: k for k, v in rn_barcode_N10.iteritems()}
rn_barcode_N10_dedup = {v: k for k, v in barcode_rn_N10.iteritems()}
len(set(rn_barcode_N10_dedup.values())) # 281120

# Reduce N21 dictionary based on the deduplication of N10 dictionary
rn_barcode_N21_1_dedup = {rn: rn_barcode_N21_1[rn] for rn in rn_barcode_N10_dedup.keys() if len(rn_barcode_N21_1[rn]) == 21 and "N" not in rn_barcode_N21_1[rn]}
len(rn_barcode_N21_1_dedup) # 276402
len(set(rn_barcode_N21_1_dedup.values())) # 276402

# Create fasta file with deduplicated N21 dictionary based on the deduplication of N10 dictionary
fasta = ''
for rn in rn_barcode_N21_1_dedup:
  fasta = fasta + ">%s\n%s\n" % (rn, rn_barcode_N21_1_dedup[rn])

ofasta = open("~/fasta/Lib139.fasta", "w")
ofasta.write(fasta)
ofasta.close()

# Create a file containing ACGT counts based on the deduplication of N10 dictionary
index_base = defaultdict(int)

for s in rn_barcode_N21_1_dedup.values():
  for i in range(0, len(s)):
    index_base[(i+1, list(s)[i])] += 1


ofile = open("~/tables/Lib139.txt", "w")

ofile.write("\t" + "\t".join(["A", "C", "G", "T"]) + "\n")

for i in range(1, 22):
  ofile.write("%s\t" % str(i) + "\t".join([str(index_base[(i, b)]) for b in ["A", "C", "G", "T"]]) + "\n")


ofile.close()  


# Common N21 barcodes
cnt = Counter()
for s in rn_barcode_N21_2.values():
  cnt[s] += 1

cnt.most_common(10) # [('GAAGGGTAAGTGGGTTGCACA', 3), ('GGTCCACCGCTCAAGGGTCGC', 3), ('AGTTAGATTCTCAGTGGGCGG', 3), ('TACCTGAACGCCGCTGTGGAG', 3), ('TCTACCTGAACGCCGCTGTAG', 3), ('AGCGCAGGGTTGAGAGTGGGA', 3), ('AGAGAGAAAGTACACTGGGGA', 3), ('TGGGTAGTAGTGGAACTGGAG', 3), ('TTCACACGAATCACCTCCCGT', 3), ('CAGCAGGAGGTACACAGATAG', 3)]

# Deduplicate N21 dictionary
barcode_rn_N21_2 = {v: k for k, v in rn_barcode_N21_2.iteritems()}
rn_barcode_N21_2_dedup = {v: k for k, v in barcode_rn_N21_2.iteritems()}
len(set(rn_barcode_N21_2_dedup.values())) # 552169

# Create fasta file with deduplicated N21 dictionary
fasta = ''
for rn in rn_barcode_N21_2_dedup:
  fasta = fasta + ">%s\n%s\n" % (rn, rn_barcode_N21_2_dedup[rn])

ofasta = open("~/fasta/Lib139_N21.fasta", "w")
ofasta.write(fasta)
ofasta.close()

# Create a file containing ACGT counts
index_base = defaultdict(int)

for s in rn_barcode_N21_2_dedup.values():
  for i in range(0, len(s)):
    index_base[(i+1, list(s)[i])] += 1


ofile = open("~/tables/Lib139_N21.txt", "w")

ofile.write("\t" + "\t".join(["A", "C", "G", "T"]) + "\n")

for i in range(1, 22):
  ofile.write("%s\t" % str(i) + "\t".join([str(index_base[(i, b)]) for b in ["A", "C", "G", "T"]]) + "\n")


ofile.close()  



##########
# Lib140 #
##########

# Define variables
d = "fastq_trimmed"
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


# Common N10 barcodes
cnt = Counter()
for s in rn_barcode_N10.values():
  cnt[s] += 1

cnt.most_common(10) # [('GGGGGGAGGG', 96), ('GGGGGAGGGG', 91), ('GGGGGGGAGG', 84), ('GGGGAGGGGG', 83), ('GGAGGGGGGG', 83), ('GGAGGGAGGG', 82), ('GGGAGGGGGA', 77), ('GGAGGGGGAG', 74), ('GGGGGGAGAG', 74), ('GGGGGAAGGG', 72)]

# Deduplicate N10 dictionary
barcode_rn_N10 = {v: k for k, v in rn_barcode_N10.iteritems()}
rn_barcode_N10_dedup = {v: k for k, v in barcode_rn_N10.iteritems()}
len(set(rn_barcode_N10_dedup.values())) # 276071

# Reduce N21 dictionary based on the deduplication of N10 dictionary
rn_barcode_N21_1_dedup = {rn: rn_barcode_N21_1[rn] for rn in rn_barcode_N10_dedup.keys() if len(rn_barcode_N21_1[rn]) == 21 and "N" not in rn_barcode_N21_1[rn]}
len(rn_barcode_N21_1_dedup) # 252457
len(set(rn_barcode_N21_1_dedup.values())) # 252457

# Create fasta file with deduplicated N21 dictionary
fasta = ''
for rn in rn_barcode_N21_1_dedup:
  fasta = fasta + ">%s\n%s\n" % (rn, rn_barcode_N21_1_dedup[rn])

ofasta = open("~/fasta/Lib140.fasta", "w")
ofasta.write(fasta)
ofasta.close()

# Create a file containing ACGT counts
index_base = defaultdict(int)

for s in rn_barcode_N21_1_dedup.values():
  for i in range(0, len(s)):
    index_base[(i+1, list(s)[i])] += 1


ofile = open("~/tables/Lib140.txt", "w")

ofile.write("\t" + "\t".join(["A", "C", "G", "T"]) + "\n")

for i in range(1, 22):
  ofile.write("%s\t" % str(i) + "\t".join([str(index_base[(i, b)]) for b in ["A", "C", "G", "T"]]) + "\n")


ofile.close()


# Common N21 barcodes
cnt = Counter()
for s in rn_barcode_N21_2.values():
  cnt[s] += 1

cnt.most_common(10) # [('TTCTTCCCTTTAATGGGCACG', 3), ('ATTAGAAGGGTGGGGAAGAAA', 3), ('TCTACCTGAACGCCGCTGTTG', 3), ('TACCTGAACGCCGCTGTAGGA', 3), ('AAGCGCGGACTGGAAAGGAGG', 3), ('GCGTTGTTGATAACTATGGGG', 3), ('GTGGGGGAAGTGAGGAGGAGG', 3), ('GGGGACAGAGTATGAAGGAGG', 3), ('AATAGTTGAGTGGAAAGCTCG', 3), ('GAAATAAGGTTGGTTGATAGC', 3)]

# Deduplicate N21 dictionary
barcode_rn_N21_2 = {v: k for k, v in rn_barcode_N21_2.iteritems()}
rn_barcode_N21_2_dedup = {v: k for k, v in barcode_rn_N21_2.iteritems()}
len(set(rn_barcode_N21_2_dedup.values())) # 500414

# Create fasta file with deduplicated N21 dictionary
fasta = ''
for rn in rn_barcode_N21_2_dedup:
  fasta = fasta + ">%s\n%s\n" % (rn, rn_barcode_N21_2_dedup[rn])

ofasta = open("~/fasta/Lib140_N21.fasta", "w")
ofasta.write(fasta)
ofasta.close()

# Create a file containing ACGT counts
index_base = defaultdict(int)

for s in rn_barcode_N21_2_dedup.values():
  for i in range(0, len(s)):
    index_base[(i+1, list(s)[i])] += 1


ofile = open("~/tables/Lib140_N21.txt", "w")

ofile.write("\t" + "\t".join(["A", "C", "G", "T"]) + "\n")

for i in range(1, 22):
  ofile.write("%s\t" % str(i) + "\t".join([str(index_base[(i, b)]) for b in ["A", "C", "G", "T"]]) + "\n")


ofile.close()  
```



### Sequence logos

```r
# Load the required packages
require(ggplot2)
require(ggseqlogo)
require(Biostrings)
library(scales)

# Create custom colour scheme
# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
# show_col(hue_pal()(4))
cs <- make_col_scheme(chars=c('A', 'C', 'G', 'T'), cols=hue_pal()(4))

##########
# Lib138 #
##########

# Load sequences
Lib138_file <- readDNAStringSet("~/fasta/Lib138.fasta")
Lib138_seqs <- as.data.frame(Lib138_file)$x
length(Lib138_seqs) # 217029, ok

# Plot sequence logo bits and prob
gg <- ggplot() +
geom_logo(Lib138_seqs, method = 'prob', col_scheme = cs) +
theme_logo() +
theme(axis.line.y = element_line(color="black"))

ggsave('~/figures/Lib138_prob.png')


##########
# Lib139 #
##########

# Load sequences
Lib139_file <- readDNAStringSet("~/fasta/Lib139_N21.fasta")
Lib139_seqs <- as.data.frame(Lib139_file)$x
length(Lib139_seqs) # 552169, ok

# Plot sequence logo bits and prob
gg <- ggplot() +
geom_logo(Lib139_seqs, method = 'prob', col_scheme = cs) +
theme_logo() +
theme(axis.line.y = element_line(color="black"))

ggsave('~/figures/Lib139_prob.png', width = 14)


##########
# Lib140 #
##########

# Load sequences
Lib140_file <- readDNAStringSet("~/fasta/Lib140_N21.fasta")
Lib140_seqs <- as.data.frame(Lib140_file)$x
length(Lib140_seqs) # 500414, ok

# Plot sequence logo bits and prob
gg <- ggplot() +
geom_logo(Lib140_seqs, method = 'prob', col_scheme = cs) +
theme_logo() +
theme(axis.line.y = element_line(color="black"))

ggsave('~/figures/Lib140_prob.png', width = 14)
```
