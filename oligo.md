
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
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)



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


### Comparing counts

```r
library(data.table)

# Enlarge the view width when printing tables
options(width = 300)

# Load libraries
lib138 <- fread("~/tables/Lib138.txt")
colnames(lib138)[1] <- "pos"

lib139 <- fread("~/tables/Lib139.txt")
colnames(lib139)[1] <- "pos"

lib140 <- fread("~/tables/Lib140.txt")
colnames(lib140)[1] <- "pos"


#####################
# Lib138 vs. Lib140 #
#####################
lib138_lib140 <- cbind(lib138[, .(G, other=A+C+T)], lib140[12:21, .(G, other=A+C+T)])
colnames(lib138_lib140) <- c("G.lib138", "other.lib138", "G.lib140", "other.lib140")

lib138_lib140[, odds := (G.lib138/other.lib138)/(G.lib140/other.lib140)]
lib138_lib140[, pval_f := fisher.test(matrix(c(G.lib138, other.lib138, G.lib140, other.lib140), nrow = 2))$p.value, by = 1:nrow(lib138_lib140)]
lib138_lib140[, pval_q := chisq.test(matrix(c(G.lib138, other.lib138, G.lib140, other.lib140), nrow = 2))$p.value, by = 1:nrow(lib138_lib140)]

lib138_lib140
#    G.lib138 other.lib138 G.lib140 other.lib140      odds       pval_f       pval_q
# 1:    86907       130122    86947       165510 1.2713750 0.000000e+00 0.000000e+00
# 2:    73677       143352    87948       164509 0.9613729 1.669730e-10 1.684780e-10
# 3:    77760       139269    88124       164333 1.0411957 4.327105e-11 4.345580e-11
# 4:    75973       141056    87471       164986 1.0158994 1.030426e-02 1.033280e-02
# 5:    76332       140697    88432       164025 1.0062883 3.071588e-01 3.081797e-01
# 6:    77800       139229    88390       164067 1.0372131 2.396764e-09 2.401747e-09
# 7:    78209       138820    88238       164219 1.0485097 9.872242e-15 9.810232e-15
# 8:    76809       140220    87739       164718 1.0283727 5.152470e-06 5.153196e-06
# 9:    78434       138595    88913       163544 1.0409411 5.233433e-11 5.212921e-11
#10:    77747       139282    88130       164327 1.0408156 6.522288e-11 6.487237e-11


#####################
# Lib139 vs. Lib140 #
#####################
lib139_lib140 <- cbind(lib139[12:21, .(G, other=A+C+T)], lib140[12:21, .(G, other=A+C+T)])
colnames(lib139_lib140) <- c("G.lib139", "other.lib139", "G.lib140", "other.lib140")

lib139_lib140[, odds := (G.lib139/other.lib139)/(G.lib140/other.lib140)]
lib139_lib140[, pval_f := fisher.test(matrix(c(G.lib139, other.lib139, G.lib140, other.lib140), nrow = 2))$p.value, by = 1:nrow(lib139_lib140)]
lib139_lib140[, pval_q := chisq.test(matrix(c(G.lib139, other.lib139, G.lib140, other.lib140), nrow = 2))$p.value, by = 1:nrow(lib139_lib140)]

lib139_lib140
#    G.lib139 other.lib139 G.lib140 other.lib140      odds       pval_f       pval_q
# 1:   100398       176004    86947       165510 1.0858558 2.145023e-46 2.290044e-46
# 2:    94743       181659    87948       164509 0.9755598 1.931301e-05 1.940241e-05
# 3:    95807       180595    88124       164333 0.9892865 6.270575e-02 6.279071e-02
# 4:    95541       180861    87471       164986 0.9963866 5.319651e-01 5.335171e-01
# 5:    96217       180185    88432       164025 0.9904528 9.684574e-02 9.725182e-02
# 6:    96132       180270    88390       164067 0.9898346 7.715034e-02 7.736764e-02
# 7:    96226       180176    88238       164219 0.9939470 2.944422e-01 2.945425e-01
# 8:    95703       180699    87739       164718 0.9943017 3.240534e-01 3.245475e-01
# 9:    97080       179322    88913       163544 0.9957849 4.640386e-01 4.655292e-01
#10:    96517       179885    88130       164327 1.0004468 9.401538e-01 9.406367e-01


#####################
# Lib138 vs. Lib139 #
#####################
lib138_lib139 <- cbind(lib138[, .(G, other=A+C+T)], lib139[12:21, .(G, other=A+C+T)])
colnames(lib138_lib139) <- c("G.lib138", "other.lib138", "G.lib139", "other.lib139")

lib138_lib139[, odds := (G.lib138/other.lib138)/(G.lib139/other.lib139)]
lib138_lib139[, pval_f := fisher.test(matrix(c(G.lib138, other.lib138, G.lib139, other.lib139), nrow = 2))$p.value, by = 1:nrow(lib138_lib139)]
lib138_lib139[, pval_q := chisq.test(matrix(c(G.lib138, other.lib138, G.lib139, other.lib139), nrow = 2))$p.value, by = 1:nrow(lib138_lib139)]

lib138_lib139
#    G.lib138 other.lib138 G.lib139 other.lib139      odds        pval_f        pval_q
# 1:    86907       130122   100398       176004 1.1708507 3.209171e-157 2.094121e-157
# 2:    73677       143352    94743       181659 0.9854576  1.554213e-02  1.560027e-02
# 3:    77760       139269    95807       180595 1.0524713  1.639867e-17  1.609480e-17
# 4:    75973       141056    95541       180861 1.0195836  1.285855e-03  1.289395e-03
# 5:    76332       140697    96217       180185 1.0159881  8.358264e-03  8.411812e-03
# 6:    77800       139229    96132       180270 1.0478651  6.580868e-15  6.633601e-15
# 7:    78209       138820    96226       180176 1.0548949  5.006841e-19  4.984608e-19
# 8:    76809       140220    95703       180699 1.0342662  2.111794e-08  2.118601e-08
# 9:    78434       138595    97080       179322 1.0453473  1.312203e-13  1.311646e-13
#10:    77747       139282    96517       179885 1.0403507  4.259812e-11  4.292397e-11
```
