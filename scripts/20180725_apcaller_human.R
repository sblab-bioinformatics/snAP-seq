#!/usr/bin/env Rscript

library(data.table)
library(edgeR)

args = commandArgs(trailingOnly=TRUE)

# Load data
print("Load data")
data <- fread(paste("zcat", args[1]))
setnames(data, c("chr", "pos", "Lib134", "Lib135", "Lib145", "Lib146", "Lib136", "Lib137", "Lib147", "Lib148"))

# Define group
group <- factor(c('ap', 'ap', 'ap', 'ap', 'input', 'input', 'input', 'input'))

# Define DGEList object
print("Define DGEList object")
y <- DGEList(counts = data[,-c(1,2)], group = group, genes = data[,c(1,2)])

# Filter and get the top 100k sites according to cpm
print("Filter and get the top 100k sites according to cpm")
if(grepl("chrM", args[1])){
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
write.table(detable_e, paste("../tables/", gsub(".gz", ".txt", args[1]), sep=""), row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


