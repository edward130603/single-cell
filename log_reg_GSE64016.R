#Differential expression with logistic regression (Ntranos, et al.)

library("MultiAssayExperiment")
library("SummarizedExperiment") #for rowData()
library("lmtest") #for lrtest()

##load data
data1 = readRDS("data/GSE63818-GPL16791.rds") 
#download from "http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE63818-GPL16791.rds"

##extract
data1_tx = experiments(data1)[["tx"]] #transcripts
cData = colData(data1) #column info
fData = rowData(data1@ExperimentList$tx) #transcript info
load(sca_name) #from MAST_GSE64016.R


tpm = t(assays(data1_tx)[["TPM"]]) #design matrix
cell_select = cData$geo_accession %in% getwellKey(sca)
tpm = tpm[cell_select, as.character(mcols(sca)$transcript)] #select genes
y = as.numeric(droplevels(cData$source_name_ch1[cell_select])) - 1 #0/1 binary response vector

y_samp = sample(y)
load(sets_name) #from MAST_GSE64016.R

regs = rep(NA, length(sets)) #init vector for storing results
#LR test for full model vs just mean
for (i in 1:length(sets)){
  regs[i] = lrtest(glm(y~.,data = data.frame(tpm[,sets[[i]]]), family = binomial))[2,5] 
} #replace y_samp with y (or vice versa)

adjusted = p.adjust(regs, 'fdr')

out1 = data.frame(gene = names(sets), p = regs, fdr = adjusted) #output
out1 = out1[out1$fdr<0.05,] #select top genes

hist(regs, main = "Distribution of logistic regression p-values", xlab = "p")

#subsample
y_samp2 = y_samp[1:60]
regs2 = rep(NA, length(sets)) #init vector for storing results
#LR test for full model vs just mean
for (i in 1:length(sets)){
  regs2[i] = lrtest(glm(y_samp2~.,data = data.frame(tpm[1:60,sets[[i]]]), family = binomial))[2,5] 
}
hist(regs2, main = "Distribution of logistic regression p-values", xlab = "p")
