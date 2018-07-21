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
load(sca_name)


tpm = t(assays(data1_tx)[["TPM"]]) #design matrix
cell_select = cData$geo_accession %in% getwellKey(sca)
tpm = tpm[cell_select, as.character(mcols(sca)$transcript)] #filter 7 week, select genes
y = as.numeric(cData$source_name_ch1[cell_select]) - 1 #0/1 binary response vector
y_samp = sample(y)
load(sets_name)



regs = rep(NA, length(sets)) #init vector for storing results
#LR test for full model vs just mean
for (i in 1:length(sets)){
  regs[i] = lrtest(glm(y~.,data = data.frame(tpm[,sets[[i]]]), family = binomial))[2,5] 
}
prop.table(table(regs < 0.05))
adjusted = p.adjust(regs, 'fdr')

out1 = data.frame(gene = names(sets), p = regs, fdr = adjusted) #output
out1 = out1[out1$fdr<0.05,] #select top genes


#venn.diagram(list(out1$gene, sigModules$set, fcHurdleSig_g$primerid), 
#             "test.png", 
#             category.names = c("Logistic regression", "MAST Transcript GSEA", "MAST Genes"),
#             margin = 0.2)

