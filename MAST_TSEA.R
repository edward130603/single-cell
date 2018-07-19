#MAST Transcript set enrichment analysis

##load libraries
library("MultiAssayExperiment")
library("SummarizedExperiment")
library(MAST)
library(parallel)
options(mc.cores = detectCores() - 1)
library(data.table)

##load data
data1 = readRDS("data/GSE63818-GPL16791.rds")
data1 = readRDS(gzcon(url("http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE63818-GPL16791.rds")))

##extract 
data1_tx = experiments(data1)[["tx"]] #transcripts
cData = colData(data1) #column info
fData = rowData(data1@ExperimentList$tx) #transcript info
mData = metadata(data1) #metadata

##prepare data
scaRaw = FromMatrix(log2(assays(data1_tx)[["TPM"]]+1), cData, fData)

filterCrit = mData$salmon_summary$percent_mapped > 60
sca = subset(scaRaw,filterCrit) #quality
sca = sca[which(freq(sca)> 0.2), ]

cdr2 = colSums(assay(sca)>0)
colData(sca)$cngeneson = scale(cdr2)

sca = subset(sca, sca$characteristics_ch1 == "developmental stage: 7 week gestation")

##run model
zlmCond <- zlm(~source_name_ch1 + cngeneson, sca)
zlm_name = paste0("tmp/zlm_", Sys.Date(), ".RData")
save(zlmCond, file = zlm_name)
load(zlm_name)

##bootstrap
boots <- bootVcov1(zlmCond, R = 50)
boots_name = paste0("tmp/boots_", Sys.Date(), ".RData")
save(boots, file = boots_name)
load(boots_name)

##create gene sets
fData2 = as.data.table(fData[,1:2])
setkey(fData2, gene)
fData2$transcript = as.character(fData2$transcript)
fData2$gene = as.character(fData2$gene)

genes = unique(fData2$gene)

sets = vector("list", length(unique(fData2$gene)))
#for (i in genes){
#  sets[[i]] = test$transcript[test$gene == i]
#}

for (i in genes){
  sets[[i]] = which(mcols(sca2)$gene == i)
}

sets = sets[names(sets) != ""]
sets = sets[lapply(sets,length)>0]
gsea = gseaAfterBoot(zlmCond, boots, sets, CoefficientHypothesis('source_name_ch1Somatic Cells'))

load("sets.RData")
load("gsea_7-1.RData")
z_stat_comb <- summary(gsea)
sigModules <- z_stat_comb[combined_adj<.05]
gseaTable <- melt(sigModules[,.(set, disc_Z, cont_Z, combined_Z)], id.vars='set')


