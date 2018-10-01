#MAST Transcript set enrichment analysis

##load libraries
library(splatter)
library(scater)
library(lmtest)
library(MAST)
library(dplyr)
library(tibble)
library(pROC)
options(mc.cores = detectCores() - 1)
library(data.table)


##simulate data
sim.groups = splatSimulate(group.prob = c(0.5, 0.5),
                           method = "groups")
sim.groups = addGeneLengths(sim.groups)
#sim.groups = normalise(sim.groups, return_log = F)

tpm(sim.groups) <- calculateTPM(sim.groups, rowData(sim.groups)$Length)

fData = rowData(sim.groups)
cData = colData(sim.groups)
table(fData$DEFacGroup1 == fData$DEFacGroup2)

load("C:/Users/Edward/Google Drive/IS_RG/single-cell/tmp/sets_2018-07-19.RData")
set_distribution = prop.table(table(sapply(sets,length)))
set.seed(123)

fData$Set = NA
j = 1
i = 1
while (i < nrow(fData)){
  set_size = sample(prob = set_distribution, 
                    x = as.numeric(names(set_distribution)),
                    replace = T, size = 1)
  fData$Set[i:min(nrow(fData),i+set_size-1)] = paste0("Set",j)
  j = j+1
  i = i+set_size-1
}

genes = unique(fData$Set)
sets = vector("list", length(genes))
names(sets) = genes

for (i in genes){
  sets[[i]] = which(fData$Set == i)
}

##logistic regression
tpm = t(assays(sim.groups)[["tpm"]]) #design matrix
y = ifelse(cData$Group == "Group1", 0, 1)  #0/1 binary response vector

regs = rep(NA, length(sets)) #init vector for storing results
#LR test for full model vs just mean
for (i in 1:length(sets)){
  regs[i] = lrtest(glm(y~.,data = data.frame(tpm[,sets[[i]]]), family = binomial))[2,5] 
}
names(regs) = names(sets)

hist(regs, main = "Distribution of logistic regression p-values", xlab = "p")

fData %>%
  as.tibble() %>%
  group_by(Set) %>%
  summarise(de = !all(DEFacGroup1 == DEFacGroup2),
            de2 = !all(DEFacGroup1/DEFacGroup2 > (1/1.5) & DEFacGroup1/DEFacGroup2 <1.5 ),
            keep = all(DEFacGroup1 == DEFacGroup2) | any(DEFacGroup1/DEFacGroup2 < (1/1.5)) | any(DEFacGroup1/DEFacGroup2 >1.5 )) %>%
  mutate(log_reg = regs[Set]) ->
  out1


roc(de~log_reg, out1, plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))


roc(de2~log_reg, out1 %>% filter(keep), plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))

with(out1 %>% filter(keep), table(de2, p.adjust(log_reg, method = "fdr") < 0.05))

##MAST-TSEA

sca = FromMatrix(log2(assays(sim.groups)[["tpm"]]+1), cData, fData)

cdr2 = colSums(assay(sca)>0)
colData(sca)$cngeneson = scale(cdr2)

##run model
zlmCond = zlm(~Group + cngeneson, sca)
##bootstrap
boots = bootVcov1(zlmCond, R = 50)
boots_name = paste0("tmp/boots_", Sys.Date(), ".RData")
save(boots, file = boots_name)
load(boots_name)


##gsea
gsea = gseaAfterBoot(zlmCond, boots, sets, CoefficientHypothesis('GroupGroup2'))
gsea_name = paste0("tmp/gsea_", Sys.Date(), ".RData")
save(gsea, file = gsea_name)
load(gsea_name)

z_stat_comb <- summary(gsea) #all genes
sigModules <- z_stat_comb[combined_adj<.05] #significant genes
hist(z_stat_comb$combined_P, main = "Distribution of MAST-TSEA p-values",
     xlab = "p")
out1$mast_tsea = z_stat_comb$combined_P[match(out1$Set, z_stat_comb$set)]
out1$mast_tsea[which(is.nan(out1$mast_tsea))] = 1
roc(de~mast_tsea, out1, plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))
roc(de2~mast_tsea, out1 %>% filter(keep), plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))

##mast-transcript
zlmCond = zlm(~Group + cngeneson, sca)
summaryCond = summary(zlmCond, doLRT='GroupGroup2') 
summaryDt = summaryCond$datatable
fcHurdle = merge(summaryDt[contrast=='GroupGroup2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                 summaryDt[contrast=='GroupGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
fcHurdle$fdr = p.adjust(fcHurdle$`Pr(>Chisq)`,method="fdr")
fcHurdle = merge(fcHurdle, fData, by.x = "primerid", by.y = "Gene")

fcHurdle %>%
  as.tibble() %>%
  group_by(Set) %>%
  summarise(n = n(),
            df = sum(GeneMean),
            sidak = 1-(1-min(Pr..Chisq.))^n,
            fisher = pchisq(-2* sum(log(Pr..Chisq.)) ,2*n, lower.tail = F),
            lancaster = pchisq(sum(qchisq(Pr..Chisq., GeneMean, lower.tail = F)) ,df, lower.tail = F),
            any = any(fdr < 0.05)) -> fcHurdle2

out1$mast_transcript = fcHurdle2$any
out1$mast_transcript2 = fcHurdle2$sidak
out1$mast_transcript3 = fcHurdle2$fisher
out1$mast_transcript4 = fcHurdle2$lancaster
with(out1 %>% filter(keep), table(de2, mast_transcript))
with(out1 %>% filter(keep), table(de2, p.adjust(mast_transcript2, method = "fdr") < 0.05))
with(out1, table(de, p.adjust(mast_transcript2, method = "fdr") < 0.05))
with(out1 %>% filter(keep), table(de2, p.adjust(mast_transcript3, method = "fdr") < 0.05))
with(out1, table(de, p.adjust(mast_transcript3, method = "fdr") < 0.05))
with(out1 %>% filter(keep), table(de2, p.adjust(mast_transcript4, method = "fdr") < 0.05))
with(out1, table(de, p.adjust(mast_transcript4, method = "fdr") < 0.05))
roc(de2~mast_transcript2, out1 %>% filter(keep), plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))
roc(de~mast_transcript2, out1, plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))
roc(de2~mast_transcript3, out1 %>% filter(keep), plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))
roc(de~mast_transcript3, out1, plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))
roc(de2~mast_transcript4, out1 %>% filter(keep), plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))

##mast-sum
tpm %>%
  t() %>%
  as.tibble() %>%
  mutate(set = rep(names(sets), sapply(sets, length))) %>%
  group_by(set) %>%
  summarise_all(sum) %>%
  column_to_rownames("set")->
  data1
cData$wellKey = cData$Cell
fData_gene = data.frame(primerid = rownames(data1))
sca = FromMatrix(log2(as.matrix(data1)+1), cData, fData_gene)
cdr2 = colSums(assay(sca)>0)
colData(sca)$cngeneson = scale(cdr2)
zlmCond = zlm(~Group + cngeneson, sca)
summaryCond = summary(zlmCond, doLRT='GroupGroup2') 
summaryDt = summaryCond$datatable
fcHurdle = merge(summaryDt[contrast=='GroupGroup2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='GroupGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
out1$mast = fcHurdle$`Pr(>Chisq)`
roc(de~mast, out1, plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))
roc(de2~mast, out1 %>% filter(keep), plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))
save(out1, file = "tmp/out1_2018-07-31.RData")

##removing freq < 0.2 genes
fData$freq = freq(sca)
fData%>%
  as.tibble() %>%
  group_by(Set) %>%
  filter(!all(freq<0.2)) %>%
  select(Set) %>%
  unique() ->
  freq

out2 = merge(out1, freq, by = "Set")

roc(de~mast, out2, plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))
roc(de~mast_tsea, out2, plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))
roc(de~mast_transcript2, out2, plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))
roc(de~log_reg, out2, plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))

##p-values
lengths = data.frame(Set = names(sets),length =  sapply(sets, length))
out1 = merge(out1, lengths, by = "Set")
out1$`Number of Transcripts` = ifelse(out1$length <= 5, as.character(out1$length), "6+")
ggplot(data = out1, aes(x = mast_transcript2, y = log_reg))+
  geom_point(aes(alpha = 0.1, color = `Number of Transcripts`), shape = 16) +
  guides(alpha = F)+
  scale_color_brewer() +
  labs(x = "MAST-Transcript", y = "Logistic Regression")
ggplot(data = out1, aes(x = mast_transcript2, y = log_reg))+
  geom_point(aes(alpha = 0.1), shape = 16) +
  guides(alpha = F)+
  labs(x = "MAST-Transcript", y = "Logistic Regression")
