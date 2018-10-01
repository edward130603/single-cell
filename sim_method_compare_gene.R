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

##logistic regression
tpm = t(assays(sim.groups)[["tpm"]]) #design matrix
y = ifelse(cData$Group == "Group1", 0, 1)  #0/1 binary response vector
set.seed(100)
y_samp = sample(y)

regs = rep(NA, ncol(tpm)) #init vector for storing results
#LR test for full model vs just mean
for (i in 1:ncol(tpm)){
  regs[i] = lrtest(glm(y~.,data = data.frame(tpm[,i]), family = binomial))[2,5] 
}
names(regs) = colnames(tpm)

hist(regs, main = "Distribution of logistic regression p-values", xlab = "p")

fData %>%
  as.tibble() %>%
  mutate(de = DEFacGroup1 != DEFacGroup2,
         log_reg = regs[as.character(Gene)]) ->
  out1


roc(de~log_reg, out1, plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))

with(out1, table(de, p.adjust(log_reg, method = "fdr") < 0.05))

##MAST

sca = FromMatrix(log2(assays(sim.groups)[["tpm"]]+1), cData, fData)

cdr2 = colSums(assay(sca)>0)
colData(sca)$cngeneson = scale(cdr2)
set.seed(100)
colData(sca)$Group = sample(colData(sca)$Group)
sca = sca[which(freq(sca)> 0.2), ]

##run model
zlmCond = zlm(~Group + cngeneson, sca)
summaryCond = summary(zlmCond, doLRT='GroupGroup2') 
summaryDt = summaryCond$datatable
fcHurdle = merge(summaryDt[contrast=='GroupGroup2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                 summaryDt[contrast=='GroupGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
out1 = merge(out1, fcHurdle[,1:2], by.x = "Gene", by.y = "primerid")
colnames(out1)[10] = "mast"

hist(out1$mast, main = "Distribution of MAST p-values", xlab = "p")
hist(out1$log_reg, main = "Distribution of logistic regression p-values", xlab = "p")


roc(de~mast, out1, plot = T,
    print.auc = T, print.auc.x = 0.1, print.auc.y = 0.1,
    xlim = c(0,1))
with(out1, table(de, p.adjust(mast, method = "fdr") < 0.05))

##t test
t_p = sapply(as.character(out1$Gene), function(x){t.test(x = tpm[y==0, x], 
                                                   y = tpm[y==1, x], var.equal = T)$p.value})
t_p[is.nan(t_p)] = 1

#plot p-values
ggplot(data = out1, aes(x = mast, y = log_reg))+
  geom_point(aes(alpha=0.1), shape = 16) +
  guides(alpha = F)+
  labs(x = "MAST", y = "Logistic Regression", title = "P-values")
