import subprocess
import re
import datetime
from subprocess import call

for i in range(10):
    shfile='MAST_'+str(datetime.date.today())+'_'+str(i)+'.sh'
    rfile ='MAST_'+str(datetime.date.today())+'_'+str(i)+'.R'

    x=open(shfile,'w+')
    x.write('#!/bin/bash\n')
    x.write('#$ -cwd\n')
    x.write('#$ -l mem_free=64G\n')
    x.write('#$ -l h_vmem=64G\n')
    x.write('#$ -q thornton.q\n')
    x.write('R CMD BATCH ' + rfile)
    x.close()

    r=open(rfile,'w+')
    r.write('library("MultiAssayExperiment")\n')
    r.write('library("SummarizedExperiment")\n')
    r.write('library(MAST)\n')
    r.write('library(parallel)\n')
    r.write('library(data.table)\n')
    r.write('load("tmp/sca_2018-07-24.RData")\n')
    r.write('colData(sca)$source_name_ch1 = sample(colData(sca)$source_name_ch1)\n')
    r.write('zlmCond = zlm(~source_name_ch1 + cngeneson, sca)\n')
    r.write('boots = bootVcov1(zlmCond, R = 50)\n')
    r.write('load("tmp/sets_2018-07-24.RData")\n')
    r.write('gsea = gseaAfterBoot(zlmCond, boots, sets, CoefficientHypothesis("source_name_ch1Somatic Cells"))\n')
    r.write('z_stat_comb <- summary(gsea) #all genes\n')
    r.write('save(z_stat_comb, file = "MAST_sample_'+str(datetime.date.today())+'_'+str(i)+'.RData")\n')
    r.close()
    toqsub="qsub " + shfile
    #call(toqsub,shell=True)
