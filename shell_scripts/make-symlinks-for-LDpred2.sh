#!/bin/bash

hsq=0.1
pcausal=0.1
newdir=../dat/hsq${hsq}_pcausal${pcausal}/LDpred2
mkdir -p ${newdir}

for chr in `seq 1 22`; do
        ln -s /net/mulan/disk2/yasheng/predictionProject/plink_file/ukb/chr${chr}.bed ${newdir}/chr${chr}.bed 
        ln -s /net/mulan/disk2/yasheng/predictionProject/plink_file/ukb/chr${chr}.bim ${newdir}/chr${chr}.bim 
      done
