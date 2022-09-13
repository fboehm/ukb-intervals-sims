#!/bin/bash

dummyfam=~/research/ukb-intervals-sims/dat/chr22.fam
newdir=~/research/ukb-intervals-sims/hapmap3/ukb
for chr in `seq 1 22`; do
        ln -s /net/mulan/disk2/yasheng/predictionProject/plink_file/ukb/chr${chr}.bed ${newdir}/chr${chr}.bed 
        ln -s /net/mulan/disk2/yasheng/predictionProject/plink_file/ukb/chr${chr}.bim ${newdir}/chr${chr}.bim
        ln -s ${dummyfam} ${newdir}/chr${chr}.fam
      done
