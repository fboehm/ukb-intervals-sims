scenarios=( I II III IV)
distributions=( laplace normal scaledt)
hsqs=(0.1 0.2 0.5)

for scenario in ${scenarios[@]}; do
  for distribution in ${distributions[@]}; do
    for hsq in ${hsqs[@]}; do
      bed=~/research/ukb-intervals-sims/hapmap3/plink_files_for_sims/chr1.bed
      bim=~/research/ukb-intervals-sims/hapmap3/plink_files_for_sims/chr1.bim
      outdir=~/research/ukb-intervals-sims/dat-quant-5fold-manual-sim/gemma/scenario${scenario}/${distribution}/hsq${hsq}
      ln -s ${bed} ${outdir}/chr1.bed 
      ln -s ${bim} ${outdir}/chr1.bim 
    done
  done
done
