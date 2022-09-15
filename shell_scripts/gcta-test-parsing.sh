#!/bin/bash

#files= (../hapmap3/snp_effects/*)

declare -a files
for file in ~/research/ukb-intervals-sims/hapmap3/snp_effects/*.txt
do
    files=("${files[@]}" "$file")
done

#echo "${files[@]}"


for i in "${files[@]}"; do
  # parse hsq from file name
  # first, cut with period as field separator
  #filename=~/research/ukb-intervals-sims/hapmap3/snp_effects/scenarioIII_laplace_hsq0.5.txt
  filename=$i
  echo $filename
  field2=$( echo "${filename}" | cut -d'.' -f2 )
  echo "$field2"
  hsq=$( echo "scale=2 ; $field2 / 10" | bc )
  echo ${hsq}
  out_suffix=$( basename $filename )
  echo ${out_suffix}
done
