# GOAL: fix directory structure in 03_subsample/binary

for p in `seq 1 25`; do
  mkdir -p ../03_subsample/binary/pheno${p}/val
  mkdir -p ../03_subsample/binary/pheno${p}/verif
  mv ../03_subsample/binary/pheno${p}/val_ukb ../03_subsample/binary/pheno${p}/val/ukb
  mv ../03_subsample/binary/pheno${p}/verif_ukb ../03_subsample/binary/pheno${p}/verif/ukb
done  
