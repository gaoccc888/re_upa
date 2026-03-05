###计算每个谱系等位基因频率
#!/bin/bash

for pop in north central southeast
do
   vcftools --gzvcf 311_filtered.recode_noscaffold.vcf.gz \
            --keep ${pop}.list \
            --freq \
            --out ${pop} &
done

wait
echo "All jobs finished."

python calculate_U.py \
    --north north.frq \
    --central central.frq \
    --southeast southeast.frq \
    --window 50000 \
    --step 10000 \
    --mode north2central
    #--threads 10
##north2central
North > 0.95
Southeast < 0.05
Central > 0.20
#mode
north2central
southeast2central
