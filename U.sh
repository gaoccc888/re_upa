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
