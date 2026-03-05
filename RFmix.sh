# 使用 conda 安装
conda install -c bioconda rfmix

####提取ref.vcf
vcftools \
        --gzvcf 311_phased_all_genome.vcf.gz \
        --keep ref_sample.list \
        --recode \
        --recode-INFO-all \
        --stdout \
        >ref.vcf
####提取central_hyb.vcf
vcftools \
        --gzvcf 311_phased_all_genome.vcf.gz \
        --keep central_hyb_sample.list \
        --recode \
        --recode-INFO-all \
        --stdout \
        >query.vcf
###make the genetic map
cut -f1,2 Upa_chromosomes.fasta.fai > chr_length.txt

awk '{
    chr=$1
    len=$2
    print chr"\t0\t0.0"
    print chr"\t"len"\t"len/1000000
}' chr_length.txt > genetic.map

####cut -f1,2 query.vcf | grep -v "#" | awk '{print $0"\t"$2/650000}' >genetic.map

# get a list of all scaffolds
cut -f1 genetic.map | sort | uniq >chr.list

####按 染色体拆分的策略  ./run_rfmix.sh 
# 单条染色体任务函数
#!/bin/bash

for SCAFF in $(cat chr.list); do
    echo "========================================="
    echo "正在处理染色体/chr: ${SCAFF}"
    mkdir -p ${SCAFF}

    # 直接在本地执行 rfmix 命令，并将输出日志保存到对应文件夹
    rfmix \
        --query-file=query.vcf \
        --reference-file=ref.vcf \
        --sample-map=ref_sample.map \
        --genetic-map=genetic.map \
        --output-basename=${SCAFF}/${SCAFF} \
        --chromosome=${SCAFF} \
        --n-threads=20 \
        --random-seed=12345 > ${SCAFF}/${SCAFF}.log 2>&1

    echo "完成: ${SCAFF}"
done

echo "所有任务已运行完毕！"
