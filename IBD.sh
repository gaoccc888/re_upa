########################提取 IBD 片段 (Refined IBD)

# 首先查看你的 VCF 里到底有哪些染色体名称
zgrep -v "^#" /home/gao/workspace/RFmix/311_phased_all_genome.vcf.gz | cut -f1 | uniq > chrom_list.txt

# 循环运行每一条染色体
for chr in $(cat chrom_list.txt); do
    echo "Processing chromosome: $chr"
    java -Xmx100g -jar refined-ibd.17Jan20.102.jar \
        gt=/home/gao/workspace/RFmix/311_phased_all_genome.vcf.gz \
        chrom=$chr \
        out=Ulmus_Refined_IBD_$chr \
        window=40.0 \
        lod=2.0 \
        nthreads=20
done

# 运行完成后，将结果合并（可选）
zcat Ulmus_Refined_IBD_*.ibd.gz | gzip > Ulmus_Refined_IBD_all.ibd.gz



