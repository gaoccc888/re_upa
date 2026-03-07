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




##################################Global Ancestry Painting (GAP)绘图

python /home/gao/workspace/RFmix/Output_GAP/MergeMultiSample.py

python /home/gao/workspace/RFmix/RFMIX2-Pipeline-to-plot-main/GAP/BedToGap.py \
  --input /home/gao/workspace/RFmix/Output_GAP/Sample_Merged_All.bed \
  --out /home/gao/workspace/RFmix/Output_GAP/Final_Plot_Input.bed
python /home/gao/workspace/RFmix/RFMIX2-Pipeline-to-plot-main/GAP/GAP.py \
  --input /home/gao/workspace/RFmix/Output_GAP/Final_Plot_Input.bed \
  --output /home/gao/workspace/RFmix/Output_GAP/Final_Ancestry_Plot.pdf

##################################Local Ancestry Painting (LAP) 绘图、
python /home/gao/workspace/RFmix/Output_GAP/split_msp.py
# ==========================================
# 1. 定义基础路径（请确保这些路径在你的机器上是正确的）
# ==========================================
BASE_PATH="/home/gao/workspace/RFmix/RFMIX2-Pipeline-to-plot-main/LAP"
SCRIPTS_DIR="/home/gao/workspace/RFmix/RFMIX2-Pipeline-to-plot-main/LAP"
SPLIT_DIR="/home/gao/workspace/RFmix/split_msp_files"
OUT_DIR="/home/gao/workspace/RFmix/Output_LAP"

# 创建输出目录
mkdir -p $OUT_DIR

# 确保你在 LAP 根目录下，因为 LAP.py 需要在当前目录读取 hg38.svg
cd $BASE_PATH

# ==========================================
# 2. 开始循环处理每个个体
# ==========================================
# 获取样本列表（提取前缀）
SAMPLES=$(ls $SPLIT_DIR | cut -d'_' -f1 | sort | uniq)

for s in $SAMPLES; do
    echo "------------------------------------------"
    echo "🚀 正在处理个体: $s"
    echo "------------------------------------------"

    # Step 1: 合并染色体 (生成 hap1.bed 和 hap2.bed)
    python "$SCRIPTS_DIR/RFMIX2ToBed.py" \
        --prefix "$SPLIT_DIR/$s" \
        --chr {1..14} \
        --output "$OUT_DIR/"

    # Step 2: 合并单倍体
    python "$SCRIPTS_DIR/BedToLAP.py" \
        --bed1 "$OUT_DIR/${s}_hap1.bed" \
        --bed2 "$OUT_DIR/${s}_hap2.bed" \
        --out "$OUT_DIR/${s}_plot_input.bed"

    # 关键点：删除 'Chr' 前缀（防止 int() 转换报错）
    sed -i 's/Chr//g' "$OUT_DIR/${s}_plot_input.bed"

    # Step 3: 绘制局部祖源图 (LAP)
    # 此时 LAP.py 会调用 rsvg-convert 将 SVG 转为 PDF
    python "LAP.py" \
        -I "$OUT_DIR/${s}_plot_input.bed" \
        -B hg38 \
        -O "$OUT_DIR/${s}_Karyogram.pdf"

    if [ $? -eq 0 ]; then
        echo "✅ 个体 $s 绘制完成，结果已存入 $OUT_DIR"
    else
        echo "❌ 个体 $s 绘制失败，请检查上方报错"
    fi
done
###全部LAI
python /home/gao/workspace/RFmix/Output_GAP/plot_population_summary.py   
