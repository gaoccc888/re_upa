# 1. 进入工作目录
cd /home/gao/workspace/ancestral_allele/Radical

# 2. 复制矩阵文件 (如果还没复制)
cp /home/gao/workspace/ancestral_allele/Radical/annovar/example/grantham.matrix ./

# 3. 格式转换 (VCF -> Annovar)
# 注意：如果输入的是标准VCF，用 vcf4 可能比 vcf4old 更稳妥，但如果确定是旧格式则保持不变
/home/gao/workspace/ancestral_allele/Radical/annovar/convert2annovar.pl \
    --format vcf4old \
    /home/gao/workspace/ancestral_allele/314_filtered_noscaffold_renamed.easySFS.recode_polarized.vcf \
    --outfile 314_filtered_noscaffold_renamed.easySFS.recode_polarized.avinput

NOTICE: Read 3265524 lines and wrote 3167932 different variants at 3265478 genomic positions (3265478 SNPs and 0 indels)
NOTICE: Among 3265478 different variants at 3265478 positions, 869095 are heterozygotes, 2298837 are homozygotes
NOTICE: Among 3265478 SNPs, 2113149 are transitions, 1152329 are transversions (ratio=1.83)


# 4. 运行注释 (核心步骤)
# 数据库目录指向你存放 SFZ.A.onlychr_refGene.txt 的地方
/home/gao/workspace/ancestral_allele/Radical/annovar/annotate_variation.pl \
    -buildver Upa \
    -outfile Upa \
    --aamatrixfile grantham.matrix \
    314_filtered_noscaffold_renamed.easySFS.recode_polarized.avinput \
    /home/gao/workspace/ancestral_allele/Radical/annovar/

# 5. 提取所有非同义突变
grep 'nonsynonymous' Upa.exonic_variant_function > Upa.grantham.nonsynonymous_variants.txt

# 6. 筛选 Grantham Score > 150 的激进突变
# 解释：awk 根据 "AAMatrix=" 切分，取后面的数字判断是否 > 150
awk -F'AAMatrix=' '{split($2, a, ","); if(a[1] > 150) print $0}' Upa.grantham.nonsynonymous_variants.txt > Upa.grantham.radical_sites_list.txt

# 7. 统计并展示结果数量
echo "--------------------------------"
echo "非同义突变总数:"
wc -l  Upa.grantham.nonsynonymous_variants.txt
echo "Radical (激进) 突变数 (>150):"
wc -l Upa.grantham.radical_sites_list.txt
echo "--------------------------------"

#非同义突变总数:
#121287 Upa.grantham.nonsynonymous_variants.txt
#Radical (激进) 突变数 (>150):
#7117 Upa.grantham.radical_sites_list.txt


# 8.提取 Radical突变位点
# 提取 Chr 和 Pos，用制表符分隔
awk '{print $5"\t"$6}' Upa.grantham.radical_sites_list.txt > Upa.grantham.radical_positions.txt

# 检查一下格式 (应该是两列：Chr01  18764)
head Upa.grantham.radical_positions.txt

vcftools \
    --vcf /home/gao/workspace/ancestral_allele/314_filtered_noscaffold_renamed.easySFS.recode_polarized.vcf \
    --positions Upa.grantham.radical_positions.txt \
    --recode --recode-INFO-all \
    --out Upa_polarized.annovar_Radical

## After filtering, kept 314 out of 314 Individuals
## Outputting VCF file...
## After filtering, kept 7117 out of a possible 3265478 Sites
## Run Time = 28.00 seconds

