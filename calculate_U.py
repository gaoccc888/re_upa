import allel
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
import functools
import logging
import os
import sys
from tqdm import tqdm

# --- 1. 基础参数设置 ---
vcf_file = '311_filtered.recode_noscaffold.vcf.gz' 
window_size = 50000  # 50 kb
step_size = 10000    # 10 kb step
outgroup_max_freq = 0.01  # 外类群最大频率 (1%)
target_min_freq = 0.10    # 目标群体最小频率 (10%)
source_fixed_freq = 1.0   # 源群体必须固定 (100%)
num_threads = 20          
log_filename = 'u_statistic_process.log'

# --- 2. 配置日志系统 ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler(log_filename, mode='w', encoding='utf-8'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def read_samples(file_path):
    if not os.path.exists(file_path):
        logger.error(f"找不到样本文件: {file_path}")
        return []
    with open(file_path, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def process_single_chrom(chrom, chroms, positions, mask_N2C, mask_SE2C, window_size, step_size):
    try:
        chrom_mask = (chroms == chrom)
        pos_chrom = positions[chrom_mask]
        u_sites_N2C_chrom = mask_N2C[chrom_mask].astype(int)
        u_sites_SE2C_chrom = mask_SE2C[chrom_mask].astype(int)
        
        # 滑窗计算
        counts_N2C, windows, _ = allel.windowed_statistic(
            pos_chrom, u_sites_N2C_chrom, statistic=np.sum, size=window_size, step=step_size
        )
        counts_SE2C, _, _ = allel.windowed_statistic(
            pos_chrom, u_sites_SE2C_chrom, statistic=np.sum, size=window_size, step=step_size
        )
        
        chrom_results = []
        for i in range(len(windows)):
            chrom_results.append({
                'CHROM': chrom,
                'START': windows[i][0],
                'END': windows[i][1],
                'U_North_to_Central': counts_N2C[i],
                'U_Southeast_to_Central': counts_SE2C[i]
            })
        return chrom_results
    except Exception as e:
        return f"Error in chrom {chrom}: {str(e)}"

def main():
    logger.info("程序开始运行 (严格固定模式: Source=100%, Outgroup<1%, Central>10%)...")
    
    north_samples = read_samples('north_pure.txt')
    central_samples = read_samples('central_all.txt')
    southeast_samples = read_samples('southeast_pure.txt')
    
    logger.info(f"正在读取 VCF 文件: {vcf_file}")
    try:
        callset = allel.read_vcf(vcf_file, fields=['calldata/GT', 'variants/CHROM', 'variants/POS', 'samples'])
        gt = allel.GenotypeArray(callset['calldata/GT'])
        chroms = callset['variants/CHROM']
        positions = callset['variants/POS']
        vcf_samples = list(callset['samples'])
    except Exception as e:
        logger.critical(f"读取 VCF 失败: {e}")
        return

    idx_north = [vcf_samples.index(s) for s in north_samples if s in vcf_samples]
    idx_central = [vcf_samples.index(s) for s in central_samples if s in vcf_samples]
    idx_southeast = [vcf_samples.index(s) for s in southeast_samples if s in vcf_samples]

    logger.info("开始计算各群体等位基因频率...")
    # 计算次等位基因频率 (假设 index 1 是我们要看的衍生等位基因)
    freq_north = gt.count_alleles(subpop=idx_north).to_frequencies()[:, 1]
    freq_central = gt.count_alleles(subpop=idx_central).to_frequencies()[:, 1]
    freq_southeast = gt.count_alleles(subpop=idx_southeast).to_frequencies()[:, 1]

    # --- 核心修改部分：严格判定逻辑 ---
    logger.info("正在执行严格位点过滤...")
    
    # 逻辑 1: North -> Central (North固定, SE几乎没有, Central有频率)
    mask_N2C = (freq_north == source_fixed_freq) & \
               (freq_southeast < outgroup_max_freq) & \
               (freq_central > target_min_freq)
    
    # 逻辑 2: Southeast -> Central (SE固定, North几乎没有, Central有频率)
    mask_SE2C = (freq_southeast == source_fixed_freq) & \
                (freq_north < outgroup_max_freq) & \
                (freq_central > target_min_freq)

    unique_chroms = np.unique(chroms)
    logger.info(f"检测到 {len(unique_chroms)} 条染色体，准备并行处理。")

    worker_func = functools.partial(
        process_single_chrom, 
        chroms=chroms, 
        positions=positions, 
        mask_N2C=mask_N2C, 
        mask_SE2C=mask_SE2C, 
        window_size=window_size,
        step_size=step_size
    )

    all_window_results = []
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = list(tqdm(executor.map(worker_func, unique_chroms), 
                           total=len(unique_chroms), 
                           desc="染色体处理中"))
        
        for res in futures:
            if isinstance(res, list):
                all_window_results.extend(res)
            else:
                logger.error(res)

    df_results = pd.DataFrame(all_window_results)
    df_results = df_results.sort_values(['CHROM', 'START']).reset_index(drop=True)
    
    output_file = 'U_statistic_strict_fixed_50k_10k.csv'
    df_results.to_csv(output_file, index=False)
    logger.info(f"计算完成! 结果已存入 {output_file}")
    logger.info(f"N2C 总有效位点数: {np.sum(mask_N2C)}, SE2C 总有效位点数: {np.sum(mask_SE2C)}")

if __name__ == '__main__':
    main()
