#!/usr/bin/env python
import pysam  # 导入 pysam 库以处理 BAM/SAM 文件
import argparse  # 用于解析命令行参数
import logging  # 导入 logging 模块以支持日志记录
from collections import defaultdict  # 从 collections 模块导入 defaultdict，用于统计计数

# 配置日志系统，设置日志级别为 INFO 并定义日志格式
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('count reads')  # 创建一个名为 'count reads' 的日志记录器

def main():
    # 设置参数解析器
    parser = argparse.ArgumentParser(description='Count Small RNA Assigned to Each Transcripts')
    # 添加输入 BAM 文件参数（必需）
    parser.add_argument('--bam', '-b', type=str, required=True, help="Input bam file, should in transcript coordinate") 
    # 添加链特异性参数（默认为 "no"）
    parser.add_argument('--strandness', '-s', type=str, default="no", help="strandness", choices=["forward", "reverse", "no"])
    # 添加输出计数文件参数（必需）
    parser.add_argument('--counts', '-c', type=str, help="Output counts", required=True)
    # 添加统计信息输出文件参数（可选）
    parser.add_argument('--stats', type=str, help="Counting statistics")
    # 添加最小映射质量参数（默认为 0）
    parser.add_argument('--min-mapping-quality', '-m', default=0, type=int, help="Min mapping quality")
    
    # 解析命令行参数
    args = parser.parse_args()
    logger.info('read input BAM/SAM file: ' + args.bam)  # 记录读取的输入 BAM/SAM 文件
    
    # 打开 BAM 文件以进行读取
    sam = pysam.AlignmentFile(args.bam, "rb")
    stats = defaultdict(int)  # 初始化统计信息字典，用于统计不同类别的读取数
    counts = defaultdict(int)  # 初始化计数字典，用于计算每个转录本的读取数
    
    # 将链特异性字符串转换为整数（0: no, 1: forward, 2: reverse）
    strandness = {'no': 0, 'forward': 1, 'reverse': 2}.get(args.strandness, 0)
    
    # 遍历 BAM 文件中的每个读取
    for read in sam:
        if read.is_unmapped:  # 检查读取是否未对齐
            stats["unmapped"] += 1  # 增加未对齐读取的计数
            continue  # 继续下一条读取
        if read.mapping_quality < args.min_mapping_quality:  # 检查映射质量是否低于阈值
            stats["MAP-filtered"] += 1  # 增加映射质量过滤的计数
            continue  # 继续下一条读取
        # 检查链特异性的有效性
        if ((strandness == 1) and read.is_reverse) or ((strandness == 2) and (not read.is_reverse)):
            stats["invalid-strand"] += 1  # 增加无效链的计数
            continue  # 继续下一条读取
        
        # 计数有效的读取
        stats["counted"] += 1  # 增加计数的读取数
        counts[read.reference_name] += 1  # 根据参考名称增加对应转录本的计数
    
    # 将计数结果写入输出文件
    with open(args.counts, "w") as f:
        for sq in sam.header['SQ']:  # 遍历 SAM 文件中的所有参考序列
            seq_id = sq['SN']  # 获取序列 ID
            f.write(f"{seq_id}\t{counts[seq_id]}\n")  # 写入序列 ID 和其对应的读取数

    # 如果指定了统计信息输出文件，则将统计信息写入
    if args.stats is not None:
        with open(args.stats, "w") as f:
            f.write("tx_id\tcount\n")  # 写入标题行
            for category in ["counted", "unmapped", "MAP-filtered", "invalid-strand"]:
                f.write(f"{category}\t{stats[category]}\n")  # 写入每种类别的统计信息

# 主程序入口
if __name__ == "__main__":
    main()  # 调用主函数
    