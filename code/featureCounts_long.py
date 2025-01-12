import pandas as pd

# 读取 metadata.txt 文件
metadata_file = '/data/2024-bioinfo-shared/lichengxi/exRNA-long/metadata.txt'
metadata = pd.read_csv(metadata_file, sep='\t')  

# 遍历 metadata 数据，构建 featureCounts 命令
commands = []
for index, row in metadata.iterrows():
    sample_id = row['sample_id']
    library = row['library']

    # 假设 BAM 文件的名称由 sample_id 决定
    bam_file = f'/data/2024-bioinfo-shared/lichengxi/exRNA-long/bam/{sample_id}.bam'  # BAM文件路径
    output_prefix = f'/data/2024-bioinfo-shared/lichengxi/output/counts/{sample_id}'  # 输出路径

    if library == 'reverse':
        strand_option = '-s 2'
    elif library == 'forward':
        strand_option = '-s 1'
    else:
        print(f"未知的 library 值: {library}")
        continue

    command = (
        f"featureCounts --countReadPairs -O -M {strand_option} -p -t exon -g gene_id "
        f"-a /data/2024-bioinfo-shared/lichengxi/genenome/gencode.v38.annotation.gff3 -o {output_prefix} {bam_file}"
    )

    commands.append(command)

# 输出生成的命令
for cmd in commands:
    print(cmd)

# 将命令写入bash脚本
with open('/data/2024-bioinfo-shared/lichengxi/run/featureCounts_long.sh', 'w') as fh:
    fh.write('#!/bin/bash\n#SBATCH -J featureCounts_long\n#SBATCH -p CN_BIOT\n#SBATCH --nodes=1\n#SBATCH --ntasks=4\n#SBATCH --output=%j.out\n#SBATCH --error=%j.err')
    fh.write('\n\n')
    fh.write('# Get the software\n')
    fh.write(f"export PATH=/data/2022-bioinfo-shared/softwares/subread/subread-2.0.3-Linux-x86_64/bin:$PATH")
    fh.write('\n\n')
    fh.write('# run code.\n')
    for cmd in commands:
        fh.write(cmd)
        fh.write('\n')
        
