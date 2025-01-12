import pandas as pd

# 读取 metadata.txt 文件
metadata_file = '/data/2024-bioinfo-shared/lichengxi/exRNA-short/metadata.txt'
metadata = pd.read_csv(metadata_file, sep='\t')  
script_path = '/data/2024-bioinfo-shared/lichengxi/count_transcripts.py'

# 遍历 metadata 数据，构建命令
commands = []
for index, row in metadata.iterrows():
    sample_id = row['sample id']

    # 假设 BAM 文件的名称由 sample_id 决定
    bam_file_mi = f'/data/2024-bioinfo-shared/lichengxi/exRNA-short/bam/{sample_id}/miRNA.bam'
    bam_file_pi = f'/data/2024-bioinfo-shared/lichengxi/exRNA-short/bam/{sample_id}/piRNA.bam'
    counts_file_mi = f'/data/2024-bioinfo-shared/lichengxi/output/counts_miRNA/{sample_id}'
    counts_file_pi = f'/data/2024-bioinfo-shared/lichengxi/output/counts_piRNA/{sample_id}'
    stats_file_mi = f'/data/2024-bioinfo-shared/lichengxi/output/stats_miRNA/{sample_id}.txt'
    stats_file_pi = f'/data/2024-bioinfo-shared/lichengxi/output/stats_piRNA/{sample_id}.txt'
    

    command_mi = (
        f'python {script_path} --bam {bam_file_mi} --counts {counts_file_mi} --stats {stats_file_mi} -s forward'
    )
    command_pi = (
        f'python {script_path} --bam {bam_file_pi} --counts {counts_file_pi} --stats {stats_file_pi} -s forward'
    )

    commands.append(command_mi)
    commands.append(command_pi)


# 输出生成的命令
for cmd in commands:
    print(cmd)
tool_path = '/data/2022-bioinfo-shared/softwares/miniconda3/envs/python-env/bin'

# 将命令写入bash脚本
with open('/data/2024-bioinfo-shared/lichengxi/run/count_trancripts.sh', 'w') as fh:
    fh.write('#!/bin/bash\n#SBATCH -J counttranscript\n#SBATCH -p CN_BIOT\n#SBATCH --nodes=1\n#SBATCH --ntasks=4\n#SBATCH --output=%j.out\n#SBATCH --error=%j.err')
    fh.write('\n\n')
    fh.write('# Get the software\n')
    fh.write(f"export PATH={tool_path}:$PATH")
    fh.write('\n\n')
    fh.write('# run code.\n')
    for cmd in commands:
        fh.write(cmd)
        fh.write('\n')
