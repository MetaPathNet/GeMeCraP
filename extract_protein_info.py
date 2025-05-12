#!/usr/bin/env python3

import sys
import re

def extract_and_rename_genes(gff_file, output_file):
    """
    从GFF文件中提取基因信息并按原始基因ID顺序重命名为contig_1, contig_2等
    
    Args:
        gff_file: 输入的GFF文件路径
        output_file: 输出文件路径
    """
    genes = []
    
    # 读取GFF文件
    with open(gff_file, 'r') as f:
        for line in f:
            # 跳过注释行
            if line.startswith('#'):
                continue
            
            # 解析GFF行
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            # 只处理类型为"gene"的行
            if fields[2] == 'gene':
                contig = fields[0]  # 染色体/contig ID
                start = int(fields[3])  # 起始位置
                end = int(fields[4])    # 终止位置
                strand = fields[6]      # 正负链
                
                # 提取原始基因ID
                original_id = ""
                attributes = fields[8].split(';')
                for attr in attributes:
                    if attr.startswith('ID='):
                        original_id = attr.split('=')[1]
                        break
                
                # 保存基因信息
                genes.append({
                    'contig': contig,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'original_id': original_id
                })
    
    # 按原始基因ID排序
    # 假设原始ID格式为"MRSxxxxxx"，我们提取数字部分并按数字大小排序
    def get_id_number(gene):
        match = re.search(r'(\d+)', gene['original_id'])
        if match:
            return int(match.group(1))
        return 0
    
    genes.sort(key=get_id_number)
    
    # 将排序后的基因写入输出文件，重命名为contig_1, contig_2等
    with open(output_file, 'w') as out:
        # 写入表头
        out.write("new_geneID\toriginalID\tchromosome\tstart\tend\tstrand\n")
        
        # 写入基因信息
        for i, gene in enumerate(genes, 1):
            new_id = f"contig_{i}"
            out.write(f"{new_id}\t{gene['original_id']}\t{gene['contig']}\t"
                     f"{gene['start']}\t{gene['end']}\t{gene['strand']}\n")
    
    print(f"共处理了 {len(genes)} 个基因，并重命名为 contig_1 至 contig_{len(genes)}")
    print(f"结果已保存到 {output_file}")

if __name__ == "__main__":
    # 如果从命令行运行脚本并提供文件名
    if len(sys.argv) > 2:
        extract_and_rename_genes(sys.argv[1], sys.argv[2])
    else:
        # 使用默认文件名
        extract_and_rename_genes("input.gff", "renamed_genes.tsv")
