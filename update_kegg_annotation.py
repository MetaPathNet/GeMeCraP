#!/usr/bin/env python3

import sys

def update_kegg_annotation(renamed_genes_file, kegg_file, output_file):
    """
    根据重命名后的基因信息，更新KEGG注释文件中的基因名称
    
    Args:
        renamed_genes_file: 重命名后的基因信息文件路径
        kegg_file: 原KEGG注释文件路径
        output_file: 更新后的KEGG注释文件输出路径
    """
    # 创建原始基因ID到新基因ID的映射
    id_mapping = {}
    
    # 读取重命名后的基因信息文件
    with open(renamed_genes_file, 'r') as f:
        # 跳过表头
        header = f.readline()
        
        # 读取每一行信息
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                new_id = fields[0]  # contig_1, contig_2, ...
                original_id = fields[1]  # MRS000001, MRS000002, ...
                
                # 建立映射关系
                id_mapping[original_id] = new_id
    
    # 读取KEGG注释文件并更新基因ID
    updated_annotations = []
    not_found = []
    
    with open(kegg_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                original_id = fields[0]  # 原始基因ID
                kegg_id = fields[1]      # KEGG ID
                
                # 查找对应的新基因ID
                if original_id in id_mapping:
                    new_id = id_mapping[original_id]
                    updated_annotations.append((new_id, kegg_id))
                else:
                    # 记录未找到对应关系的基因
                    not_found.append(original_id)
                    # 保留原始ID
                    updated_annotations.append((original_id, kegg_id))
    
    # 将更新后的注释写入输出文件
    with open(output_file, 'w') as out:
        for gene_id, kegg_id in updated_annotations:
            out.write(f"{gene_id}\t{kegg_id}\n")
    
    # 打印处理结果
    print(f"注释文件更新完成，共处理了 {len(updated_annotations)} 条注释记录")
    if not_found:
        print(f"警告: 有 {len(not_found)} 个基因在重命名文件中未找到对应关系")
        print("未找到对应关系的基因: " + ", ".join(not_found[:5]) + 
              ("..." if len(not_found) > 5 else ""))
    print(f"更新后的注释已保存到 {output_file}")

if __name__ == "__main__":
    if len(sys.argv) > 3:
        update_kegg_annotation(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        # 使用默认文件名
        update_kegg_annotation("renamed_genes.tsv", "kegg_annotation.txt", "updated_kegg_annotation.txt")
