#!/usr/bin/env python3
import sys
import re
from collections import defaultdict

def parse_metabolite(line):
    """解析代谢物字符串，返回retention time和mass"""
    rt, mass_str = line.strip().split('_')
    # 使用正则表达式提取数字部分作为mass
    mass = float(re.match(r'(\d+\.?\d*)', mass_str).group(1))
    return rt, mass, line.strip()

def calculate_ppm(mass1, mass2):
    """计算两个质量之间的ppm差异"""
    return abs(mass1 - mass2) / min(mass1, mass2) * 1e6

def group_metabolites(filename, ppm_threshold=10):
    """将代谢物按照mass分组"""
    # 读取所有代谢物
    metabolites = []
    with open(filename) as f:
        for line in f:
            if line.strip():
                rt, mass, original = parse_metabolite(line)
                metabolites.append((mass, original))
    
    # 按质量排序
    metabolites.sort()
    
    # 分组
    groups = defaultdict(list)
    current_group_mass = None
    
    for mass, original in metabolites:
        if current_group_mass is None:
            current_group_mass = mass
            groups[mass].append(original)
        else:
            if calculate_ppm(mass, current_group_mass) <= ppm_threshold:
                groups[current_group_mass].append(original)
            else:
                current_group_mass = mass
                groups[mass].append(original)
    
    return groups

def main():
    if len(sys.argv) != 2:
        print("example: python group_metabolites.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    groups = group_metabolites(input_file)
    
    # 输出结果
    for mass, metabolites in sorted(groups.items()):
        print(f"{mass}\t{','.join(metabolites)}")

if __name__ == "__main__":
    main()
