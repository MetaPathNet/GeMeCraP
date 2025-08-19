#!/usr/bin/env python3
import argparse
import json
import re

def load_adducts(adduct_file):
    """
    Load adduct masses from a file with lines like '+H 1.007825'
    """
    adducts = {}
    with open(adduct_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = re.split(r'\s+', line)
            adducts[parts[0]] = float(parts[1])
    return adducts

def neutral_mass(val, adducts):
    """
    Convert a string like '188.0707+H' to its neutral mass using adduct corrections.
    +H means subtract the mass of H; -H means add it.
    """
    match = re.match(r'([\d\.]+)([+\-][A-Za-z0-9]+)?', val)
    if not match:
        raise ValueError(f"Invalid value: {val}")
    mass = float(match.group(1))
    ad = match.group(2)
    if ad:
        try:
            ad_mass = adducts[ad]
        except KeyError:
            raise KeyError(f"Unknown adduct: {ad}")
        if ad.startswith('+'):
            mass -= ad_mass
        else:
            mass += ad_mass
    return mass

def ppm_diff(m1, m2):
    """Return the difference in ppm between two masses."""
    return abs(m1 - m2) / m1 * 1e6

def parse_clusters(cluster_file):
    """
    Parse clusters separated by lines of dashes (----). Returns a list of blocks,
    each block is a list of non-empty lines.
    """
    clusters = []
    block = []
    with open(cluster_file, 'r') as f:
        for line in f:
            stripped = line.rstrip()
            if stripped.startswith('-') and set(stripped) == {'-'}:
                if block:
                    clusters.append(block)
                    block = []
            else:
                if stripped:
                    block.append(stripped)
        if block:
            clusters.append(block)
    return clusters

def count_symbols(s):
    """Count occurrences of '+' or '-' in a string."""
    return s.count('+') + s.count('-')

def compute_block_masses(block, adducts):
    """
    For a block, extract lines with a 'source' key, compute their neutral masses,
    and count total '+'/'-' symbols.
    """
    pairs = [json.loads(line) for line in block if '"source"' in line]
    masses = []
    sym_count = 0
    for p in pairs:
        s, t = p['source'], p['target']
        m1 = neutral_mass(s, adducts)
        m2 = neutral_mass(t, adducts)
        masses.append((m1, m2))
        sym_count += count_symbols(s) + count_symbols(t)
    return masses, sym_count

def deduplicate(clusters, adducts, ppm_tol):
    """
    Deduplicate clusters: for each new block, compare its masses to existing unique blocks.
    If within ppm_tol on all pairs, consider redundant. Keep block with fewer symbols.
    """
    unique = []
    unique_data = []  # list of (masses, symbol_count)
    for block in clusters:
        masses, sym = compute_block_masses(block, adducts)
        found = False
        for idx, (umasses, usym) in enumerate(unique_data):
            if len(masses) == len(umasses) and all(
                ppm_diff(m1, u1) <= ppm_tol and ppm_diff(m2, u2) <= ppm_tol
                for (m1, m2), (u1, u2) in zip(masses, umasses)
            ):
                found = True
                # Keep the one with fewer '+'/'-'
                if sym < usym:
                    unique[idx] = block
                    unique_data[idx] = (masses, sym)
                break
        if not found:
            unique.append(block)
            unique_data.append((masses, sym))
    return unique

def write_clusters(clusters, out_file):
    """Write deduplicated clusters back to a file, separated by '----'."""
    with open(out_file, 'w') as f:
        for block in clusters:
            for line in block:
                f.write(line + '\n')
            f.write('----------------------------------------\n')

def main():
    parser = argparse.ArgumentParser(
        description='Deduplicate cluster entries in clust.txt based on 20 ppm tolerance.'
    )
    parser.add_argument('cluster_file', help='Path to cluster file (clust.txt)')
    parser.add_argument('adduct_file', help='Path to adduct mass file (adduct_file.txt)')
    parser.add_argument('-o', '--output', default='deduped_clusters.txt', help='Output file')
    parser.add_argument('--ppm', type=float, default=20.0, help='PPM tolerance (default: 20)')
    args = parser.parse_args()

    adducts = load_adducts(args.adduct_file)
    clusters = parse_clusters(args.cluster_file)
    unique = deduplicate(clusters, adducts, args.ppm)
    write_clusters(unique, args.output)

if __name__ == '__main__':
    main()

