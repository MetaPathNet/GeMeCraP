import argparse
import json
import re
import itertools
import math
from pathlib import Path
from typing import Dict, List, Tuple, Set


def parse_cluster_file(path: Path) -> List[Dict]:
    """
    Parse the cluster file into a list of cluster entries.
    Each entry is a dict with keys: 'entries' (list of JSON dicts) and 'genes' (str).
    """
    clusters = []
    raw = path.read_text().splitlines()
    sep_pattern = re.compile(r"^-{5,}")
    current = []
    for line in raw:
        if sep_pattern.match(line):
            if current:
                clusters.append(current)
                current = []
        elif line.strip():
            current.append(line.strip())
    if current:
        clusters.append(current)

    parsed = []
    for block in clusters:
        entry = {'entries': [], 'genes': ''}
        for line in block:
            if line.startswith('{') and 'cluster_genes' not in line:
                entry['entries'].append(json.loads(line))
            elif 'cluster_genes' in line:
                data = json.loads(line)
                entry['genes'] = data.get('cluster_genes', '')
        parsed.append(entry)
    return parsed


def load_group_mapping(path: Path) -> Dict[str, List[str]]:
    """
    Load central_group.txt or mz_group.txt into a mapping from cluster_metabolite -> list of pre_metabolites
    """
    mapping: Dict[str, List[str]] = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = re.split(r"\s+", line.strip(), maxsplit=1)
        key = parts[0]
        values = [v.strip() for v in parts[1].split(',')] if len(parts) > 1 else []
        mapping[key] = values
    return mapping


def load_msms(path: Path) -> Dict[str, List[Tuple[float, float]]]:
    """
    Parse the MSMS.msp file into a dict mapping cleaned_name -> list of (mz, intensity)
    """
    content = path.read_text()
    # split entries by blank line
    entries = re.split(r"\n\s*\n", content.strip())
    msms_data: Dict[str, List[Tuple[float, float]]] = {}
    for entry in entries:
        lines = entry.splitlines()
        name_line = next((l for l in lines if l.startswith('Name:')), None)
        if not name_line:
            continue
        match = re.search(r"\(([^)]+)\)", name_line)
        if not match:
            continue
        raw = match.group(1)
        cleaned = clean_metabolite(raw)
        # extract numeric peak lines
        peaks: List[Tuple[float, float]] = []
        for l in lines:
            if re.match(r"^[0-9]+\.?[0-9]*\s+", l):
                parts = l.split()
                try:
                    mz_val = float(parts[0])
                    intensity = float(parts[1])
                    peaks.append((mz_val, intensity))
                except (ValueError, IndexError):
                    continue
        msms_data[cleaned] = peaks
    return msms_data


def get_top_peaks(peaks: List[Tuple[float, float]]) -> List[float]:
    """
    Return the list of m/z values of the top 50% peaks by intensity, preserving floats.
    """
    if not peaks:
        return []
    sorted_peaks = sorted(peaks, key=lambda x: x[1], reverse=True)
    cutoff = math.ceil(len(sorted_peaks) / 2)
    return [mz for mz, _ in sorted_peaks[:cutoff]]


def clean_metabolite(s: str) -> str:
    """
    Remove any trailing letters or symbols, keeping only digits, underscores, and dots.
    """
    return re.sub(r"[^0-9_.]+$", "", s)


def find_pre_metabolites(node: str, central: Dict[str, List[str]], mz: Dict[str, List[str]]) -> List[str]:
    key = re.match(r"^[0-9.]+", node)
    key_str = key.group(0) if key else node
    if node == key_str:
        return central.get(key_str, [])
    return mz.get(key_str, [])


def ppm_match(list1: List[float], list2: List[float], ppm: float = 20.0) -> Set[float]:
    """
    Return set of m/z from list1 that have at least one match in list2 within given ppm tolerance.
    """
    matches: Set[float] = set()
    for mz1 in list1:
        for mz2 in list2:
            if abs(mz1 - mz2) <= mz1 * ppm / 1e6:
                matches.add(mz1)
                break
    return matches


def find_common_fragments_ppm(cleaned_path: List[str], msms: Dict[str, List[Tuple[float, float]]]) -> List[Set[float]]:
    """
    For each sliding window of 3 cleaned metabolites, find common peaks within 20 ppm tolerance.
    """
    # get top 50% peaks for each metabolite
    top_peaks = []
    for met in cleaned_path:
        peaks = msms.get(met, [])
        top = get_top_peaks(peaks)
        top_peaks.append(top)
    # sliding window of 3
    commons = []
    for i in range(len(top_peaks) - 2):
        # match between first and second
        match12 = ppm_match(top_peaks[i], top_peaks[i+1])
        # match those results with third
        match123 = ppm_match(list(match12), top_peaks[i+2])
        commons.append(match123)
    return commons


def main():
    parser = argparse.ArgumentParser(
        description="Analyze cluster paths and MSMS fragment overlaps within 20 ppm."
    )
    parser.add_argument("--cluster", type=Path, required=True)
    parser.add_argument("--central", type=Path, required=True)
    parser.add_argument("--mz", type=Path, required=True)
    parser.add_argument("--msms", type=Path, required=True)
    args = parser.parse_args()

    clusters = parse_cluster_file(args.cluster)
    central_map = load_group_mapping(args.central)
    mz_map = load_group_mapping(args.mz)
    msms_data = load_msms(args.msms)

    for cl in clusters:
        entries = cl['entries']
        genes = cl['genes']
        metabolites = [entries[0]['source']] + [e['target'] for e in entries]
        print(f"cluster_genes: {genes}\nshared fragment <= 20ppm")
        print("->".join(metabolites))
        pre_lists = [find_pre_metabolites(m, central_map, mz_map) for m in metabolites]
        multi = [m for m, lst in zip(metabolites, pre_lists) if len(lst) > 1]
        if multi:
            print(f"Note: The grouped metabolites below correspond to multiple actual metabolites with retention times: {multi}")
        # each real path combination
        for rp in itertools.product(*pre_lists):
            cleaned_rp = [clean_metabolite(m) for m in rp]
            commons = find_common_fragments_ppm(cleaned_rp, msms_data)
            segments = []
            for i, com in enumerate(commons):
                trio = cleaned_rp[i:i+3]
                segments.append(f"{'->'.join(trio)} : {sorted(com)}")
            print("; ".join(segments))
        print()

if __name__ == "__main__":
    main()

