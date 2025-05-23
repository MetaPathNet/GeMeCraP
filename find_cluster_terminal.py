import argparse
import json
import pandas as pd
import ast
from typing import List, Dict, Set, Any, Tuple, Optional
from collections import defaultdict
import itertools
import logging

# Set logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_network(network_file: str) -> List[List[Dict[str, Any]]]:
    """Parse the network file with pathways separated by '----------------------------------------'."""
    pathways = []
    pathway = []
    separator = "----------------------------------------"

    try:
        with open(network_file, 'r', encoding='utf-8') as file:
            for line in file:
                line = line.strip()
                if line == separator:
                    if pathway:
                        pathways.append(pathway)
                        pathway = []
                elif line.startswith("{") and line.endswith("}"):
                    try:
                        line_dict = ast.literal_eval(line)
                        pathway.append(line_dict)
                    except (ValueError, SyntaxError) as e:
                        logging.warning(f"Skipping invalid line: {line}\nError: {e}")
            if pathway:
                pathways.append(pathway)
    except FileNotFoundError:
        logging.error(f"Network file not found: {network_file}")
        exit(1)
    except Exception as e:
        logging.error(f"Error reading network file: {e}")
        exit(1)
    
    logging.info(f"Parsed {len(pathways)} pathways from network file.")
    return pathways


def load_reaction_mapping(reaction_file: str) -> Dict[str, List[str]]:
    """
    Load reaction to Orthology (K numbers) mapping from a tab-separated file.
    Now supports duplicate ENTRY rows by appending K numbers.
    """
    reaction_to_k: Dict[str, List[str]] = defaultdict(list)
    try:
        reaction_df = pd.read_csv(reaction_file, sep="\t", dtype=str)
    except FileNotFoundError:
        logging.error(f"Reaction mapping file not found: {reaction_file}")
        exit(1)
    except Exception as e:
        logging.error(f"Error reading reaction mapping file: {e}")
        exit(1)

    for _, row in reaction_df.iterrows():
        entry = row.get("ENTRY")
        orthology = row.get("Orthology") or ''
        if pd.notnull(entry) and pd.notnull(orthology):
            ks = [k.strip() for k in orthology.split(",") if k.strip()]
            reaction_to_k[entry].extend(ks)

    for entry, ks in reaction_to_k.items():
        reaction_to_k[entry] = list(dict.fromkeys(ks))

    logging.info(f"Loaded reaction to K mapping for {len(reaction_to_k)} entries.")
    return dict(reaction_to_k)


def load_k_to_gene_mapping(kegg_file: str) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
    """
    Load KEGG file mapping K numbers to genes and capture annotation detail.
    Expected columns: Gene, K_number, Annotation
    """
    k_to_gene: Dict[str, Set[str]] = defaultdict(set)
    gene_annotations: Dict[str, str] = {}

    try:
        with open(kegg_file, 'r', encoding='utf-8') as file:
            for line in file:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    gene, k_number, annotation = parts[0].strip(), parts[1].strip(), parts[-1].strip()
                    gene_annotations[gene] = annotation
                    if k_number.startswith("K"):
                        k_to_gene[k_number].add(gene)
                else:
                    logging.debug(f"Skipping malformed line: {line}")
    except FileNotFoundError:
        logging.error(f"KEGG mapping file not found: {kegg_file}")
        exit(1)
    except Exception as e:
        logging.error(f"Error reading KEGG mapping file: {e}")
        exit(1)

    logging.info(f"Loaded K to gene mapping for {len(k_to_gene)} K numbers.")
    return k_to_gene, gene_annotations


def load_gene_list(gene_list_file: str) -> Set[str]:
    """Load transcriptome differential genes from a file."""
    try:
        with open(gene_list_file, 'r', encoding='utf-8') as file:
            genes = set(line.strip() for line in file)
    except FileNotFoundError:
        logging.error(f"Gene list file not found: {gene_list_file}")
        exit(1)
    except Exception as e:
        logging.error(f"Error reading gene list file: {e}")
        exit(1)
    
    logging.info(f"Loaded {len(genes)} differential genes.")
    return genes


def validate_cluster(
    cluster: List[str],
    differential_count: int,
    cluster_edges: Set[int],
    total_edges: int
) -> Optional[Dict[str, Any]]:
    """
    Validate a cluster based on size, differential genes, and full coverage of edges.
    """
    size = len(cluster)
    covered = len(cluster_edges) == total_edges
    if 2 <= size <= 3 and differential_count >= 1 and covered:
        return {"genes": cluster.copy()}
    if size > 3 and differential_count >= 2 and covered:
        return {"genes": cluster.copy()}
    return None


def expand_cluster_with_gaps(cluster_genes: List[str]) -> List[str]:
    """
    Expand a cluster by filling in any missing numeric suffix genes between the
    minimum and maximum suffix for the same prefix.
    """
    if not cluster_genes:
        return []
    prefixes = []
    suffixes = []
    for gene in cluster_genes:
        if "_" in gene:
            prefix, suffix_str = gene.rsplit("_", 1)
            try:
                suffix = int(suffix_str)
            except ValueError:
                continue
            prefixes.append(prefix)
            suffixes.append(suffix)
    if not suffixes:
        return cluster_genes.copy()
    prefix = prefixes[0]
    min_suf, max_suf = min(suffixes), max(suffixes)
    expanded = [f"{prefix}_{i}" for i in range(min_suf, max_suf + 1)]
    # Sort and return unique
    return sort_genes(set(expanded))


def identify_clusters(
    sorted_genes: List[str],
    differential_genes: Set[str],
    gene_to_edges: Dict[str, Set[int]],
    total_edges: int,
    max_gap: int,
    gap_count: int
) -> List[Dict[str, Any]]:
    """
    Identify clusters based on adjacency, then expand gaps and validate.
    """
    clusters = []
    current_cluster: List[str] = []
    differential_count = 0
    current_edges: Set[int] = set()
    current_gap = 0

    def process_current():
        nonlocal clusters
        if not current_cluster:
            return
        # Expand cluster to include gap genes
        expanded = expand_cluster_with_gaps(current_cluster)
        # Recalculate differential count and edges
        diff_count_expanded = sum(1 for gene in expanded if gene in differential_genes)
        edges_expanded: Set[int] = set()
        for gene in expanded:
            edges_expanded.update(gene_to_edges.get(gene, set()))
        # Validate expanded cluster
        result = validate_cluster(expanded, diff_count_expanded, edges_expanded, total_edges)
        if result:
            clusters.append(result)

    for gene in sorted_genes:
        if not current_cluster:
            current_cluster = [gene]
            differential_count = 1 if gene in differential_genes else 0
            current_edges = gene_to_edges.get(gene, set()).copy()
            current_gap = 0
            continue

        prev = current_cluster[-1]
        prefix_prev, suffix_prev = prev.rsplit("_", 1) if "_" in prev else (prev, '0')
        prefix_cur, suffix_cur = gene.rsplit("_", 1) if "_" in gene else (gene, '0')

        same_prefix = (prefix_prev == prefix_cur)
        gap = 0
        if same_prefix:
            try:
                gap = int(suffix_cur) - int(suffix_prev)
            except ValueError:
                gap = 0

        if same_prefix and gap == 1:
            current_cluster.append(gene)
            if gene in differential_genes:
                differential_count += 1
            current_edges.update(gene_to_edges.get(gene, set()))
        elif same_prefix and 2 <= gap <= max_gap and (gap_count == 1 or current_gap < 1):
            current_cluster.append(gene)
            if gene in differential_genes:
                differential_count += 1
            current_edges.update(gene_to_edges.get(gene, set()))
            current_gap += 1
        else:
            # Process and reset
            process_current()
            current_cluster = [gene]
            differential_count = 1 if gene in differential_genes else 0
            current_edges = gene_to_edges.get(gene, set()).copy()
            current_gap = 0

    # Final cluster
    process_current()
    logging.info(f"Identified {len(clusters)} clusters.")
    return clusters


def process_pathways(
    pathways: List[List[Dict[str, Any]]],
    reaction_to_k: Dict[str, List[str]],
    k_to_gene: Dict[str, Set[str]],
    gene_annotations: Dict[str, str],
    differential_genes: Set[str],
    max_gap: int,
    gap_count: int
) -> List[Dict[str, Any]]:
    """
    Process pathways to include activator genes and identify valid clusters.
    """
    results = []
    all_k_in_reactions = set(itertools.chain.from_iterable(reaction_to_k.values()))

    for idx, pathway in enumerate(pathways):
        gene_to_edges: Dict[str, Set[int]] = defaultdict(set)
        pathway_genes: Set[str] = set()

        for e_idx, edge in enumerate(pathway):
            for reaction in edge.get("diff", []):
                for k in reaction_to_k.get(reaction, []):
                    for gene in k_to_gene.get(k, []):
                        pathway_genes.add(gene)
                        gene_to_edges[gene].add(e_idx)

        total_edges = len(pathway)

        for k, genes in k_to_gene.items():
            if k not in all_k_in_reactions:
                for gene in genes:
                    ann = gene_annotations.get(gene, '').lower()
                    if 'activator' in ann or 'activating' in ann:
                        if gene not in pathway_genes:
                            pathway_genes.add(gene)
                            gene_to_edges[gene] = set()

        if pathway_genes & differential_genes:
            logging.info(f"Processing pathway {idx} with {len(pathway_genes)} genes.")
            sorted_genes = sort_genes(pathway_genes)
            clusters = identify_clusters(
                sorted_genes,
                differential_genes,
                gene_to_edges,
                total_edges,
                max_gap,
                gap_count
            )
            if clusters:
                results.append({"pathway_index": idx, "pathway": pathway, "clusters": clusters})

    logging.info(f"Total pathways processed: {len(results)}")
    return results


def sort_genes(genes: Set[str]) -> List[str]:
    """Sort genes by prefix and numeric suffix."""
    def key_fn(g: str):
        parts = g.rsplit("_", 1)
        prefix = parts[0]
        try:
            suffix = int(parts[1])
        except (IndexError, ValueError):
            suffix = 0
        return prefix, suffix
    return sorted(genes, key=key_fn)


def write_output(output_file: str, processed_pathways: List[Dict[str, Any]]) -> None:
    separator = "----------------------------------------"
    try:
        with open(output_file, 'w', encoding='utf-8') as file:
            for info in processed_pathways:
                for edge in info["pathway"]:
                    file.write(json.dumps({"source": edge.get("source"),
                                           "target": edge.get("target"),
                                           "diff": edge.get("diff")}, ensure_ascii=False) + "\n")
                for cl in info["clusters"]:
                    file.write(json.dumps({"cluster_genes": ", ".join(cl["genes"])}, ensure_ascii=False) + "\n")
                file.write(separator + "\n")
        logging.info(f"Output written to {output_file}")
    except Exception as e:
        logging.error(f"Error writing output file: {e}")
        exit(1)


def main():
    parser = argparse.ArgumentParser(description="Process metabolic network pathways and extract related gene clusters.")
    parser.add_argument("--network", required=True, help="Path to the metabolic network file.")
    parser.add_argument("--reaction", required=True, help="Path to the file containing all reactions.")
    parser.add_argument("--kegg", required=True, help="Path to the KEGG annotation file.")
    parser.add_argument("--gene_list", required=True, help="Path to the transcriptome differential gene list.")
    parser.add_argument("--output", required=True, help="Path to save the output results.")
    parser.add_argument("--max_gap", type=int, default=2, help="Maximum allowable gap between adjacent genes in a cluster."
                             "For example, max_gap=2 allows up to one gene between adjacent genes.default max_gap = 2.")
    parser.add_argument("--gap_count", type=int, choices=[0,1], default=0, help="Set to 1 to allow multiple gaps within a cluster or 0 to allow at most one gap.Default = 0.")
    args = parser.parse_args()

    pathways = parse_network(args.network)
    reaction_to_k = load_reaction_mapping(args.reaction)
    k_to_gene, gene_annotations = load_k_to_gene_mapping(args.kegg)
    differential_genes = load_gene_list(args.gene_list)

    processed = process_pathways(
        pathways, reaction_to_k, k_to_gene, gene_annotations,
        differential_genes, args.max_gap, args.gap_count
    )

    write_output(args.output, processed)

if __name__ == "__main__":
    main()

