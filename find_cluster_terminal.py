import argparse
import json
import pandas as pd
import ast
from typing import List, Dict, Set, Any, Tuple, Optional
from collections import defaultdict
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
    Each ENTRY may correspond to multiple K numbers, separated by commas.
    """
    reaction_to_k = {}
    try:
        reaction_df = pd.read_csv(reaction_file, sep="\t", dtype={'ENTRY': str, 'Orthology': str})
    except FileNotFoundError:
        logging.error(f"Reaction mapping file not found: {reaction_file}")
        exit(1)
    except Exception as e:
        logging.error(f"Error reading reaction mapping file: {e}")
        exit(1)

    for _, row in reaction_df.iterrows():
        entry = row["ENTRY"]
        orthology_field = row["Orthology"]
        if pd.notnull(entry) and pd.notnull(orthology_field):
            # Split K numbers by comma and strip any whitespace
            k_numbers = [k.strip() for k in orthology_field.split(",")]
            reaction_to_k[entry] = k_numbers
            logging.debug(f"Mapped ENTRY {entry} to K numbers: {k_numbers}")
    
    logging.info(f"Loaded reaction to K mapping for {len(reaction_to_k)} entries.")
    return reaction_to_k


def load_k_to_gene_mapping(kegg_file: str) -> Tuple[Dict[str, Set[str]], Set[str]]:
    """
    check k_number and genes
    """
    k_to_gene: Dict[str, Set[str]] = defaultdict(set)
    unknown_genes: Set[str] = set()
    try:
        with open(kegg_file, 'r', encoding='utf-8') as file:
            for line in file:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    gene = parts[0].strip()
                    k_number = parts[1].strip()
                    # check gene annotations
                    if k_number.startswith("K"):
                        k_to_gene[k_number].add(gene)
                    else:
                        # unknown genes set
                        unknown_genes.add(gene)
                        logging.debug(f"Gene {gene} is unknown due to invalid k_number: {k_number}")
    except FileNotFoundError:
        logging.error(f"KEGG mapping file not found: {kegg_file}")
        exit(1)
    except Exception as e:
        logging.error(f"Error reading KEGG mapping file: {e}")
        exit(1)
    
    logging.info(f"Loaded K to gene mapping for {len(k_to_gene)} K numbers.")
    logging.info(f"Identified {len(unknown_genes)} unknown genes.")
    return k_to_gene, unknown_genes


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
    total_edges: int,
    allow_unknown: bool,
    unknown_genes: Set[str],
    gap_count: int
) -> Optional[Dict[str, Any]]:
    """
    Validate the current cluster based on:
      1) Coverage and differential gene counts:
         - If 2 <= cluster size <= 3: must contain at least one differential gene, and cover all edges.
         - If cluster size > 3: must contain at least two differential genes, and cover all edges.
      2) Check for unknown genes if the following conditions are met:
         - allow_unknown == True
         - gap_count == 0
         - total_edges > 3
         - (total_edges - len(cluster_edges)) == 1
         In this case, look for unknown genes in the positions immediately to the left, right,
         and any missing positions in the middle.
      3) Return a dictionary {"genes": [...], "unknown_genes": [...]} if valid; otherwise return None.
    """

    # 1) Coverage + differential genes check
    if 2 <= len(cluster) <= 3:
        # Must have at least one differential gene and cover all edges
        valid = (differential_count >= 1) and (len(cluster_edges) == total_edges)
        if valid:
            return {"genes": cluster.copy(), "unknown_genes": []}

    elif len(cluster) > 3:
        # Must have at least two differential genes and cover all edges
        valid = (differential_count >= 2) and (len(cluster_edges) == total_edges)
        if valid:
            return {"genes": cluster.copy(), "unknown_genes": []}

    # 2) If not covered by the above coverage checks, try the unknown logic
    #    (some cases might not care about coverage if specifically looking for unknown patterns).
    if allow_unknown and gap_count == 0 and total_edges > 3 and (total_edges - len(cluster_edges) == 1):
        # Identify any unknown genes from left, right, or missing positions in the middle
        prefix = cluster[0].rsplit("_", 1)[0]  # prefix is consistent if they are in the same cluster
        suffixes = []
        for g in cluster:
            try:
                suffixes.append(int(g.rsplit("_", 1)[1]))
            except ValueError:
                # If suffix is not numeric, skip
                continue

        min_suf, max_suf = min(suffixes), max(suffixes)
        unknown_found = []

        # Left side
        left_gene = f"{prefix}_{min_suf - 1}"
        if left_gene in unknown_genes:
            unknown_found.append(left_gene)

        # Right side
        right_gene = f"{prefix}_{max_suf + 1}"
        if right_gene in unknown_genes:
            unknown_found.append(right_gene)

        # Missing positions in the middle
        for missing_suf in range(min_suf + 1, max_suf):
            if missing_suf not in suffixes:
                mid_gene = f"{prefix}_{missing_suf}"
                if mid_gene in unknown_genes:
                    unknown_found.append(mid_gene)

        # Regardless of whether unknown_found is empty or not, we consider this cluster "valid" under unknown logic.
        # If we only mark a cluster valid when unknown genes do exist, remove the comment block below:
        # if not unknown_found:
        #     return None
        if unknown_found:
            return {
                "genes": cluster.copy(),
                "unknown_genes": unknown_found  # could be an empty list if none found
            }

    # If none of the above conditions are met, the cluster is not valid
    return None

def identify_clusters(
    sorted_genes: List[str],
    differential_genes: Set[str],
    gene_to_edges: Dict[str, Set[int]],
    total_edges: int,
    max_gap: int,
    gap_count: int,
    allow_unknown: bool,
    unknown_genes: Set[str]
) -> List[Dict[str, Any]]:
    """
    Identify clusters of genes based on:
      - Gene adjacency (controlled by max_gap and gap_count).
      - The presence of differential genes.
      - Coverage of edges within the pathway.
      - Optional unknown gene detection if allow_unknown is True and gap_count=0.

    Returns a list of cluster dictionaries, each of the form:
      {
         "genes": [...],
         "unknown_genes": [...]
      }
    """
    clusters = []
    current_cluster = []
    differential_count = 0
    current_cluster_edges = set()
    current_gap_count = 0

    for gene in sorted_genes:
        if not current_cluster:
            # Start a new cluster
            current_cluster = [gene]
            differential_count = 1 if gene in differential_genes else 0
            current_cluster_edges = gene_to_edges.get(gene, set()).copy()
            current_gap_count = 0
            logging.debug(f"Started a new cluster with gene {gene}")
            continue

        prev_gene = current_cluster[-1]
        prev_parts = prev_gene.rsplit("_", 1)
        current_parts = gene.rsplit("_", 1)

        # If they have the same prefix, compute gap
        if prev_parts[0] == current_parts[0]:
            try:
                gap = int(current_parts[1]) - int(prev_parts[1])
            except ValueError:
                gap = 0  # If suffix is not numeric, treat as no gap

            if gap == 1:
                # Tightly connected (gap=1)
                current_cluster.append(gene)
                if gene in differential_genes:
                    differential_count += 1
                current_cluster_edges.update(gene_to_edges.get(gene, set()))
            elif 2 <= gap <= max_gap:
                # Allowed gap
                if gap_count == 1:
                    # Multiple gaps allowed
                    current_cluster.append(gene)
                    current_gap_count += 1
                    if gene in differential_genes:
                        differential_count += 1
                    current_cluster_edges.update(gene_to_edges.get(gene, set()))
                elif gap_count == 0:
                    # Only one gap allowed in total
                    if current_gap_count < 1:
                        current_cluster.append(gene)
                        current_gap_count += 1
                        if gene in differential_genes:
                            differential_count += 1
                        current_cluster_edges.update(gene_to_edges.get(gene, set()))
                    else:
                        # Cannot allow this gap; validate the current cluster and start a new one
                        cluster_check = validate_cluster(
                            current_cluster,
                            differential_count,
                            current_cluster_edges,
                            total_edges,
                            allow_unknown,
                            unknown_genes,
                            gap_count
                        )
                        if cluster_check is not None:
                            clusters.append(cluster_check)

                        # Start a new cluster
                        current_cluster = [gene]
                        differential_count = 1 if gene in differential_genes else 0
                        current_cluster_edges = gene_to_edges.get(gene, set()).copy()
                        current_gap_count = 0
            else:
                # gap is too large or negative, validate current cluster and start a new one
                cluster_check = validate_cluster(
                    current_cluster,
                    differential_count,
                    current_cluster_edges,
                    total_edges,
                    allow_unknown,
                    unknown_genes,
                    gap_count
                )
                if cluster_check is not None:
                    clusters.append(cluster_check)

                current_cluster = [gene]
                differential_count = 1 if gene in differential_genes else 0
                current_cluster_edges = gene_to_edges.get(gene, set()).copy()
                current_gap_count = 0

        else:
            # Different prefix, validate current cluster and start a new one
            cluster_check = validate_cluster(
                current_cluster,
                differential_count,
                current_cluster_edges,
                total_edges,
                allow_unknown,
                unknown_genes,
                gap_count
            )
            if cluster_check is not None:
                clusters.append(cluster_check)

            current_cluster = [gene]
            differential_count = 1 if gene in differential_genes else 0
            current_cluster_edges = gene_to_edges.get(gene, set()).copy()
            current_gap_count = 0

    # Check the last cluster after loop
    if current_cluster:
        cluster_check = validate_cluster(
            current_cluster,
            differential_count,
            current_cluster_edges,
            total_edges,
            allow_unknown,
            unknown_genes,
            gap_count
        )
        if cluster_check is not None:
            clusters.append(cluster_check)

    logging.info(f"Identified {len(clusters)} clusters.")
    return clusters

def process_pathways(
    pathways: List[List[Dict[str, Any]]],
    reaction_to_k: Dict[str, List[str]],
    k_to_gene: Dict[str, Set[str]],
    unknown_genes: Set[str],
    differential_genes: Set[str],
    max_gap: int,
    gap_count: int,
    allow_unknown: bool
) -> List[Dict[str, Any]]:
    """
    Process pathways to identify clusters and filter pathways based on the conditions.
    Each element of the returned list has:
      {
         "pathway_index": int,
         "pathway": [...original edges...],
         "clusters": [
            { "genes": [...], "unknown_genes": [...] },
            ...
         ]
      }
    """

    results = []
    for pathway_idx, pathway in enumerate(pathways):
        pathway_genes = set()
        gene_to_edges: Dict[str, Set[int]] = defaultdict(set)

        # Collect genes for the pathway
        for edge_idx, edge in enumerate(pathway):
            edge_genes = set()
            for reaction in edge.get("diff", []):
                k_numbers = reaction_to_k.get(reaction, [])
                for k in k_numbers:
                    genes = k_to_gene.get(k, set())
                    if genes:
                        edge_genes.update(genes)

            for gene in edge_genes:
                pathway_genes.add(gene)
                gene_to_edges[gene].add(edge_idx)

        total_edges = len(pathway)

        # Only process pathways if they contain at least one differential gene
        if pathway_genes & differential_genes:
            logging.info(f"Processing pathway {pathway_idx} with {len(pathway_genes)} genes and {total_edges} edges.")
            # Sort genes
            sorted_genes = sort_genes(pathway_genes)

            # Identify clusters
            clusters = identify_clusters(
                sorted_genes,
                differential_genes,
                gene_to_edges,
                total_edges,
                max_gap,
                gap_count,
                allow_unknown,
                unknown_genes
            )

            if clusters:
                results.append({
                    "pathway_index": pathway_idx,
                    "pathway": pathway,
                    "clusters": clusters
                })
                logging.info(f"Pathway {pathway_idx}: Identified {len(clusters)} clusters.")
        else:
            logging.info(f"Skipping pathway {pathway_idx} because it has no differential genes.")

    logging.info(f"Total pathways processed: {len(results)}")
    return results

def sort_genes(genes: Set[str]) -> List[str]:
    """
    Sort genes by their prefix and numeric suffix.
    Example: gene1_2 < gene1_3 < gene1_10 < gene2_1
    """
    def gene_key(gene: str):
        parts = gene.rsplit("_", 1)
        prefix = parts[0]
        try:
            suffix = int(parts[1])
        except (IndexError, ValueError):
            suffix = 0
        return prefix, suffix

    sorted_list = sorted(genes, key=gene_key)
    logging.debug(f"Sorted genes: {sorted_list}")
    return sorted_list

def write_output(output_file: str, processed_pathways: List[Dict[str, Any]]) -> None:
    separator = "----------------------------------------"
    try:
        with open(output_file, 'w', encoding='utf-8') as file:
            for pathway_info in processed_pathways:
                # 1) edges
                for edge in pathway_info["pathway"]:
                    source = edge.get("source", "")
                    target = edge.get("target", "")
                    diff = edge.get("diff", [])
                    edge_dict = {"source": source, "target": target, "diff": diff}
                    file.write(json.dumps(edge_dict, ensure_ascii=False) + "\n")
                
                # 2) clusters
                for cluster_obj in pathway_info["clusters"]:
                    cluster_str = ", ".join(cluster_obj["genes"])
                    
                    cluster_dict = {
                        "cluster_genes": cluster_str
                    }
                    
                    # if exist unknown_genes
                    if cluster_obj.get("unknown_genes"):
                        unknown_str = ", ".join(cluster_obj["unknown_genes"])
                        cluster_dict["unknown_genes"] = unknown_str

                    file.write(json.dumps(cluster_dict, ensure_ascii=False) + "\n")
                
                file.write(separator + "\n")
        logging.info(f"Output written to {output_file}")
    except Exception as e:
        logging.error(f"Error writing output file: {e}")
        exit(1)


def main():
    parser = argparse.ArgumentParser(description="Process metabolic network pathways and extract related gene clusters.")
    parser.add_argument("--network", required=True, help="Path to the metabolic network file.")
    parser.add_argument("--reaction", required=True, help="Path to the file containing all reactions.")
    parser.add_argument("--kegg", required=True, help="Path to the KEGG mapping file.")
    parser.add_argument("--gene_list", required=True, help="Path to the transcriptome differential gene list.")
    parser.add_argument("--output", required=True, help="Path to save the output results.")
    parser.add_argument("--max_gap", type=int, default=2, 
                        help="Maximum allowable gap between adjacent genes in a cluster when the cluster size is greater than 3. "
                             "For example, max_gap=2 allows up to one gene between adjacent genes.")
    parser.add_argument("--gap_count", type=int, choices=[0, 1], default=0,
                        help="Set to 1 to allow multiple gaps within a cluster or 0 to allow at most one gap.")
    parser.add_argument("--allow_unknown", action='store_true', 
                        help="Enable inclusion of unknown genes in clusters. This option is effective only when gap_count is set to 0.")
    
    args = parser.parse_args()


    # Load data
    pathways = parse_network(args.network)
    reaction_to_k = load_reaction_mapping(args.reaction)
    k_to_gene, unknown_genes = load_k_to_gene_mapping(args.kegg)
    differential_genes = load_gene_list(args.gene_list)

    # Process pathways
    processed_pathways = process_pathways(
        pathways,
        reaction_to_k,
        k_to_gene,
        unknown_genes,
        differential_genes,
        args.max_gap,
        args.gap_count,
        args.allow_unknown
    )

    # Write output
    write_output(args.output, processed_pathways)


if __name__ == "__main__":
    main()

