This is an example dataset for *Clostridium sporogenes*.

This directory includes:
- Gene annotation information  
- Gene name revised: `gene_rename.txt`
- Upregulated genes from transcriptomic data: `gene_list.txt`  
- Metabolic gene clusters: `cs_cluster.txt`  

The metabolic network file has not been uploaded due to its large size.

In the metabolic gene cluster file, entries like `188.0707+H` represent positive ion mode. The actual neutral mass should be calculated as `188.0707 - 1.007825`, where 1.007825 is the exact molecular weight of H.

"cluster_genes": "contig2_280, contig2_281, contig2_282, contig2_283, contig2_284, contig2_285, contig2_286, contig2_287, contig2_288, contig2_289" is target cluster.

## Usage Example

### Step 1: Build a Metabolite Network 

1. **Group the `mz` values from the metabolomics**

   ```bash
   python grouped_metabolites.py mz.txt
   ```
  
   ```bash
   python grouped_metabolites.py mz.txt | awk '{print $1}' > mz_filter.txt
   ```
   Some metabolites were confirmed to be in the positive ion mode. After calculating their accurate m/z values, they were incorporated into the central.txt file for subsequent clustering analysis.
   ```bash
   python grouped_metabolites.py central.txt | awk '{print $1}' > central_filter.txt
   ```
3. **Construct the metabolite network**

   ```bash
   python network_construct_adduct.py \
       --start_weight 189.0789 \
       --end_weight 203.0577 \
       central_filter.txt mz_filter.txt all_reaction.xls adduct_file.txt --output cs_network.txt
   ```
---

### Step 2: Map Metabolic Network to Gene Clusters

1. **Identify potential gene clusters involved in the pathway**

   ```bash
   python find_cluster_terminal.py \
       --network cs_network.txt \
       --reaction all_reaction.xls \
       --kegg cs.kegg \
       --gene_list gene_list.txt \
       --max_gap 4 \
       --output cs_cluster.txt
   ```
    Since the revised gene names of *Clostridium sporogenes* were fully conformed to the format required for cluster identification, allowing this step to be performed directly.
   
   This will generate `cs_cluster.txt` with predicted gene clusters associated with metabolite pathways.

---
