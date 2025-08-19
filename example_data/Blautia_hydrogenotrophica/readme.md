This is an example dataset for *Blautia hydrogenotrophica*.

Note that gene names in *B. hydrogenotrophica* follow a pattern typically separated by underscores and ending with numeric identifiers, so no renaming is required.

Due to the large size of the constructed metabolic network, it has not been uploaded.

Instead, we provide the result file `bh_cluster.txt`, which contains gene clusters validated through transcriptomic data.

In the metabolic gene cluster file, an entry like “273.0633+NH4” represents a positive ion mode, where the actual neutral mass should be calculated as “273.0633 - 18.034374”, with 18.034374 being the exact molecular weight of NH4.

## Usage Example

### Step 1: Build a Metabolite Network 

1. **Group the `mz` values from the metabolomics**

   ```bash
   python grouped_metabolites.py central.txt > central_group.txt
   python grouped_metabolites.py mz.txt > mz_group.txt
   ```

   This generates a grouped list of similar `mz` values. Extract the first column for downstream use:

   ```bash
   awk '{print $1}' central_group.txt >central_filter.txt
   awk '{print $1}' mz_group.txt > mz_filter.txt
   ```
2. **Construct the metabolite network**

   ```bash
   python network_construct_adduct.py \
       --start_weight 175.0637 \
       --end_weight 204.0898 \
       central_filter.txt mz_filter.txt all_reaction.xls adduct_file.txt --output bh_network.txt
   ```
---

### Step 2: Map Metabolic Network to Gene Clusters

1. **Identify potential gene clusters involved in the pathway**

   ```bash
   python find_cluster_terminal.py \
       --network bh_network.txt \
       --reaction all_reaction.xls \
       --kegg bh.kegg \
       --gene_list gene_list.txt \
       --output bh_cluster.txt
   ```
    Since the genome of *Blautia hydrogenotrophica* was assembled in-house and annotated using KEGG, the gene names fully conform to the format required for cluster identification (e.g., contig_2_1553, contig_2_1554, etc.), allowing this step to be performed directly.
   
   This will generate `bh_cluster.txt` with predicted gene clusters associated with metabolite pathways.

2. **checking shared MS/MS fragments**

   ```bash
   python cluster_MSMS.py \
       --cluster bh_cluster.txt \
       --central central_group.txt \
       --mz mz_group.txt \
       --msms MSMS.msp
   ```
    This step will print on screen the shared fragment information of adjacent metabolites.
---
