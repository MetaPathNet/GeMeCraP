# GeMeCraP: Gene-Metabolite Cross-Referencing Assignment to Pathway

**GeMeCraP** (Gene-Metabolite Cross-Referencing Assignment to Pathway) is a novel computational framework designed to enable unambiguous, de novo elucidation of complete microbial metabolic pathways.

The gut microbiota produces hundreds of small-molecule metabolites that act as crucial signals for host-microbe interactions. While there is growing interest in harnessing these microbial metabolites for next-generation therapies, the biosynthetic pathways for many remain poorly characterized. Traditional approaches rely heavily on comparative genomics to identify genes and enzymes responsible for specific metabolite production. GeMeCraP offers an integrative, gene-metabolite network-based method to fill this gap.

---

## System Requirements

### Software Dependencies

- Python 3.12.3  
- `pandas` 2.2.3  
- `numpy` 2.2.1  
- `typing_extensions` 4.13.1  
- `jsonpatch` 1.33  
- `jsonpointer` 3.0.0  

### Hardware Recommendations

- CPU: 8+ cores  
- Memory: 64 GB RAM or more  
- Disk space: At least 250 GB  

Resource usage will vary depending on input dataset size.

---

## Usage Example

### Step 1: Build a Metabolite Network (Example: *Weissella paramesenteroides*)

1. **Cluster the metabolite `mz` values from the metabolomics dataset**

   ```bash
   python grouped_metabolites.py mz.txt
   ```

   This generates a clustered list of similar `mz` values. Extract the first column for downstream use:

   ```bash
   python grouped_metabolites.py mz.txt | awk '{print $1}' > mz_filter.txt
   ```

2. **Prepare the required input files**

   - `adduct_file.txt`: Ion adduct forms (e.g., +H, +Na, etc.)
   - `central.txt`: Central m/z values to define pathway anchors
   - `all_reaction.xls`: Reaction database

3. **Construct the metabolite network**

   ```bash
   python network_construct_adduct.py \
       --start_weight 118.0635457 \
       --end_weight 117.079 \
       central.txt mz_filter.txt all_reaction.xls adduct_file.txt > wp_network.txt
   ```

---

### Step 2: Map Metabolite Network to Gene Clusters

1. **Standardize gene names based on genome annotation**

   Gene clusters are identified based on gene proximity. By default, genes on the same contig (e.g., `contig_1`, `contig_2`) are considered neighbors, allowing for up to 3 intervening genes. If the gene annotations use alternative naming (e.g., `MRS000002`), you can reformat them as follows:

   ```bash
   python extract_protein_info.py wp.gff wp_gene_start_end.txt
   ```

2. **Update KEGG annotations with new gene names**

   ```bash
   python update_kegg_annotation.py wp_gene_start_end.txt wp.kegg wp_update.kegg
   ```

3. **Prepare the list of differentially expressed genes**

   - `gene_list.txt`: A list of genes identified as differentially expressed in transcriptome analysis.

4. **Identify potential gene clusters involved in the pathway**

   ```bash
   python find_cluster_terminal.py \
       --network wp_network.txt \
       --reaction all_reaction.xls \
       --kegg wp_update.kegg \
       --gene_list gene_list.txt \
       --output wp_cluster.txt
   ```

   This will generate `wp_cluster.txt` with predicted gene clusters associated with metabolite pathways.

---

## License

This project is open-source and available under the MIT License.

---

## Contact

For questions, suggestions, or contributions, please feel free to open an [issue](https://github.com/yourusername/GeMeCraP/issues) or contact the maintainer.
