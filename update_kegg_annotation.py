#!/usr/bin/env python3

import sys

def update_kegg_annotation(renamed_genes_file, kegg_file, output_file):
    """
    Update gene IDs in KEGG annotation file based on renamed genes information
    
    Args:
        renamed_genes_file: Path to the file containing renamed gene information
        kegg_file: Path to the original KEGG annotation file
        output_file: Path to the output file with updated KEGG annotations
    """
    # Create mapping from original gene IDs to new gene IDs
    id_mapping = {}
    
    # Read the renamed genes file
    with open(renamed_genes_file, 'r') as f:
        # Skip header
        header = f.readline()
        
        # Read each line
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                new_id = fields[0]  # contig_1, contig_2, ...
                original_id = fields[1]  # MRS000001, MRS000002, ...
                
                # Build the mapping relationship
                id_mapping[original_id] = new_id
    
    # Read KEGG annotation file and update gene IDs
    updated_annotations = []
    not_found = []
    
    with open(kegg_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                original_id = fields[0]  # Original gene ID
                kegg_id = fields[1]      # KEGG ID
                
                # Find corresponding new gene ID
                if original_id in id_mapping:
                    new_id = id_mapping[original_id]
                    # Preserve all original columns, only replace the first column ID
                    updated_fields = [new_id] + fields[1:]
                    updated_annotations.append(updated_fields)
                else:
                    # Record genes without mapping
                    not_found.append(original_id)
                    # Keep original ID and all columns
                    updated_annotations.append(fields)
    
    # Write updated annotations to output file
    with open(output_file, 'w') as out:
        for fields in updated_annotations:
            out.write("\t".join(fields) + "\n")
    
    # Print processing results
    print(f"Annotation file update completed, processed {len(updated_annotations)} annotation records")
    if not_found:
        print(f"Warning: {len(not_found)} genes were not found in the renaming file")
        print("Genes without mapping: " + ", ".join(not_found[:5]) + 
              ("..." if len(not_found) > 5 else ""))
    print(f"Updated annotations saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) > 3:
        update_kegg_annotation(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        # Use default filenames
        update_kegg_annotation("renamed_genes.tsv", "kegg_annotation.txt", "updated_kegg_annotation.txt")
