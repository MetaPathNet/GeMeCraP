import pandas as pd
import numpy as np
from collections import deque
import bisect
import argparse
import tqdm
import sys

def load_data(adduct_file, central_file, mz_file, diff_file):
    # Load data from the adduct_file
    # Adduct names and masses
    adduct_data = pd.read_csv(adduct_file, header=None, sep='\t', names=['adduct', 'mass'])

    # Load central file
    # Central compound masses
    central_data = pd.read_csv(central_file, header=None, names=['weight'])
    central_weights = central_data['weight'].tolist()

    # Initialize the ms_data dictionary
    # Initialize a dictionary to store central compound masses
    ms_data = {}
    # Add central.txt data to ms_data
    # Key: central compound mass, Value: central compound mass
    for weight in central_weights:
        key = f"{weight}"
        ms_data[key] = weight

    # Load mz_file
    mz_data = pd.read_csv(mz_file, header=None, names=['original_weight'])
    mz_weights = mz_data['original_weight'].tolist()
    # Generate corresponding adducts for mz_weights using adduct_data
    # Central is fixed, mz +/- adduct
    for original_weight in mz_weights:
        for _, row in adduct_data.iterrows():
            adduct_with_sign = row['adduct']  # "+H" or "-H"
            mass = row['mass']
            
            # Separate the sign and the actual ion name
            sign = adduct_with_sign[0]   # '+' or '-'
            ion = adduct_with_sign[1:]
            
            if sign == '+':
                # adjusted_weight = original_weight - mass
                adjusted_weight = original_weight - mass
                key = f"{original_weight}+{ion}"
            else:
                # If the sign is negative
                # adjusted_weight = original_weight + mass
                adjusted_weight = original_weight + mass
                key = f"{original_weight}-{ion}"
            
            ms_data[key] = adjusted_weight

    # Load diff_file
    # diff_values for target matching, diff_other_values for associated reactions
    diff_df = pd.read_csv(diff_file, sep='\t')
    diff_values = diff_df['diff_mass'].values
    diff_other_values = diff_df['ENTRY'].values

    return ms_data, (diff_values, diff_other_values)

'''Check if two nodes represent the same molecule (difference < 10ppm)'''
def is_same_molecule(weight1, weight2, threshold_ppm=10):
    return abs(weight1 - weight2) / weight1 * 1e6 < threshold_ppm

'''Get the actual molecule weight (removing the effect of adducts)'''
def get_actual_molecule_weight(key, ms_data):
    # For keys without adducts, return their value directly
    if '+' not in key and '-' not in key:
        return ms_data[key]
    # For keys with adducts, get their adjusted mass through ms_data
    return ms_data[key]

def find_matches(start_weight, end_weight, ms_data, diff_data, max_depth):
    diff_values, diff_other_values = diff_data

    '''Pre-calculate and cache all possible difference matches'''
    diff_mapping = {}
    diff_tolerance = 0.0065
    sorted_bounds = []
    for i, diff_val in enumerate(diff_values):
        lower_bound = diff_val - diff_tolerance
        upper_bound = diff_val + diff_tolerance
        key = (lower_bound, upper_bound)
        if key not in diff_mapping:
            diff_mapping[key] = []
            sorted_bounds.append(key)
        diff_mapping[key].append(diff_other_values[i])
    
    # Sort the boundaries for binary search
    sorted_bounds.sort(key=lambda x: x[0])  # Sort by lower bound
    lower_bounds = [bound[0] for bound in sorted_bounds]
    upper_bounds = [bound[1] for bound in sorted_bounds]
    
    queue = deque([(f"{start_weight}", 0, [], set(), set())])  # Initialize queue, add a set to store actual molecule weights
    results = []

    # Split the keys of ms_data and convert them to a numpy array of strings
    ms_keys = np.array(list(ms_data.keys()))
    ms_weights = np.array(list(ms_data.values()))

    candidate_base_weights = np.array([k.split('+')[0].split('-')[0] for k in ms_keys])
    
    # Create progress bar, set to dynamic mode
    pbar = tqdm.tqdm(dynamic_ncols=True, desc="Search Progress", leave=True)
    
    processed_nodes = 0
    update_interval = 100  # Update progress bar every 100 nodes
    max_queue_len = len(queue)
    found_paths = 0
    
    while queue:
        # Update maximum queue length
        max_queue_len = max(max_queue_len, len(queue))
        
        current_key, current_depth, path, used_nodes, used_molecules = queue.popleft()
        processed_nodes += 1
        
        # Periodically update progress bar
        if processed_nodes % update_interval == 0:
            pbar.set_postfix({
                'Processed': processed_nodes, 
                'Queue Size': len(queue), 
                'Max Queue': max_queue_len, 
                'Paths Found': found_paths
            })
            pbar.update(update_interval)

        if current_depth >= max_depth:
            continue

        # Extract the base weight of the current node
        if '+' in current_key:
            base_weight = current_key.split('+')[0]
        elif '-' in current_key:
            base_weight = current_key.split('-')[0]
        else:
            base_weight = current_key

        # Copy the set of used nodes to prevent sharing between paths
        used_nodes = used_nodes.copy()
        used_nodes.add(base_weight)

        current_weight_value = ms_data[current_key]
        '''Copy the set of used actual molecule weights'''
        used_molecules = used_molecules.copy()
        used_molecules.add(current_weight_value)

        # Exclude used base weights
        used_mask = np.isin(candidate_base_weights, list(used_nodes))
        valid_indices = np.where(~used_mask)[0]

        # Get valid candidate node information
        valid_weights = ms_weights[valid_indices]
        valid_keys = ms_keys[valid_indices]

        diffs = np.abs(valid_weights - current_weight_value)

        for i in range(len(valid_indices)):
            diff_ms = diffs[i]
            key = valid_keys[i]
            weight = valid_weights[i]

            '''Check if this new node's actual weight has already appeared in the path (within 10ppm tolerance)'''
            is_duplicate = False
            for used_mol_weight in used_molecules:
                if is_same_molecule(weight, used_mol_weight):
                    is_duplicate = True
                    break
            if is_duplicate:
                continue

            # Check if the difference is close to any value in diff mass
            diff_subset = []
            # Accelerate retrieval with binary search
            pos = bisect.bisect_right(lower_bounds, diff_ms)
            # Check ranges to the left
            for j in range(pos-1, -1, -1):
                if diff_ms >= lower_bounds[j] and diff_ms <= upper_bounds[j]:
                    diff_subset.extend(diff_mapping[(lower_bounds[j], upper_bounds[j])])
                # If current lower bound is too small, stop searching left
                if lower_bounds[j] < diff_ms - 2*diff_tolerance:
                    break
            # Check ranges to the right
            for j in range(pos, len(sorted_bounds)):
                if diff_ms >= lower_bounds[j] and diff_ms <= upper_bounds[j]:
                    diff_subset.extend(diff_mapping[(lower_bounds[j], upper_bounds[j])])
                # If current lower bound is too large, stop searching right
                if lower_bounds[j] > diff_ms + 2*diff_tolerance:
                    break
            if diff_subset:
                new_step = {'source': str(current_key), 'target': str(key), 'diff': diff_subset}
                new_path = path + [new_step]

                # check 20ppm
                if abs(weight - end_weight) / weight * 1e6 < 20:
                    results.append(new_path)
                    found_paths += 1
                else:
                    if current_depth + 1 < max_depth:
                        queue.append((key, current_depth + 1, new_path, used_nodes, used_molecules))
    
    # Update final progress
    pbar.set_postfix({
        'Processed': processed_nodes, 
        'Queue Size': len(queue), 
        'Max Queue': max_queue_len, 
        'Paths Found': found_paths
    })
    pbar.update(processed_nodes % update_interval)
    pbar.close()
    
    return results

def main():
    parser = argparse.ArgumentParser(description="Analyze metabolic networks from mass spectrometry data.")
    parser.add_argument("central_file", help="Path to the central data file.")
    parser.add_argument("mz_file", help="Path to the mz data file.")
    parser.add_argument("diff_file", help="Path to the reaction database file.")
    parser.add_argument("adduct_file", help="Path to the adduct file.")
    parser.add_argument("--start_weight", type=float, default=189.0789, help="Starting weight for matches.")
    parser.add_argument("--end_weight", type=float, default=203.0577, help="End weight to determine when to stop.")
    parser.add_argument("--max_depth", type=int, default=5, help="Maximum depth to search in the network, default depth = 5.")
    parser.add_argument("--output", type=str, default="output.txt", help="Output file path, default name output.txt.")
    
    args = parser.parse_args()

    print(f"Loading data files...")
    ms_data, diff_data = load_data(args.adduct_file, args.central_file, args.mz_file, args.diff_file)
    print(f"Data loading complete, starting path search...")
    
    paths = find_matches(args.start_weight, args.end_weight, ms_data, diff_data, args.max_depth)
    
    print(f"Found {len(paths)} paths, writing to {args.output}...")
    with open(args.output, 'w') as f:
        for path in paths:
            for step in path:
                f.write(str(step) + '\n')
            f.write("-" * 40 + '\n')  # Separator line between paths
    
    print(f"Successfully wrote {len(paths)} paths to {args.output}")

if __name__ == '__main__':
    main()
