import pandas as pd
import argparse
from collections import deque
import numpy as np

def load_data(central_file, mz_file, diff_file, adduct_file):
    # Load data from the adduct_file
    adduct_data = pd.read_csv(adduct_file, header=None, sep='\t', names=['adduct', 'mass'])
    #adduct_data['mass'] = adduct_data['mass'].astype(float)
    
    # load central
    central_data = pd.read_csv(central_file, header=None, names=['weight'])
    central_weights = central_data['weight'].tolist()
    
    # Initialize the ms_data dictionary
    ms_data = {}
    
    # Add central.txt data to ms_data
    for weight in central_weights:
        key = f"{weight}"
        ms_data[key] = weight
    
    mz_data = pd.read_csv(mz_file, header=None, names=['original_weight'])
    mz_weights = mz_data['original_weight'].tolist()
    
    # Generate corresponding adducts for mz_weights using adduct_data
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

    # load diff_file
    diff_df = pd.read_csv(diff_file, sep='\t')
    diff_values = diff_df['diff_mass'].values
    diff_other_values = diff_df['ENTRY'].values
    
    return ms_data, (diff_values, diff_other_values)


def find_matches(start_weight, end_weight, ms_data, diff_data, max_depth):
    diff_values, diff_other_values = diff_data

    queue = deque([(f"{start_weight}", 0, [], set())])  # Initialize the queue
    results = []

    # Split the keys of ms_data and convert them to a numpy array of strings
    ms_keys = np.array(list(ms_data.keys()))
    ms_weights = np.array(list(ms_data.values()))

    candidate_base_weights = np.array([k.split('+')[0].split('-')[0] for k in ms_keys])

    while queue:
        current_key, current_depth, path, used_nodes = queue.popleft()
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

            # Check if the difference is close to any value in diff mass
            close_mask = np.abs(diff_values - diff_ms) <= 0.005
            if np.any(close_mask):
                diff_subset = diff_other_values[close_mask].tolist()
                new_step = {'source': str(current_key), 'target': str(key), 'diff': diff_subset}
                new_path = path + [new_step]

                # check 20ppm
                if abs(weight - end_weight) / weight * 1e6 < 20:
                    results.append(new_path)
                else:
                    if current_depth + 1 < max_depth:
                        queue.append((key, current_depth + 1, new_path, used_nodes))
    return results


def main():
    parser = argparse.ArgumentParser(description="Analyze metabolic networks from mass spectrometry data.")
    parser.add_argument("central_file", help="Path to the central data file.")
    parser.add_argument("mz_file", help="Path to the mz data file.")
    parser.add_argument("diff_file", help="Path to the differential formula file.")
    parser.add_argument("adduct_file", help="Path to the adduct file.")
    parser.add_argument("--start_weight", type=float, default=175.0634, help="Starting weight for matches.")
    parser.add_argument("--end_weight", type=float, default=204.0905, help="End weight to determine when to stop.")
    parser.add_argument("--max_depth", type=int, default=5, help="Maximum depth to search in the network.")
    
    args = parser.parse_args()

    ms_data, diff_data = load_data(args.central_file, args.mz_file, args.diff_file, args.adduct_file)
    paths = find_matches(args.start_weight, args.end_weight, ms_data, diff_data, args.max_depth)
    
    for path in paths:
        for step in path:
            print(step)
        print("-" * 40)  # Separator line between paths

if __name__ == '__main__':
    main()

