import os
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from pymol import cmd
import matplotlib.pyplot as plt
from pymol import stored
import copy
import tqdm
import re

def is_druggable(score, distances, max_dim):
    """
    Check if the primary site is druggable based on the criteria:
    1. The primary site must have at least 16 probe clusters.
    2. Distances between primary and secondary sites must be <= 8 Å.
    3. The maximum dimension of the binding site must be >= 10 Å.
    Additionally returns medium scoring (primary site has at least 13 probe clusters)
    and borderline (primary site has between 13 and 16 probe clusters)
    """
    #high_scoring = len(primary_site) >= 16
    high_scoring = score >= 16
    ccd = any(d[1] <= 8 for d in distances)
    maximum_distance = max_dim >= 10
    medium_scoring = score >= 13
    borderline = 16 > score >= 13
    
    return high_scoring, ccd, maximum_distance, medium_scoring, borderline

def calculate_center_of_mass(selection):
    """
    Calculate the center of mass (COM) for a given selection of atoms.
    """
    model = cmd.get_model(selection)
    total_mass = 0
    com = np.array([0.0, 0.0, 0.0])
    
    for atom in model.atom:
        mass = atom.get_mass()  # Get atomic mass
        total_mass += mass
        com += np.array(atom.coord) * mass
    
    if total_mass > 0:
        com /= total_mass  # Divide by total mass to get the COM
    return com

def calculate_centroid(selection):
    """
    Calculate the centroid for a given selection of atoms.
    """
    model = cmd.get_model(selection)
    centroid = np.array([0.0, 0.0, 0.0])
    num_atoms = 0

    for atom in model.atom:
        num_atoms += 1
        centroid += np.array(atom.coord)
    
    if num_atoms > 0:
        centroid /= num_atoms  # Divide by num atoms to get the centroid
    return centroid

def sites_within_cutoff_com(primary_site, secondary_sites, cutoff):
    """
    Measure the distances between the center of mass (COM) of the primary site
    and the secondary consensus sites. 
    """
    valid_sites = []
    
    # Calculate centroid for the primary site
    primary_com = calculate_center_of_mass(primary_site)
    
    for site in secondary_sites:
        if primary_site == site:
            continue
        else:
            secondary_com = calculate_center_of_mass(site)
            # Calculate Euclidean distance between the centroid
            dist = np.linalg.norm(primary_com - secondary_com)
            if dist <= cutoff:
                valid_sites.append(site)
    
    return valid_sites

def sites_within_cutoff_centroid(primary_site, secondary_sites, cutoff):
    """
    Measure the distances between the centroid of the primary site
    and the secondary consensus sites. 
    """
    valid_sites = []
    
    # Calculate centroid for the primary site
    primary_com = calculate_centroid(primary_site)
    
    for site in secondary_sites:
        if primary_site == site:
            continue
        else:
            secondary_com = calculate_centroid(site)
            # Calculate Euclidean distance between the centroid
            dist = np.linalg.norm(primary_com - secondary_com)
            if dist <= cutoff:
                valid_sites.append(site)
    
    return valid_sites

def measure_distances_com(primary_site, secondary_sites):
    """
    Measure the distances between the center of mass (COM) of the primary site
    and the secondary consensus sites. Heavy atoms, iteratively, more than one primary hotspot
    """
    distances = []
    
    # Calculate center of mass for the primary site
    primary_com = calculate_center_of_mass(primary_site)
    
    for site in secondary_sites:
        if primary_site == site:
            continue
        else:
            secondary_com = calculate_center_of_mass(site)
            # Calculate Euclidean distance between the centers of mass
            dist = np.linalg.norm(primary_com - secondary_com)
            distances.append((site, dist))
    
    return distances

def measure_distances_centroid(primary_site, secondary_sites):
    """
    Measure the distances between the centroid of the primary site
    and the secondary consensus sites. Heavy atoms, iteratively, more than one primary hotspot
    """
    distances = []
    
    # Calculate centroid for the primary site
    primary_com = calculate_centroid(primary_site)
    
    for site in secondary_sites:
        if primary_site == site:
            continue
        else:
            secondary_com = calculate_centroid(site)
            # Calculate Euclidean distance between the centroid
            dist = np.linalg.norm(primary_com - secondary_com)
            distances.append((site, dist))
    
    return distances


def calculate_max_dim(secondary_sites):
    """
    Calculate the maximum dimension (distance between the two farthest points) of the binding site
    based on the center of mass of all secondary consensus sites.
    """
    all_coords = []
    for site in secondary_sites:
        model = cmd.get_model(site)
        for atom in model.atom:
            all_coords.append(np.array(atom.coord))

    # If less than two atoms, we can't calculate a distance
    if len(all_coords) < 2:
        return 0.0

    # Calculate the maximum distance between all pairs of coordinates
    max_dim = 0.0
    for i, coord1 in enumerate(all_coords):
        for coord2 in all_coords[i + 1:]:
            dist = np.linalg.norm(coord1 - coord2)  # Euclidean distance
            if dist > max_dim:
                max_dim = dist

    return max_dim

def get_consensus_sites():
    """
    Extract the consensus sites from the FTMap output.
    This function will return a list of consensus site objects.
    """
    return [obj for obj in cmd.get_object_list() if obj.startswith('consensus.')]

def get_all_binding_sites():
    """
    Extract all available binding sites from the session.
    """
    return [obj for obj in cmd.get_object_list() if obj.startswith('binding_site.')]

def get_binding_sites(binding_site_number):
    """
    Extract the specific binding site based on user input.
    """
    return [obj for obj in cmd.get_object_list() if obj.startswith(f'binding_site.{binding_site_number}')]

def find_primary_consensus_site(consensus_sites, binding_site):
    """
    Find the primary consensus site for a given binding site.
    The primary site will be the one with the most probe clusters.
    """
    primary_consensus_site = None
    highest_rank = float('inf')
    for consensus_site in consensus_sites:
        #print(measure_distances(binding_site, [consensus_site]))
        if cmd.select(f'{consensus_site} within 2 of {binding_site}'):
            rank = int(consensus_site.split('.')[1])
            if rank < highest_rank:
                highest_rank = rank
                primary_consensus_site = consensus_site

    return primary_consensus_site

def write_output(dists, output_file):
    """
    Write the druggability analysis results to a CSV file.
    """
    df = pd.DataFrame(dists).T
    df.to_csv(output_file)
    print(f"Results written to {output_file}")


def calculate_all_criteria_percentages(df, criteria_list):
    """
    Calculate the percentage of binding sites meeting each specified criterion independently
    and store each percentage as a separate column in a new DataFrame.
    
    Parameters:
        df (DataFrame): The DataFrame with binding site data.
        criteria_list (list): List of criteria to analyze, each either a single string or a tuple of strings.
        
    Returns:
        DataFrame: A DataFrame with binding sites as rows and each criterion's percentage as columns.
    """
    # Create a new DataFrame for storing percentages
    percentage_df = pd.DataFrame(index=df['binding_site'].unique())
    #print(percentage_df)

    for criterion in criteria_list:
        # Create label based on whether the criterion is a tuple (multiple conditions) or a single string
        if isinstance(criterion, tuple):
            label = " & ".join(criterion)
        else:
            label = criterion
            criterion = (criterion,)

        # Calculate the percentage for the current criterion per binding site
        percentages = (
            df.groupby('binding_site')
              .apply(lambda x: (x[list(criterion)].all(axis=1).mean()) * 100)
        )
        
        # Add this percentage data to the DataFrame
        percentage_df[label] = percentages

    # Fill any NaN values that may have been created due to missing binding sites with 0
    percentage_df = percentage_df.fillna(0)
    return percentage_df

def plot_criteria_percentages_bar(args):
    """
    Plot the percentage of binding sites meeting each criterion using a bar chart.
    
    Parameters:
        percentage_df (DataFrame): DataFrame where each column is the percentage of binding sites 
                                   meeting that specific criterion.
        top_n (int, optional): If provided, only the first N binding sites will be plotted.
    """
    percentage_df = pd.read_csv(args.percentage_file)
    # If top_n is specified, select only the first N binding sites
    if args.top_n:
        percentage_df = percentage_df.iloc[:args.top_n]

    percentage_df.plot(kind='bar', figsize=(12, 8), width=0.8)

    plt.xlabel('Binding Site')
    plt.ylabel('Percentage of Binding Sites Meeting Criteria (%)')
    plt.title('Percentage of Binding Sites Meeting Each Criterion')
    plt.xticks(rotation=45, ha='right')
    plt.ylim(0, 100)
    plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
    plt.tight_layout()
    plt.show()


def plot_single_binding_site(args):
    """
    Plot the percentage of a single binding site meeting each criterion as a bar chart.
    
    Parameters:
        percentage_df (DataFrame): DataFrame where each column is the percentage of binding sites 
                                   meeting that specific criterion.
        binding_site (str): The specific binding site to plot.
    """
    percentage_df = pd.read_csv(args.percentage_file, index_col=0)
    binding_site = f"binding_site.{args.binding_site}"
    if binding_site not in percentage_df.index:
        print(f"Binding site '{binding_site}' not found.")
        return
    
    single_site_df = percentage_df.loc[binding_site]

    single_site_df.plot(kind='bar', figsize=(10, 6), width=0.6, color='skyblue')

    plt.xlabel('Criteria')
    plt.ylabel('Percentage of Binding Site Meeting Criteria (%)')
    plt.title(f'Percentage of Binding Site {binding_site} Meeting Each Criterion')
    plt.xticks(rotation=45, ha='right')
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.show()

def safe_parse(s):
    if isinstance(s, float) or isinstance(s, int):
        return ("unknown", float(s))
    
    # If it's a stringified float
    try:
        return ("unknown", float(s))
    except:
        pass
    
    # If it's a tuple with np.float64
    match = re.match(r"\('(.*?)',\s*np\.float64\((.*?)\)\)", s)
    if match:
        key = match.group(1)
        value = float(match.group(2))
        return (key, value)
    
    # Final fallback if it's a tuple with just a float (no np.float64)
    match = re.match(r"\('(.*?)',\s*(.*?)\)", s)
    if match:
        key = match.group(1)
        value = float(match.group(2))
        return (key, value)

    # Nothing matched
    raise ValueError(f"Unparsable string: {s}")

def plot_criteria_histogram_all(args):
    """
    Plot a histogram of the specified druggability criterion across all binding sites.
    Optionally, limit to the top_n binding sites.
    """
    df = pd.read_csv(args.results_file)
    index = str(args.results_file).find("_drug_analysis_output.csv")
    target = str(args.results_file)[0:index]

    if args.state == "bound":
        df = df[df['bound state'] == 'bound']
    elif args.state == "unbound":
        df = df[df['bound state'] == 'unbound']

    if args.criteria not in df.columns:
        print(f"Criterion '{args.criteria}' not found in the results file.")
        return

    bins = None
    if args.criteria == "high_scoring":
        column = "score"
        bins = np.arange(0, 41, 4)
    elif args.criteria == "ccd":
        column = "min_ccd"
        df['min_ccd'] = df['min_ccd'].apply(safe_parse)
        bins = np.arange(0, 31, 2)
    elif args.criteria == "maximum_distance":
        column = "max_dim"
        bins = np.arange(0, 31, 2)
    else:
        print(f"Invalid criterion '{args.criteria}' specified.")
        return

    if args.top_n:
        df = df.head(args.top_n)

    plt.figure(figsize=(10, 6))
    if args.criteria == "high_scoring" or args.criteria == "maximum_distance":
        plt.hist(df[column], bins=bins, edgecolor='black', rwidth=0.8)
    elif args.criteria == "ccd":
        values = df[column].apply(lambda x: x[1])
        plt.hist(values, bins=bins, edgecolor='black', rwidth=0.8)
    plt.xticks(bins)
    if not args.state:
        plt.title(f"Distribution of '{args.criteria}' for {target} All Binding Sites")
    else:
        plt.title(f"Distribution of '{args.criteria}' for {target} All Binding Sites ({args.state})")
    plt.xlabel(f"'{args.criteria}' Status")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.show()

def plot_criteria_histogram_single(args):
    """
    Plot a histogram of the specified druggability criterion for a single binding site.
    """
    df = pd.read_csv(args.results_file)
    index = str(args.results_file).find("_drug_analysis_output.csv")
    target = str(args.results_file)[0:index]

    if args.state == "bound":
        df = df[df['bound state'] == 'bound']
    elif args.state == "unbound":
        df = df[df['bound state'] == 'unbound']

    if args.criteria not in df.columns:
        print(f"Criterion '{args.criteria}' not found in the results file.")
        return

    bins = None
    if args.criteria == "high_scoring":
        column = "score"
        bins = np.arange(0, 41, 4)
    elif args.criteria == "ccd":
        column = "min_ccd"
        df['min_ccd'] = df['min_ccd'].apply(safe_parse)
        bins = np.arange(0, 31, 2)
    elif args.criteria == "maximum_distance":
        column = "max_dim"
        bins = np.arange(0, 31, 2)
    else:
        print(f"Invalid criterion '{args.criteria}' specified.")
        return


    binding_site = f"binding_site.{args.binding_site}"
    df_single = df[df['binding_site'] == binding_site]
    if df_single.empty:
        print(f"Binding site '{args.binding_site}' not found in the results file.")
        return

    plt.figure(figsize=(6, 4))
    if args.criteria == "high_scoring" or args.criteria == "maximum_distance":
        plt.hist(df_single[column], bins=bins, edgecolor='black', rwidth=0.8)
    elif args.criteria == "ccd":
        values = df_single[column].apply(lambda x: x[1])
        plt.hist(values, bins=bins, edgecolor='black', rwidth=0.8)
    plt.xticks(bins)
    if not args.state:
        plt.title(f"Distribution of '{args.criteria}' for {target} {args.binding_site}")
    else:
        plt.title(f"Distribution of '{args.criteria}' for {target} {args.binding_site} ({args.state})")
    plt.xlabel(f"'{args.criteria}' Status")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.show()

def bound_labeling_single(bound_states, pdbid, bs_sesh_path, binding_sites):
    cmd.reinitialize()
    if pdbid not in bound_states:
        bound_states[pdbid] = {}
    for binding_site in binding_sites:
        pdb = pdbid.split(".")[0] if '.' in pdbid else pdbid.split("_")[0]
        chain = pdbid.split(".")[1] if '.' in pdbid else pdbid.split("_")[1]
        #query = f"{pdb}_{chain}"

        #stored.results = []
        stored.res_ppi = []
        stored.res_mol = []

        bs_sesh_path = (args.bs_sesh_path)
        bs = (args.bs)

        cmd.load(bs_sesh_path)
        cmd.fetch(f"{pdb}")
        cmd.remove("solvent")
        cmd.align(f"{pdb} and chain {chain}", "ref_structure")
        cmd.select("ppi",f"{pdb} and not chain {chain} within 2 of binding_site.{bs}")
        cmd.select("lig",f"hetatm and not inorganic within 2 of binding_site.{bs}")
        cmd.iterate("ppi","stored.res_ppi.append(resn)")
        cmd.iterate("lig","stored.res_mol.append(resn)")

        # cmd.load(bs_sesh_path)
        # cmd.fetch(query)
        # cmd.remove("solvent")
        # cmd.align(query, "ref_structure")
        # cmd.select("lig", f"hetatm and not inorganic within 2 of binding_site.{binding_site.split('.')[1]}")
        # cmd.iterate("lig", "stored.results.append(resn)")

        kind = "bound" if stored.res_ppi or stored.res_mol else "unbound"
        if binding_site not in bound_states[pdbid]:
            bound_states[pdbid][binding_site] = kind
    


def analyze_single(args):

    # Read the scores and initialize output data structure
    df = pd.read_csv(args.scores)
    dists = {}

    # Load the PyMOL session and extract the specific binding site
    cmd.reinitialize()
    cmd.load(args.bs_sesh_path)
    
    if args.all_sites:
        binding_sites = get_all_binding_sites()  # Analyze all sites
    elif args.bs:
        binding_sites = get_binding_sites(args.bs)  # Analyze specific site
    else:
        raise ValueError("You must specify either --bs or --all_sites")

    if args.bound_states:
        bound_states = dict()
        for index, row in df.iterrows():
            pdb_id = row['Structure']
            bound_labeling_single(bound_states, pdb_id, args.bs_sesh_path, binding_sites)
            os.remove(f"{pdb_id.split('.')[0]}.cif") if '.' in pdb_id else os.remove(f"{pdb_id.split('_')[0]}.cif")

        cmd.reinitialize()
        cmd.load(args.bs_sesh_path)

    for index, row in df.iterrows():
        pdb_id = row['Structure']
        ftmap_path = args.ftmap_path
        print(f"Loading {ftmap_path}/{pdb_id}_aligned_ftmap.pdb")
        if Path(f'{ftmap_path}/{pdb_id}_aligned_ftmap.pdb').is_file():
            cmd.load(f"{ftmap_path}/{pdb_id}_aligned_ftmap.pdb")
        else:
            continue
        # Extract consensus sites
        consensus_sites = get_consensus_sites()

        for binding_site in binding_sites:
            site_num = 'Binding Site ' + str(binding_site.split('.')[1])
            bound_state = bound_states[pdb_id][binding_site] if args.bound_states else None
            score = row[site_num]
            key = str(pdb_id) + '_' + str(binding_site)
            if key not in dists:
                dists[key] = {}
            print(f"Analyzing binding site: {binding_site}")

            # Find the primary consensus site for the specified binding site
            primary_consensus_site = find_primary_consensus_site(consensus_sites, binding_site)
            print(f"Primary consensus site: {primary_consensus_site}")

            # Get secondary consensus sites within proximity
            new_binding_sites = cmd.get_object_list(f"consensus.* within 2 of {binding_site}")
            #distances = measure_distances_com(primary_consensus_site, new_binding_sites)
            distances = measure_distances_centroid(primary_consensus_site, new_binding_sites)
            #consensus_sites_within_binding_site = set(new_binding_sites)

            consensus_set = set(consensus_sites)
            consensus_set.discard(primary_consensus_site)
            #working_set = set(sites_within_cutoff_com(primary_consensus_site, list(consensus_set), 8))
            working_set = set(sites_within_cutoff_centroid(primary_consensus_site, list(consensus_set), 8))
            consensus_sites_within_binding_site = copy.deepcopy(working_set)
            consensus_sites_within_binding_site.add(primary_consensus_site)
            if primary_consensus_site in working_set:
                working_set.remove(primary_consensus_site)

            to_be_added = set()
            while len(working_set) != 0:
                for consensus_site in working_set:
                    if int(consensus_site.split('.')[2]) > 3:
                        working_consensus_set = set(consensus_sites)
                        working_consensus_set.discard(consensus_site)
                        #new_sites = set(sites_within_cutoff_com(consensus_site, list(working_consensus_set),8))
                        new_sites = set(sites_within_cutoff_centroid(consensus_site, list(working_consensus_set),8))
                        if len(new_sites) != 0:
                            consensus_sites_within_binding_site.update(new_sites)
                            to_be_added.update(new_sites)
                working_set = to_be_added
                to_be_added.clear()

            consensus_sites_within_binding_site = list(consensus_sites_within_binding_site)
                    
            if len(consensus_sites_within_binding_site) == 0:
                print(f"No consensus sites found near binding site {binding_site}")
                dists[key]['pdb_id'] = pdb_id
                dists[key]['bound state'] = bound_state
                dists[key]['binding_site'] = binding_site
                dists[key]['score'] = 0
                dists[key]['min_ccd'] = 0.0
                dists[key]['max_dim'] = 0.0
                dists[key]['druggability'] = False
                dists[key]['high_scoring'] = False
                dists[key]['ccd'] = False
                dists[key]['maximum_distance'] = False
                dists[key]['medium_scoring'] = False
                dists[key]['borderline'] = False
            else:
                print(f"Consensus sites within binding site {binding_site}: {consensus_sites_within_binding_site}")
                max_dim = calculate_max_dim(consensus_sites_within_binding_site)

                print(f"Measured distances: {distances}")
                print(f"Max dimension of binding site: {max_dim}")

                # Check if the site is druggable
                high_scoring, ccd, maximum_distance, medium_scoring, borderline = is_druggable(score, distances, max_dim)
                druggability = high_scoring and ccd and maximum_distance
                dists[key]['pdb_id'] = pdb_id
                dists[key]['bound state'] = bound_state
                dists[key]['binding_site'] = binding_site
                dists[key]['score'] = score 
                dists[key]['min_ccd'] = min(distances, key=lambda x: x[1]) if distances else 0.0
                dists[key]['max_dim'] = max_dim
                dists[key]['druggability'] = druggability
                dists[key]['high_scoring'] = high_scoring
                dists[key]['ccd'] = ccd
                dists[key]['maximum_distance'] = maximum_distance
                dists[key]['medium_scoring'] = medium_scoring
                dists[key]['borderline'] = borderline
                print(f"Druggability result for {pdb_id} in binding site {binding_site}: {druggability}")

        cmd.delete('consensus.*')

    # Write results to a CSV file
    output_file = f"{args.target}_drug_analysis_output.csv"
    write_output(dists, output_file)

def analyze_multi(args):
    subfolders = [f.name for f in os.scandir(args.data_folder) if f.is_dir()]

    subfolders_set = set(subfolders)
    subfolders_set.discard('results')
    subfolders = list(subfolders_set)

    for target in tqdm.tqdm(subfolders):
        scores = f"{args.data_folder}/{target}/{target}_scores.csv"
        bs_sesh_path = f"{args.data_folder}/{target}/binding_site_session.pse"
        target = target
        # Read the scores and initialize output data structure
        df = pd.read_csv(scores)

        dists = {}

        # Load the PyMOL session and extract the specific binding site
        cmd.reinitialize()
        cmd.load(bs_sesh_path)


        binding_sites = get_all_binding_sites()
        if args.bound_states:
            bound_states = dict()
            for index, row in df.iterrows():
                pdb_id = row['Structure']
                bound_labeling_single(bound_states, pdb_id, bs_sesh_path, binding_sites)
                os.remove(f"{pdb_id.split('.')[0]}.cif") if '.' in pdb_id else os.remove(f"{pdb_id.split('_')[0]}.cif")

            cmd.reinitialize()
            cmd.load(bs_sesh_path)

        for index, row in df.iterrows():
            pdb_id = row['Structure']
            ftmap_path = f"{args.data_folder}/{target}/ftmap_files"
            #print(f"Loading {ftmap_path}/{pdb_id}_aligned_ftmap.pdb")
            if Path(f'{ftmap_path}/{pdb_id}_aligned_ftmap.pdb').is_file():
                cmd.load(f"{ftmap_path}/{pdb_id}_aligned_ftmap.pdb")
            else:
                continue
            # Extract consensus sites
            consensus_sites = get_consensus_sites()

            for binding_site in binding_sites:  
                #bound_state = bound_labeling_single(pdb_id, bs_sesh_path, binding_site.split('.')[1])
                site_num = 'Binding Site ' + str(binding_site.split('.')[1])
                bound_state = bound_states[pdb_id][binding_site] if args.bound_states else None
                score = row[site_num]
                key = str(pdb_id) + '_' + str(binding_site)
                if key not in dists:
                    dists[key] = {}
                #print(f"Analyzing binding site: {binding_site}")

                # Find the primary consensus site for the specified binding site
                primary_consensus_site = find_primary_consensus_site(consensus_sites, binding_site)
                #print(f"Primary consensus site: {primary_consensus_site}")

                # Get secondary consensus sites within proximity
                new_binding_sites = cmd.get_object_list(f"consensus.* within 2 of {binding_site}")
                #distances = measure_distances_com(primary_consensus_site, new_binding_sites)
                distances = measure_distances_centroid(primary_consensus_site, new_binding_sites)
                #consensus_sites_within_binding_site = set(new_binding_sites)

                consensus_set = set(consensus_sites)
                consensus_set.discard(primary_consensus_site)
                #working_set = set(sites_within_cutoff_com(primary_consensus_site, list(consensus_set), 8))
                working_set = set(sites_within_cutoff_centroid(primary_consensus_site, list(consensus_set), 8))
                consensus_sites_within_binding_site = copy.deepcopy(working_set)
                consensus_sites_within_binding_site.add(primary_consensus_site)
                if primary_consensus_site in working_set:
                    working_set.remove(primary_consensus_site)

                to_be_added = set()
                while len(working_set) != 0:
                    for consensus_site in working_set:
                        if int(consensus_site.split('.')[2]) > 3:
                            working_consensus_set = set(consensus_sites)
                            working_consensus_set.discard(consensus_site)
                            #new_sites = set(sites_within_cutoff_com(consensus_site, list(working_consensus_set),8))
                            new_sites = set(sites_within_cutoff_centroid(consensus_site, list(working_consensus_set),8))
                            if len(new_sites) != 0:
                                consensus_sites_within_binding_site.update(new_sites)
                                to_be_added.update(new_sites)
                    working_set = to_be_added
                    to_be_added.clear()

                consensus_sites_within_binding_site = list(consensus_sites_within_binding_site)

                if len(consensus_sites_within_binding_site) == 0:
                    #print(f"No consensus sites found near binding site {binding_site}")
                    dists[key]['pdb_id'] = pdb_id
                    dists[key]['bound state'] = bound_state
                    dists[key]['binding_site'] = binding_site
                    dists[key]['score'] = 0
                    dists[key]['min_ccd'] = 0.0
                    dists[key]['max_dim'] = 0.0
                    dists[key]['druggability'] = False
                    dists[key]['high_scoring'] = False
                    dists[key]['ccd'] = False
                    dists[key]['maximum_distance'] = False
                    dists[key]['medium_scoring'] = False
                    dists[key]['borderline'] = False
                else:
                    #print(f"Consensus sites within binding site {binding_site}: {consensus_sites_within_binding_site}")
                    max_dim = calculate_max_dim(consensus_sites_within_binding_site)

                    #print(f"Measured distances: {distances}")
                    #print(f"Max dimension of binding site: {max_dim}")

                    # Check if the site is druggable
                    high_scoring, ccd, maximum_distance, medium_scoring, borderline = is_druggable(score, distances, max_dim)
                    druggability = high_scoring and ccd and maximum_distance
                    dists[key]['pdb_id'] = pdb_id
                    dists[key]['bound state'] = bound_state
                    dists[key]['binding_site'] = binding_site
                    dists[key]['score'] = score 
                    dists[key]['min_ccd'] = min(distances, key=lambda x: x[1]) if distances else 0.0
                    dists[key]['max_dim'] = max_dim
                    dists[key]['druggability'] = druggability
                    dists[key]['high_scoring'] = high_scoring
                    dists[key]['ccd'] = ccd
                    dists[key]['maximum_distance'] = maximum_distance
                    dists[key]['medium_scoring'] = medium_scoring
                    dists[key]['borderline'] = borderline
                    print(f"Druggability result for {pdb_id} in binding site {binding_site}: {druggability}")

            cmd.delete('consensus.*')

        # Write results to a CSV file
        if not os.path.exists(f"{args.data_folder}/results"):
            os.makedirs(f"{args.data_folder}/results")
        output_file = f"{args.data_folder}/results/{target}_drug_analysis_output.csv"
        write_output(dists, output_file)
        dists.clear()

def generate_percentages(args):
    # parser = argparse.ArgumentParser(description="Analyze druggability of binding sites based on FTMap output.")
    # parser.add_argument("--results_file", help="Path to the folder containing the data files.", required=True)
    # args = parser.parse_args()
    filename_str = args.results_file
    if Path(filename_str).is_file() and filename_str.endswith(".csv"):
        if args.state == "all":
            df = pd.read_csv(filename_str)
            df_bound = df[df['bound state'] == 'bound']
            df_unbound = df[df['bound state'] == 'unbound']
        elif args.state == "bound":
            df = pd.read_csv(filename_str)
            df = df[df['bound state'] == 'bound']
        elif args.state == "unbound":
            df = pd.read_csv(filename_str)
            df = df[df['bound state'] == 'unbound']
        else:
            df = pd.read_csv(filename_str)
        #df = pd.read_csv(filename_str)
        index = str(filename_str).find("_drug_analysis_output.csv")
        target = str(filename_str)[0:index]

        criteria_list = [
        ('high_scoring'),
        ('ccd'),
        ('maximum_distance'),
        ('medium_scoring'),
        ('borderline'),
        ('high_scoring', 'ccd'),
        ('high_scoring', 'maximum_distance'),
        ('ccd', 'maximum_distance'),
        ('high_scoring', 'ccd', 'maximum_distance')
        ]
        percentage_df = calculate_all_criteria_percentages(df, criteria_list) if not df.empty else pd.DataFrame()
        if not os.path.exists(f"percentages/"):
            os.makedirs(f"percentages/")
        if args.state == "all":
            percentage_df_bound = calculate_all_criteria_percentages(df_bound, criteria_list) if not df_bound.empty else pd.DataFrame()
            percentage_df_unbound = calculate_all_criteria_percentages(df_unbound, criteria_list) if not df_unbound.empty else pd.DataFrame()
            percentage_df_bound.to_csv(f"percentages/{target}_percentage_bound.csv")
            percentage_df_unbound.to_csv(f"percentages/{target}_percentage_unbound.csv")
            percentage_df.to_csv(f"percentages/{target}_percentage.csv")
        elif args.state == "bound":
            percentage_df.to_csv(f"{args.results_folder}/percentages/{target}_percentage_bound.csv")
        elif args.state == "unbound":
            percentage_df.to_csv(f"{args.results_folder}/percentages/{target}_percentage_unbound.csv")
        else:
            percentage_df.to_csv(f"percentages/{target}_percentage.csv")
        #percentage_df.to_csv(f"percentages/{target}_percentage.csv")

def generate_percentages_multi(args):
    # parser = argparse.ArgumentParser(description="Analyze druggability of binding sites based on FTMap output.")
    # parser.add_argument("--results_folder", help="Path to the folder containing the data files.", required=True)
    # args = parser.parse_args()
    for filename in os.listdir(args.results_folder):
        print(filename)
        filename_str = args.results_folder + '/' + filename
        if Path(filename_str).is_file() and filename_str.endswith(".csv"):
            if args.state == "all":
                df = pd.read_csv(filename_str)
                df_bound = df[df['bound state'] == 'bound']
                df_unbound = df[df['bound state'] == 'unbound']
            elif args.state == "bound":
                df = pd.read_csv(filename_str)
                df = df[df['bound state'] == 'bound']
            elif args.state == "unbound":
                df = pd.read_csv(filename_str)
                df = df[df['bound state'] == 'unbound']
            else:
                df = pd.read_csv(filename_str)
            index = str(filename).find("_drug_analysis_output.csv")
            target = str(filename)[0:index]

            criteria_list = [
            ('high_scoring'),
            ('ccd'),
            ('maximum_distance'),
            ('medium_scoring'),
            ('borderline'),
            ('high_scoring', 'ccd'),
            ('high_scoring', 'maximum_distance'),
            ('ccd', 'maximum_distance'),
            ('high_scoring', 'ccd', 'maximum_distance')
            ]
            percentage_df = calculate_all_criteria_percentages(df, criteria_list)  if not df.empty else pd.DataFrame()
            if not os.path.exists(f"{args.results_folder}/percentages"):
                os.makedirs(f"{args.results_folder}/percentages")
            if args.state == "all":
                percentage_df_bound = calculate_all_criteria_percentages(df_bound, criteria_list) if not df_bound.empty else pd.DataFrame()
                percentage_df_unbound = calculate_all_criteria_percentages(df_unbound, criteria_list) if not df_unbound.empty else pd.DataFrame()
                percentage_df_bound.to_csv(f"{args.results_folder}/percentages/{target}_percentage_bound.csv")
                percentage_df_unbound.to_csv(f"{args.results_folder}/percentages/{target}_percentage_unbound.csv")
                percentage_df.to_csv(f"{args.results_folder}/percentages/{target}_percentage.csv")
            elif args.state == "bound":
                percentage_df.to_csv(f"{args.results_folder}/percentages/{target}_percentage_bound.csv")
            elif args.state == "unbound":
                percentage_df.to_csv(f"{args.results_folder}/percentages/{target}_percentage_unbound.csv")
            else:
                percentage_df.to_csv(f"{args.results_folder}/percentages/{target}_percentage.csv")
            #percentage_df.to_csv(f"{args.results_folder}/percentages/{target}_percentage.csv")

def main():
    parser = argparse.ArgumentParser(description="Analyze druggability of binding sites based on FTMap output.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    
    parser_single = subparsers.add_parser("analyze_single", help="Analyze a single target")
    parser_single.add_argument("--scores", required=True, help="Path to the scores CSV file.")
    parser_single.add_argument("--target", required=True, help="Target name for output file.")
    parser_single.add_argument("--bs_sesh_path", required=True, help="Path to binding site session PSE file.")
    parser_single.add_argument("--ftmap_path", required=True, help="Path to FTMap files.")
    parser_single.add_argument("--bs", help="Binding site number to analyze.")
    parser_single.add_argument("--all_sites", action='store_true', help="Analyze all binding sites.")
    parser_single.add_argument("--bound_states", action='store_true', help="Calculate bound and unbound label for each structure")
    parser_single.set_defaults(func=analyze_single)
    
    parser_batch = subparsers.add_parser("analyze_multi", help="Analyze multiple targets in a directory")
    parser_batch.add_argument("--data_folder", required=True, help="Path to the data folder.")
    parser_batch.add_argument("--bound_states", action='store_true', help="Calculate bound and unbound label for each structure")
    parser_batch.set_defaults(func=analyze_multi)
    
    parser_percent = subparsers.add_parser("generate_percentages", help="Generate percentages from a results file")
    parser_percent.add_argument("--results_file", required=True, help="Path to the CSV results file.")
    parser_percent.add_argument("--state", choices=["bound", "unbound", "all"], 
                                help="Filter by 'bound' or 'unbound'. If not specified, use entire dataset.")
    parser_percent.set_defaults(func=generate_percentages)
    
    parser_percent_multi = subparsers.add_parser("generate_percentages_multi", help="Generate percentages from multiple results files")
    parser_percent_multi.add_argument("--results_folder", required=True, help="Path to the folder containing results files.")
    parser_percent_multi.add_argument("--state", choices=["bound", "unbound", "all"], 
                                help="Filter by 'bound' or 'unbound'. If not specified, use entire dataset.")
    parser_percent_multi.set_defaults(func=generate_percentages_multi)

    parser_plot = subparsers.add_parser("plot_criteria_percentages_bar", help="Plot the percentage of binding sites meeting each criterion")
    parser_plot.add_argument("--percentage_file", required=True, help="Path to the percentages CSV file.")
    parser_plot.add_argument("--top_n", type=int, help="Number of binding sites to plot.")
    parser_plot.set_defaults(func=plot_criteria_percentages_bar)

    parser_plot_single = subparsers.add_parser("plot_single_binding_site", help="Plot the percentage of a single binding site meeting each criterion")
    parser_plot_single.add_argument("--percentage_file", required=True, help="Path to the percentages CSV file.")
    parser_plot_single.add_argument("--binding_site", required=True, help="Binding site to plot.")
    parser_plot_single.set_defaults(func=plot_single_binding_site)

    parser_hist_all = subparsers.add_parser("plot_criteria_histogram_all", help="Plot histogram of a criterion across all binding sites")
    parser_hist_all.add_argument("--results_file", required=True, help="Path to the drug analysis results CSV file.")
    parser_hist_all.add_argument("--criteria", required=True, help="The druggability criterion to plot.")
    parser_hist_all.add_argument("--top_n", type=int, help="Limit to the top N binding sites.")
    parser_hist_all.add_argument("--state", choices=["bound", "unbound"],
                                 help="Filter by 'bound' or 'unbound'. If not specified, use entire dataset.")
    parser_hist_all.set_defaults(func=plot_criteria_histogram_all)

    parser_hist_single = subparsers.add_parser("plot_criteria_histogram_single", help="Plot histogram of a criterion for a single binding site")
    parser_hist_single.add_argument("--results_file", required=True, help="Path to the drug analysis results CSV file.")
    parser_hist_single.add_argument("--criteria", required=True, help="The druggability criterion to plot.")
    parser_hist_single.add_argument("--binding_site", required=True, help="Specific binding site to plot (e.g., binding_site.3).")
    parser_hist_single.add_argument("--state", choices=["bound", "unbound"],
                                 help="Filter by 'bound' or 'unbound'. If not specified, use entire dataset.")
    parser_hist_single.set_defaults(func=plot_criteria_histogram_single)
    
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
