# Druggability Analysis Tool

This Python program analyzes the druggability of binding sites using FTMap output and structural data. It provides functions to analyze single or multiple targets, generate percentage reports, and visualize results.

## Requirements

Ensure you have the following dependencies installed before running the program:

```sh
pip install numpy pandas tqdm matplotlib pymol-open-source
```

## Data Availability

All data used for analysis can be found in the data folder. pdb_ids.xlsx contains the name of the target and its corresponding reference PDB ID used for analysis.
The targets results used for analysis are broken down by the source data set with each target folder containing its respective analysis information. For each target, its folder contains the FTMap results/outputs, FTMove binding site meshes in pse format, and hotspot scores for each target binding site in csv format.

## Usage

The program is designed to be run from the command line using various subcommands.

---

### 1. Analyze a Single Target

```sh
python druggability_analysis.py analyze_single --scores path/to/scores.csv --target TARGET_NAME \
--bs_sesh_path path/to/binding_site_session.pse --ftmap_path path/to/ftmap/files \
[--bs BINDING_SITE] [--all_sites] [--bound_states]
```

**Arguments:**
- `--scores`: Path to the CSV file containing binding site scores.
- `--target`: Target name for output file.
- `--bs_sesh_path`: Path to the PyMOL session file (.pse) for binding sites.
- `--ftmap_path`: Path to the FTMap files.
- `--bs`: (Optional) Specify a single binding site number to analyze.
- `--all_sites`: (Optional) Analyze all binding sites.
- `--bound_states`: (Optional) Label structures as 'bound' or 'unbound' based on the presence of ligand or interface residues.

---

### 2. Analyze Multiple Targets

```sh
python druggability_analysis.py analyze_multi --data_folder path/to/data/folder [--bound_states]
```

**Arguments:**
- `--data_folder`: Path to the main folder containing multiple targets (see structure below).
- `--bound_states`: (Optional) Label structures as 'bound' or 'unbound'.

**Expected Folder Structure for `analyze_multi`:**

```
data_folder/
├── Target1/
│   ├── Target1_scores.csv
│   ├── binding_site_session.pse
│   └── ftmap_files/
│       ├── pdbid1_aligned_ftmap.pdb
│       ├── pdbid2_aligned_ftmap.pdb
│       └── ...
├── Target2/
│   ├── Target2_scores.csv
│   ├── binding_site_session.pse
│   └── ftmap_files/
│       └── ...
└── results/
    └── (created automatically to store output CSVs)
```

---

### 3. Generate Percentages from a Single Results File

```sh
python druggability_analysis.py generate_percentages --results_file path/to/results.csv [--state {bound,unbound,all}]
```

**Arguments:**
- `--results_file`: Path to the druggability results CSV file.
- `--state`: (Optional) Filter percentages for only 'bound', 'unbound', or 'all' structures. 'all' will generate percentages for bound, unbound, and all structures.

---

### 4. Generate Percentages from Multiple Results Files

```sh
python druggability_analysis.py generate_percentages_multi --results_folder path/to/results/folder [--state {bound,unbound,all}]
```

**Arguments:**
- `--results_folder`: Path to the folder containing multiple druggability result files.
- `--state`: (Optional) Filter percentages for only 'bound', 'unbound', or 'all' structures. 'all' will generate percentages for bound, unbound, and all structures.

---

### 5. Plot the Percentage of Binding Sites Meeting Criteria (Bar Chart)

```sh
python druggability_analysis.py plot_criteria_percentages_bar --percentage_file path/to/percentages.csv [--top_n N]
```

**Arguments:**
- `--percentage_file`: Path to the percentages CSV file.
- `--top_n`: (Optional) Plot the top n binding site percentages 

---

### 6. Plot the Percentage of a Single Binding Site Meeting Criteria

```sh
python druggability_analysis.py plot_single_binding_site --percentage_file path/to/percentages.csv --binding_site BINDING_SITE_NUMBER (3 digit)
```

**Arguments:**
- `--percentage_file`: Path to the percentages CSV file.
- `--binding_site`: Number of the binding site to plot

---

### 7. Plot Histogram of a Druggability Criterion Across All Binding Sites

```sh
python druggability_analysis.py plot_criteria_histogram_all --results_file path/to/results.csv --criteria CRITERION [--state {bound,unbound}]
```

**Arguments:**
- `--results_file`: Path to the druggability analysis results file to plot.
- `--criteria`: Specify which druggability criteria to plot, options are `high_scoring`, `ccd`, and `maximum distance`
- `--state`: (Optional) Filter results to only 'bound' or 'unbound' for plotting

---

### 8. Plot Histogram of a Druggability Criterion for a Single Binding Site

```sh
python druggability_analysis.py plot_criteria_histogram_single --results_file path/to/results.csv --criteria CRITERION --binding_site BINDING_SITE [--state {bound,unbound}]
```

**Arguments:**
- `--results_file`: Path to the druggability analysis results file to plot.
- `--criteria`: Specify which druggability criteria to plot, options are `high_scoring`, `ccd`, and `maximum distance`
- `--binding_site`: Number of the binding site to plot
- `--state`: (Optional) Filter results to only 'bound' or 'unbound' for plotting

---

## Output

- Druggability results will be saved as `.csv` files (e.g., `TARGET_drug_analysis_output.csv`).
- Percentages will be saved in a `percentages/` folder, the `percentages/` folder will be in the results folder if working with multiple targets.
- Plots will be displayed using `matplotlib`.

---

## Example Commands

Analyze a single binding site:

```sh
python druggability_analysis.py analyze_single --scores data/proteinX/proteinX_scores.csv --target proteinX \
--bs_sesh_path data/proteinX/binding_site_session.pse --ftmap_path data/proteinX/ftmap_files --bs 3
```

Analyze all binding sites:

```sh
python druggability_analysis.py analyze_single --scores data/proteinX/proteinX_scores.csv --target proteinX \
--bs_sesh_path data/proteinX/binding_site_session.pse --ftmap_path data/proteinX/ftmap_files --all_sites
```

Batch analyze all targets:

```sh
python druggability_analysis.py analyze_multi --data_folder data/all_targets --bound_states
```

Generate percentages for a single results file:

```sh
python druggability_analysis.py generate_percentages --results_file results/proteinX_drug_analysis_output.csv --state all
```

Generate for multiple:

```sh
python druggability_analysis.py generate_percentages_multi --results_folder results/
```

Bar chart plot:

```sh
python druggability_analysis.py plot_criteria_percentages_bar --percentage_file percentages/proteinX_percentage.csv
```

Single binding site plot:

```sh
python druggability_analysis.py plot_single_binding_site --percentage_file percentages/proteinX_percentage.csv --binding_site 003
```

Histogram across all:

```sh
python druggability_analysis.py plot_criteria_histogram_all --results_file results/proteinX_drug_analysis_output.csv --criteria high_scoring
```

Histogram for one:

```sh
python druggability_analysis.py plot_criteria_histogram_single --results_file results/proteinX_drug_analysis_output.csv --criteria ccd --binding_site 003
```
