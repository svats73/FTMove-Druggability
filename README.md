# Druggability Analysis Tool

This Python program analyzes the druggability of binding sites using FTMap output and structural data. It provides functionalities to analyze single or multiple targets, generate percentage reports, and visualize results.

## Requirements

Ensure you have the following dependencies installed before running the program:

```sh
pip install numpy pandas seaborn tqdm matplotlib pymol-open-source
```

## Usage

The program is designed to be run from the command line using various subcommands.

### 1. Analyze a Single Target

```sh
python druggability_analysis.py analyze_single --scores path/to/scores.csv --target TARGET_NAME \
--bs_sesh_path path/to/binding_site_session.pse --ftmap_path path/to/ftmap/files \
[--bs BINDING_SITE] [--all_sites]
```

#### Arguments:
- `--scores`: Path to the CSV file containing binding site scores.
- `--target`: Target name for output file.
- `--bs_sesh_path`: Path to the PyMOL session file (.pse) for binding sites.
- `--ftmap_path`: Path to the FTMap files.
- `--bs`: (Optional) Specify a binding site number to analyze.
- `--all_sites`: (Optional) Analyze all binding sites if specified.

### 2. Analyze Multiple Targets

```sh
python druggability_analysis.py analyze_multi --data_folder path/to/data/folder
```

#### Arguments:
- `--data_folder`: Path to the folder containing multiple targets.

### 3. Generate Percentages from a Single Results File

```sh
python druggability_analysis.py generate_percentages --results_file path/to/results.csv
```

#### Arguments:
- `--results_file`: Path to the results CSV file.

### 4. Generate Percentages from Multiple Results Files

```sh
python druggability_analysis.py generate_percentages_multi --results_folder path/to/results/folder
```

#### Arguments:
- `--results_folder`: Path to the folder containing multiple results files.

### 5. Plot the Percentage of Binding Sites Meeting Criteria (Bar Chart)

```sh
python druggability_analysis.py plot_criteria_percentages_bar --percentage_file path/to/percentages.csv [--top_n N]
```

#### Arguments:
- `--percentage_file`: Path to the CSV file containing percentage data.
- `--top_n`: (Optional) Limit the number of binding sites plotted to top N.

### 6. Plot the Percentage of a Single Binding Site Meeting Criteria

```sh
python druggability_analysis.py plot_single_binding_site --percentage_file path/to/percentages.csv --binding_site BINDING_SITE
```

#### Arguments:
- `--percentage_file`: Path to the CSV file containing percentage data.
- `--binding_site`: Specific binding site to plot.

## Output

- The results of the druggability analysis will be saved as CSV files.
- Plots will be displayed using `matplotlib`.

## Example Commands

Analyze a single binding site:
```sh
python script.py analyze_single --scores data/scores.csv --target proteinX \
--bs_sesh_path data/proteinX_session.pse --ftmap_path data/ftmap --bs 3
```

Analyze all binding sites for a single target:
```sh
python script.py analyze_single --scores data/scores.csv --target proteinX \
--bs_sesh_path data/proteinX_session.pse --ftmap_path data/ftmap --all_sites
```

Analyze multiple targets in a dataset:
```sh
python script.py analyze_multi --data_folder data/all_targets
```

Generate percentages for a single results file:
```sh
python script.py generate_percentages --results_file results/proteinX_drug_analysis_output.csv
```

Generate percentages for multiple results files:
```sh
python script.py generate_percentages_multi --results_folder results/
```

Plot the percentage of binding sites meeting criteria:
```sh
python script.py plot_criteria_percentages_bar --percentage_file percentages/proteinX_percentage.csv
```

Plot the percentage of a single binding site:
```sh
python script.py plot_single_binding_site --percentage_file percentages/proteinX_percentage.csv --binding_site binding_site.3
```

---


