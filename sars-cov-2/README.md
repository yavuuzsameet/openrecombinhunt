# SARS-CoV-2 Analysis Workflow

This directory contains the Jupyter Notebooks used to process the SARS-CoV-2 dataset and run RecombinHunt experiments.

## Pipeline Steps

The notebooks are designed to be run in the following order, assuming the initial `metadata.tsv` file is present in this directory.

### 1. Data Preparation

* **Notebook:** `prep.ipynb`
* **Input:** `metadata.tsv` (raw metadata)
* **Process:** Performs data quality checks, cleans the data based on QC metrics, and applies necessary transformations to the mutation columns.
* **Output:** `metadata_cleaned.tsv`

### 2. Environment and Sample Creation

* **Notebooks:** `env.ipynb` and `create_input_files.ipynb`
* **Input:** `metadata_cleaned.tsv`
* **Process:** Generates the required environment files (e.g., `lc_df.parquet`, `change2lineage_probability.parquet`) and the per-lineage sample files needed for the experiments.
* **Output:** Environment files and sample files.

### 3. Experiment Execution

* **Notebook:** `recombinant_cases_nextstrain.ipynb`
* **Input:** The environment and sample files created in the previous step.
* **Process:** Runs the RecombinHunt experiments on the prepared data.
* **Output:** Experiment results, analysis files, and reports.

## How to Use

To replicate the analysis, execute the notebooks sequentially, starting with `prep.ipynb`.