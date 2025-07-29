# app/streamlit_app.py

import streamlit as st
import os
from pathlib import Path
import pandas as pd
import yaml
from streamlit_option_menu import option_menu
import re
import io
import json
from bs4 import BeautifulSoup
import sys
import plotly.graph_objects as go

# --- Path Setup for Imports ---
# Add the project root to the Python path to allow imports from 'src'
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.append(str(PROJECT_ROOT))

try:
    from src.utils.constants import *
except ImportError as e:
    print(f"CRITICAL ERROR: Could not import from 'src/utils/constants.py'. {e}")
    sys.exit(1)

# --- Configuration & Path Setup ---
CONFIG_PATH = Path("config/config.yaml")

try:
    with open(CONFIG_PATH, 'r') as f:
        config = yaml.safe_load(f)
    RESULTS_DIR_BASE = Path(config.get(PATHS).get(RESULTS))
    PROCESSED_DATA_DIR_BASE = Path(config.get(PATHS).get(PROCESSED_DATA))
except FileNotFoundError:
    st.error(f"Configuration file not found at {CONFIG_PATH}. Please ensure the file exists.")
    st.stop()
except Exception as e:
    st.error(f"Error loading or parsing configuration file: {e}")
    st.stop()

# --- Helper Functions ---

@st.cache_data
def get_available_viruses(results_dir: Path) -> list:
    """Scans the results directory to find which viruses have processed results."""
    virus_list = []
    recombinhunt_path = results_dir / "recombinhunt_output"
    if recombinhunt_path.is_dir():
        for entry in os.scandir(recombinhunt_path):
            if entry.is_dir():
                virus_list.append(entry.name)
    return sorted(virus_list)

def load_css(file_path: Path):
    """Loads a CSS file and injects it into the Streamlit app."""
    try:
        with open(file_path) as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    except FileNotFoundError:
        st.error(f"CSS file not found at {file_path}. Please ensure it exists.")

@st.cache_data
def load_analysis_data(virus_name: str, param_set: str) -> tuple:
    """Loads all necessary data files for a given virus and parameter set."""
    report_dir = RESULTS_DIR_BASE / "recombinhunt_output" / virus_name / param_set
    analysis_dir = RESULTS_DIR_BASE / "analysis_results" / virus_name / param_set
    metadata_path = PROCESSED_DATA_DIR_BASE / virus_name / "metadata.tsv"

    # Load extensive summary table for the main dashboard
    extensive_summary_path = report_dir / "recombinant_summary_extensive.tsv"
    extensive_df = pd.DataFrame()
    if extensive_summary_path.exists():
        extensive_df = pd.read_csv(extensive_summary_path, sep='\t')

    # Load filtering summary and merge with metadata for the explorer
    filtering_summary_path = report_dir / "recombinant_summary.tsv"
    explorer_df = pd.DataFrame()
    if filtering_summary_path.exists():
        filtering_df = pd.read_csv(filtering_summary_path, sep='\t')
        if metadata_path.exists():
            metadata_df = pd.read_csv(metadata_path, sep='\t', usecols=['Virus name', 'Collection date', 'Location'])
            # Use the first genomeID for merging
            filtering_df['merge_id'] = filtering_df['genomeIDs'].str.split(',').str[0]
            explorer_df = pd.merge(filtering_df, metadata_df, left_on='merge_id', right_on='Virus name', how='left').drop(columns=['merge_id', 'Virus name'])
        else:
            explorer_df = filtering_df # Proceed without metadata
    
    # Paths to plot images
    plots = {
        "recomb_plot": analysis_dir / "distribution_by_lineage_recombinants_only.png",
        "date_plot": analysis_dir / "distribution_by_date.png",
        "location_plot": analysis_dir / "distribution_by_country_heatmap.png"
    }
    
    return extensive_df, explorer_df, plots

@st.cache_data
def load_case_details(case_folder_path_str: str) -> dict:
    """Loads all structured report files for a single case."""
    case_folder_path = Path(case_folder_path_str)
    details = {}
    
    # Load summary JSON
    summary_path = case_folder_path / "summary.json"
    if summary_path.exists():
        with open(summary_path, 'r') as f:
            details['summary'] = json.load(f)
            
    # Load region tables
    details['region_tables'] = []
    for f in sorted(case_folder_path.glob("region_*_table.csv")):
        details['region_tables'].append(pd.read_csv(f))
        
    # Load Plotly JSON
    plot_per_region_path = case_folder_path / "plot_per_region.json"
    if plot_per_region_path.exists():
        with open(plot_per_region_path, 'r') as f:
            details['plot_per_region'] = go.Figure(json.load(f))

    plot_whole_genome_path = case_folder_path / "plot_whole_genome.json"
    if plot_whole_genome_path.exists():
        with open(plot_whole_genome_path, 'r') as f:
            details['plot_whole_genome'] = go.Figure(json.load(f))
            
    return details

# --- Main Application ---
def main():
    st.set_page_config(page_title="OpenRecombinHunt Results", page_icon="🧬", layout="wide")
    load_css(Path(__file__).parent / "style.css")
    st.title("OpenRecombinHunt Pipeline Results")

    with st.sidebar:
        st.header("Virus Selection")
        available_viruses = get_available_viruses(RESULTS_DIR_BASE)
        if not available_viruses:
            st.warning("No virus results found."); st.stop()

        menu_options = ["— Welcome —"] + [v.replace('-', ' ').title() for v in available_viruses]
        menu_icons = ["house-door-fill"] + ["virus"] * len(available_viruses)
        selected_option = option_menu(
            menu_title=None, options=menu_options, icons=menu_icons, menu_icon="cast", default_index=0,
            styles={
                "container": {"padding": "0!important", "background-color": "#fafafa"},
                "icon": {"color": "#4A90E2", "font-size": "18px"}, 
                "nav-link": {"font-size": "16px", "text-align": "left", "margin":"0px", "--hover-color": "#eee"},
                "nav-link-selected": {"background-color": "#4A90E2"},
            }
        )

    if selected_option == "— Welcome —":
        st.header("Welcome to the OpenRecombinHunt Analysis Dashboard")
        st.markdown("This application presents results from the automated recombination analysis pipeline.")
        st.info("Please select a virus from the navigation menu on the left to view its detailed analysis.")
    else:
        selected_virus = selected_option.replace(' ', '-').lower()
        st.header(f"Analysis for: {selected_option}")
        
        virus_results_path = RESULTS_DIR_BASE / "recombinhunt_output" / selected_virus
        if virus_results_path.is_dir():
            param_sets = [d.name for d in virus_results_path.iterdir() if d.is_dir()]
            if param_sets:
                latest_param_set = param_sets[0]
                
                extensive_df, explorer_df, plots = load_analysis_data(selected_virus, latest_param_set)
                
                if extensive_df.empty and explorer_df.empty:
                    st.warning("No recombinant summary data was found for this virus run.")
                    st.stop()

                tab1, tab2 = st.tabs(["📊 Summary Dashboard", "🔬 Recombinant Explorer"])

                with tab1:
                    st.subheader("Key Metrics")
                    total_recombinants = len(extensive_df)
                    total_sequences = extensive_df['n_sequences'].sum() if 'n_sequences' in extensive_df.columns else "N/A"
                    
                    col1, col2, col3 = st.columns(3)
                    col1.metric("Total Sequences Analyzed", f"{total_sequences:,}" if isinstance(total_sequences, int) else total_sequences)
                    col2.metric("Recombinant Events Detected", f"{total_recombinants:,}")
                    if isinstance(total_sequences, int) and total_sequences > 0:
                        col3.metric("Recombination Rate", f"{(total_recombinants / total_sequences):.2%}")

                    st.subheader("Extensive Summary of Recombinant Events")
                    st.dataframe(extensive_df)
                    
                    st.subheader("Distribution Plots")
                    if plots["recomb_plot"].exists():
                        st.image(str(plots["recomb_plot"]), caption="Distribution of Recombinant Events by Lineage")
                    if plots["date_plot"].exists():
                        st.image(str(plots["date_plot"]), caption="Temporal Distribution of Recombinant Events")
                    if plots["location_plot"].exists():
                        st.image(str(plots["location_plot"]), caption="Geographic Distribution of Recombinant Lineages")

                with tab2:
                    st.subheader("Recombinant Genome Explorer")
                    if not explorer_df.empty:
                        # --- Sidebar Filters ---
                        st.sidebar.header("Explorer Filters")
                        lineage_options = ["All"] + sorted(explorer_df['original_lineage'].unique())
                        selected_lineage = st.sidebar.selectbox("Filter by Original Lineage:", lineage_options)
                        
                        bp_options = ["All"] + sorted(explorer_df['breakpoint_count'].unique())
                        selected_bp = st.sidebar.selectbox("Filter by Breakpoint Count:", bp_options)
                        
                        # --- Filtering Logic ---
                        filtered_df = explorer_df.copy()
                        if selected_lineage != "All":
                            filtered_df = filtered_df[filtered_df['original_lineage'] == selected_lineage]
                        if selected_bp != "All":
                            filtered_df = filtered_df[filtered_df['breakpoint_count'] == selected_bp]
                        
                        st.dataframe(filtered_df[['genomeIDs', 'original_lineage', 'recombinant_parents', 'breakpoint_count', 'Collection date', 'Location']], hide_index=True)
                        
                        st.markdown("---")
                        st.subheader("View Case Details")
                        
                        # Create a list of options for the selectbox
                        options_list = [f"{row.genomeIDs.split(',')[0]} ({row.original_lineage})" for index, row in filtered_df.iterrows()]
                        
                        if not options_list:
                            st.info("No recombinants match the current filter criteria.")
                        else:
                            selected_case_display = st.selectbox("Select a recombinant to view details:", options=["-"] + options_list)
                            
                            if selected_case_display != "-":
                                selected_index = options_list.index(selected_case_display)
                                selected_row = filtered_df.iloc[selected_index]
                                
                                case_details = load_case_details(selected_row['case_report_folder'])
                                
                                if case_details.get('summary'):
                                    st.json(case_details['summary'])
                                if case_details.get('region_tables'):
                                    for i, region_df in enumerate(case_details['region_tables']):
                                        st.markdown(f"**Region {i+1} Table**")
                                        st.dataframe(region_df)
                                if case_details.get('plot_per_region'):
                                    st.plotly_chart(case_details['plot_per_region'], use_container_width=True)
                                    pass
                                if case_details.get('plot_whole_genome'):
                                    st.plotly_chart(case_details['plot_whole_genome'], use_container_width=True)
                    else:
                        st.info("No data available for the Recombinant Explorer.")
            else:
                st.error(f"No parameter sets found in the results directory for {selected_virus}.")
        else:
            st.error(f"No results directory found for {selected_virus}.")

if __name__ == "__main__":
    main()
