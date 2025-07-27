# app/streamlit_app.py

import streamlit as st
import os
from pathlib import Path
import pandas as pd
import yaml
from streamlit_option_menu import option_menu
import re
import io
from bs4 import BeautifulSoup
import sys

# --- Path Setup for Imports ---
# Add the project root to the Python path to allow imports from 'src'
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.append(str(PROJECT_ROOT))

try:
    from src.utils.constants import *
except ImportError as e:
    # This is a fallback for st.error which might not be available yet
    print(f"CRITICAL ERROR: Could not import from 'src/utils/constants.py'. {e}")
    print("Please ensure the file exists and you are running the app from the project root.")
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
def find_latest_result_files(virus_name: str, param_set: str) -> dict:
    """Finds the most recent set of result files (txt, html, and figures)."""
    try:
        report_dir = RESULTS_DIR_BASE / "recombinhunt_output" / virus_name / param_set
        analysis_dir = RESULTS_DIR_BASE / "analysis_results" / virus_name / param_set
        
        txt_files = list(report_dir.glob("experiment-output.*.txt"))
        if not txt_files: return {}
        latest_txt = max(txt_files, key=lambda f: f.stat().st_mtime)
        
        base_name = latest_txt.name.replace("experiment-output.", "").replace(".txt", "")
        
        return {
            "txt_report": latest_txt,
            "html_report": report_dir / f"experiment-detail.{base_name}.html",
            "recomb_plot": analysis_dir / "distribution_by_lineage_recombinants_only.png",
            "date_plot": analysis_dir / "distribution_by_date.png",
            "location_plot": analysis_dir / "distribution_by_country_heatmap.png"
        }
    except FileNotFoundError:
        return {}

@st.cache_data
def parse_report_file(report_filepath: Path) -> (pd.DataFrame, dict):
    """Parses the summary table and recombinant detections from the report file."""
    summary_lines = []
    detection_lines = []
    in_summary = False
    in_detections = False
    
    summary_df = pd.DataFrame()
    recombinant_detections = {}

    try:
        with open(report_filepath, 'r') as f:
            for line in f:
                stripped = line.strip()
                if "Overall Summary Table" in stripped:
                    in_summary = True
                    in_detections = False
                    continue
                if "Recombinant Genomes Detected" in stripped:
                    in_summary = False
                    in_detections = True
                    continue
                
                if in_summary and stripped:
                    summary_lines.append(stripped)
                elif in_detections and stripped:
                    detection_lines.append(stripped)

        if summary_lines:
            table_string = "\n".join(summary_lines)
            summary_df = pd.read_csv(io.StringIO(table_string), sep=r'\s{2,}', engine='python')
            summary_df.columns = [col.strip() for col in summary_df.columns]

        if detection_lines:
            current_lineage = None
            for line in detection_lines:
                if line.startswith("Lineage:"):
                    current_lineage = line.split(":", 1)[1].strip()
                    recombinant_detections[current_lineage] = []
                elif current_lineage:
                    ids = [i.strip() for i in line.split(',')]
                    recombinant_detections[current_lineage].extend(ids)
    except Exception as e:
        st.error(f"Failed to parse report file {report_filepath}: {e}")
    
    return summary_df, recombinant_detections

@st.cache_data
def get_recombinant_details_df(detections: dict, html_filepath: Path, metadata_filepath: Path) -> pd.DataFrame:
    """
    Parses the HTML report to get details for each recombinant and merges with metadata.
    """
    if not html_filepath.exists():
        return pd.DataFrame()

    with open(html_filepath, 'r', encoding='utf-8') as f:
        soup = BeautifulSoup(f, 'html.parser')

    details_list = []
    for header in soup.find_all('h3'):
        header_text = header.get_text()
        # Corrected regex to capture the case name (genome ID) after the colon
        match = re.search(r'Case \d+ \(.*?\): (.*)', header_text)
        if match:
            genome_id = match.group(1).split(',')[0].strip()
            
            table = header.find_next('table')
            if table:
                # Corrected logic to find the cell containing "BC:"
                bc_td = table.find('td', string=re.compile(r'BC:'))
                if bc_td:
                    bc_string = bc_td.get_text()
                    parents_part = bc_string.split(':', 1)[1].strip()
                    bp_count = parents_part.count('+') # 1 '+' means 1BP (2 parents)
                    details_list.append({'genomeID': genome_id, 'BP_Count': f"{bp_count}BP"})

    if not details_list:
        return pd.DataFrame()

    details_df = pd.DataFrame(details_list)
    lineage_map_list = []
    for lineage, ids in detections.items():
        for genome_id in ids:
            lineage_map_list.append({'genomeID': genome_id, 'Detected Lineage': lineage})
    
    lineage_map_df = pd.DataFrame(lineage_map_list)

    try:
        metadata_df = pd.read_csv(metadata_filepath, sep='\t', low_memory=False)
        merged_df = pd.merge(lineage_map_df, metadata_df[['Virus name', 'Collection date', 'Location']], left_on='genomeID', right_on='Virus name', how='left')
    except Exception:
        merged_df = lineage_map_df
        merged_df['Collection date'] = 'N/A'; merged_df['Location'] = 'N/A'

    final_df = pd.merge(merged_df, details_df, on='genomeID', how='left')
    final_df['BP_Count'] = final_df['BP_Count'].fillna('N/A')
    return final_df

@st.cache_data
def extract_case_html(html_filepath: Path, genome_id: str) -> str:
    """Extracts the specific HTML block for a single case from the main report."""
    if not html_filepath.exists():
        return "<p>HTML report not found.</p>"
    
    with open(html_filepath, 'r', encoding='utf-8') as f:
        soup = BeautifulSoup(f, 'html.parser')
    
    # Corrected regex to find the header that ENDS with the genome ID
    header = soup.find('h3', string=re.compile(f': {re.escape(genome_id)}$'))
    if header:
        case_html = str(header)
        for sibling in header.find_next_siblings():
            if sibling.name == 'h3': break
            case_html += str(sibling)
        return case_html
    return f"<p>Details for {genome_id} not found in the HTML report.</p>"

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
                
                result_files = find_latest_result_files(selected_virus, latest_param_set)
                
                if result_files and result_files["txt_report"].exists():
                    summary_df, detections = parse_report_file(result_files["txt_report"])
                    
                    tab1, tab2, tab3 = st.tabs(["📊 Summary Dashboard", "🧬 Recombinant IDs", "🔬 Recombinant Explorer"])

                    with tab1:
                        st.subheader("Key Metrics")
                        if not summary_df.empty:
                            total_recombinants = summary_df['1BP'].sum() + summary_df['2BP'].sum()
                            total_sequences = summary_df['samples_run'].sum()
                            col1, col2, col3 = st.columns(3)
                            col1.metric("Total Sequences Analyzed", f"{total_sequences:,}")
                            col2.metric("Recombinant Events Detected", f"{total_recombinants:,}")
                            if total_sequences > 0:
                                col3.metric("Recombination Rate", f"{(total_recombinants / total_sequences):.2%}")
                        
                        st.subheader("Overall Summary")
                        if not summary_df.empty:
                            st.dataframe(summary_df)
                        else:
                            st.warning("Could not parse the summary table from the report file.")
                        
                        st.subheader("Distribution Plots")
                        if result_files["recomb_plot"].exists():
                            st.image(str(result_files["recomb_plot"]), caption="Distribution of Recombinant Events by Lineage")
                        if result_files["date_plot"].exists():
                            st.image(str(result_files["date_plot"]), caption="Temporal Distribution of Recombinant Events")
                        if result_files["location_plot"].exists():
                            st.image(str(result_files["location_plot"]), caption="Geographic Distribution of Recombinant Lineages")

                    with tab2:
                        st.subheader("Detected Recombinant Genome IDs")
                        if detections:
                            for lineage, ids in detections.items():
                                with st.expander(f"Lineage: {lineage} ({len(ids)} genomes)"):
                                    st.text_area("", ", ".join(ids), height=150)
                        else:
                            st.info("No recombinant genomes were listed in the report file.")

                    with tab3:
                        st.subheader("Recombinant Genome Explorer")
                        metadata_path = PROCESSED_DATA_DIR_BASE / selected_virus / "metadata.tsv"
                        recombinant_details_df = get_recombinant_details_df(detections, result_files["html_report"], metadata_path)

                        if not recombinant_details_df.empty:
                            recombinant_lineages = sorted(recombinant_details_df['Detected Lineage'].unique())
                            selected_lineage = st.selectbox("Select a lineage:", options=recombinant_lineages)
                            
                            if selected_lineage:
                                st.markdown(f"**Recombinant Genomes for Lineage: {selected_lineage}**")
                                lineage_specific_df = recombinant_details_df[recombinant_details_df['Detected Lineage'] == selected_lineage]
                                st.dataframe(lineage_specific_df[['genomeID', 'BP_Count', 'Collection date', 'Location']], hide_index=True)
                                
                                genome_to_view = st.selectbox("Select a Genome ID to view its detailed report:", options=lineage_specific_df['genomeID'])
                                
                                if genome_to_view:
                                    case_html = extract_case_html(result_files["html_report"], genome_to_view)
                                    st.subheader(f"Detailed Report for: {genome_to_view}")
                                    st.components.v1.html(case_html, height=800, scrolling=True)
                        else:
                            st.info("No recombinant details could be parsed from the HTML report for interactive exploration.")
                else:
                    st.error(f"No result files found for {selected_virus}.")
            else:
                st.error(f"No parameter sets found for {selected_virus}.")
        else:
            st.error(f"No results directory found for {selected_virus}.")

if __name__ == "__main__":
    main()
