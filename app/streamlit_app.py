import sys
from pathlib import Path
import yaml
import pandas as pd
import streamlit as st
from streamlit_option_menu import option_menu
import plotly.express as px
import plotly.graph_objects as go
from geopy.geocoders import Nominatim

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
except FileNotFoundError:
    st.error(f"Configuration file not found at {CONFIG_PATH}. Please ensure the file exists.")
    st.stop()
except Exception as e:
    st.error(f"Error loading or parsing configuration file: {e}")
    st.stop()

def visualize(virus):
    mapping = {
        "sars-cov-2": "SARS-CoV-2",
        "yellow-fever": "Yellow Fever",
        "zika": "Zika",
        "rsv-a": "RSV-A",
        "rsv-b": "RSV-B",
        "monkeypox": "Monkeypox",
        "influenza": "Influenza (H5N1)"
    }
    return mapping.get(virus.lower())

@st.cache_data
def discover_viruses():
    """discover available viruses by scanning recombinhunt_output/ directory"""
    recombinhunt_output_path = RESULTS_DIR_BASE / RECOMBINHUNT_OUTPUT

    if not recombinhunt_output_path.exists():
        st.warning(f"No recombinhunt output found at {recombinhunt_output_path}.")
        return []

    virus_dirs = [d for d in recombinhunt_output_path.iterdir() if d.is_dir()]
    if not virus_dirs:
        st.warning(f"No virus directories found in {recombinhunt_output_path}.")
        return []

    # return only the names, not paths
    return [d.name for d in virus_dirs]

@st.cache_data
def load_master_data(virus):
    """load and merge master data for the specified virus"""

    # load param set
    try:
        params = config.get(VIRUSES).get(virus).get(PARAMETERS).get(HAPLOCOV)
        dist = params.get(DIST)
        size = params.get(SIZE)
        paramset = f"dist{dist}size{size}"
    except Exception as e:
        st.error(f"Error loading parameters for {virus}: {e}")
        return None

    # load source data
    try:
        if virus.lower() == "sars-cov-2":
            source_file = RESULTS_DIR_BASE / NEXTSTRAIN_OUTPUT / virus / "nextstrain_reformatted.tsv"
            columns = ["genomeID", "Collection date", "Submission date", "Location", "pangoLin"]
        else:
            source_file = RESULTS_DIR_BASE / HAPLOCOV_OUTPUT / virus / paramset / "haplocov_reformatted.tsv"
            columns = ["genomeID", "collectionD", "continent", "country", "pangoLin"]

        if source_file.exists():
            source_df = pd.read_csv(source_file, sep="\t", usecols=columns)
        else:
            st.warning(f"Source file not found: {source_file}")
    except Exception as e:
        st.error(f"Error loading master data for {virus}: {e}")
        return None
    
    # load recombinant summary
    try:
        recombinant_summary_file = RESULTS_DIR_BASE / RECOMBINHUNT_OUTPUT / virus / paramset / "recombinant_summary.tsv"
        if recombinant_summary_file.exists():
            recombinant_summary_df = pd.read_csv(recombinant_summary_file, sep="\t")
            recombinant_summary_df.rename(columns={"genomeIDs": "genomeID"}, inplace=True)
        else:
            st.warning(f"Recombinant summary file not found: {recombinant_summary_file}")
    except Exception as e:
        st.error(f"Error loading recombinant summary for {virus}: {e}")
        return None

    # merge dataframes
    try:
        merged_df = pd.merge(source_df, recombinant_summary_df, on="genomeID", how="left")
    except Exception as e:
        st.error(f"Error merging dataframes for {virus}: {e}")
        return None

    # add is_recombinant column to merged_df
    # if breakpoint_count.notnull()
    merged_df["is_recombinant"] = merged_df["breakpoint_count"].notnull()

    # rename some columns of merged_df
    mapping = {
        "collectionD": "collection_date",
        "Collection date": "collection_date",
        "Submission date": "submission_date",
    }
    merged_df.rename(columns=mapping, inplace=True)

    return merged_df

def apply_time_filter(df):
    """
    applies time based filtering to the dataframe
    user can select between three options:
        -no time filtering (returns all data)
        -filter by specific date selections (choose start and end dates: end date defaults to today)
        -filter by the number of latest X sequences (user inputs an integer number, we sort the dataframe by collection date and return the latest X rows)
    """
    # select filter type
    filter_type = st.pills(
        "Choose filter type:",
        options = ["No filtering", "Filter by Date Range", "Filter by Latest X Sequences"],
        selection_mode = "single"
    )

    # select filter value
    filter_value = None
    if filter_type == "Filter by Date Range":
        col1, col2 = st.columns(2)
        with col1:
            start_date = st.date_input("Start Date", value=pd.to_datetime("2020-01-01"))
        with col2:
            end_date = st.date_input("End Date", value=pd.to_datetime("today"))
        filter_value = (start_date, end_date)
    elif filter_type == "Filter by Latest X Sequences":
        filter_value = st.number_input("Number of Latest Sequences", min_value=1, value=10)

    # apply filtering to the df
    filtered_df = df
    if filter_type == "Filter by Date Range" and filter_value:
        start_date, end_date = filter_value
        # convert string df['collection_date'] to datetime object
        collection_date = pd.to_datetime(df["collection_date"]).dt.date

        filtered_df = df[(collection_date >= start_date) & (collection_date <= end_date)]
    elif filter_type == "Filter by Latest X Sequences" and filter_value:
        filtered_df = df.sort_values("collection_date", ascending=False).head(filter_value)

    return filtered_df

def create_key_metrics(df):
    """
    creates key metrics as cards for the summary dashboard
        -# Total Sequences
        -# Recombinant Events
        -Recombination Rate
        -Top Recombinant Lineage
        -Most Common Parents
    """
    st.title("🔑 Key Metrics")

    total_sequences = len(df)
    total_recombinants = df["is_recombinant"].sum()
    num_1BP = df[df["breakpoint_count"] == "1BP"].shape[0]
    num_2BP = df[df["breakpoint_count"] == "2BP"].shape[0]
    recombination_rate = (total_recombinants / total_sequences * 100) if total_sequences > 0 else 0
    top_recombinant_lineage = df[df["is_recombinant"]]["pangoLin"].mode().values[0] if total_recombinants > 0 else "N/A"
    most_common_parents = df[df["is_recombinant"]]["recombinant_parents"].value_counts().idxmax() if total_recombinants > 0 else "N/A"

    a, b = st.columns(2)
    c, d, e = st.columns(3)
    [f] = st.columns(1) 
    [g] = st.columns(1)

    a.metric("Total Sequences", total_sequences, border=True)
    b.metric("Recombination Rate", f"{recombination_rate:.2f}%", border=True)

    c.metric("Recombination Events", total_recombinants, border=True)
    d.metric("1BP", num_1BP, border=True)
    e.metric("2BP", num_2BP, border=True)

    f.metric("Top Recombinant Lineage", top_recombinant_lineage, border=True)

    g.metric("Most Common Parents", most_common_parents, border=True)

def create_summary_tables(df):
    "create summary and hotspots tables"
    st.title("📋 Summary Tables")

    st.subheader("Lineage Breakdown")
    with st.expander("Show Lineage Breakdown"):
        # create a table
        # group by each lineage (pangoLin)
        # count 1BP in breakpoint_count as 1BP Count
        # rate of 1BP Count / total sequences
        # count 2BP in breakpoint_count as 2BP Count
        # rate of 2BP Count / total sequences
        # rest is No Recombination
        # total sequences in each lineage
        lineage_breakdown = df.groupby("pangoLin").agg(
            BP1_Count=("breakpoint_count", lambda x: (x == "1BP").sum()),
            BP1_Rate=("breakpoint_count", lambda x: f"{(x == '1BP').sum() / len(x) * 100:.2f}%" if len(x) > 0 else "0.00%"),
            BP2_Count=("breakpoint_count", lambda x: (x == "2BP").sum()),
            BP2_Rate=("breakpoint_count", lambda x: f"{(x == '2BP').sum() / len(x) * 100:.2f}%" if len(x) > 0 else "0.00%"),
            No_Recombination=("breakpoint_count", lambda x: ((x != "1BP") & (x != "2BP")).sum()),
            Total_Sequences=("breakpoint_count", "size")
        ).reset_index()

        lineage_breakdown.rename(columns={
            "pangoLin": "Lineage",
            "BP1_Count": "1BP Count",
            "BP1_Rate": "1BP Rate",
            "BP2_Count": "2BP Count",
            "BP2_Rate": "2BP Rate",
            "No_Recombination": "No Recombination",
            "Total_Sequences": "Total Sequences"
        }, inplace=True)

        lineage_breakdown.set_index("Lineage", inplace=True)

        st.write(lineage_breakdown)

    st.subheader("Recombination Hotspots")
    with st.expander("Show Recombination Hotspots"):
        # group by recombinant_parents
        # display frequency
        recombination_hotspots = df.groupby("recombinant_parents").size().reset_index(name="Frequency")
        recombination_hotspots.rename(columns={"recombinant_parents": "Recombinant Parents"}, inplace=True)
        recombination_hotspots.set_index("Recombinant Parents", inplace=True)
        recombination_hotspots.sort_values(by="Frequency", ascending=False, inplace=True)
        st.write(recombination_hotspots)

def create_temporal_plot(df):
    if "collection_date" not in df.columns:
        st.warning("Collection date information is not available.")
        return

    df = df.dropna(subset=["collection_date"])
    df["collection_date"] = pd.to_datetime(df["collection_date"])

    date_range = df["collection_date"].max() - df["collection_date"].min()

    if date_range > pd.Timedelta(days=180):
        freq = "M" 
        freq_label = "Month"  
    else:
        freq = "W" 
        freq_label = "Week" 

    df["year-month"] = df["collection_date"].dt.to_period(freq)

    monthly_data = df.groupby("year-month").agg(
        {
            "is_recombinant": ["count", "sum"]
        }
    )

    monthly_data.columns = [
        "total_sequences", "recombinations"
    ]

    monthly_data = monthly_data.reset_index()
    monthly_data["year-month"] = monthly_data["year-month"].astype(str)

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=monthly_data["year-month"],
            y=monthly_data["total_sequences"],
            fill="tonexty",
            mode="none",
            name=f"Total Sequences per {freq_label}",
            fillcolor="rgba(74, 144, 226, 0.15)",
            line=dict(width=0)
        )
    )

    fig.add_trace(
        go.Scatter(
            x=monthly_data["year-month"],
            y=monthly_data["recombinations"],
            mode="lines+markers",
            name="Recombination Events",
            line=dict(color="#4A90E2",
                             width=3,
                             shape="spline",
                             smoothing=1.3),
            marker=dict(size=6, color="#4A90E2")
        )
    )

    fig.update_layout(
        title="Temporal Distribution of Recombination Events",
        xaxis_title="Time",
        yaxis_title="Number of Sequences",
        hovermode="x unified",
        height=500,
        showlegend=True
    )

    st.plotly_chart(fig, width="strecth")

def create_geographic_map(df, virus):
    # only keep is_recombinant
    df = df[df["is_recombinant"] == True]
    if df.empty:
        st.warning("No recombinant sequences found.")
        return

    if virus == "sars-cov-2":
        df['country'] = df["Location"].apply(lambda x: x.split("/")[1].strip() if isinstance(x, str) else x)

    df['country'] = df['country'].str.strip()
    geo_data = df["country"].value_counts().reset_index()
    geo_data.columns = ["country", "count"]

    lat_lon_df = pd.read_csv("app/country.csv")

    geo_data = geo_data.merge(lat_lon_df, on="country", how="left")
    
    missing_countries = geo_data[geo_data["latitude"].isna()]["country"].tolist()

    if missing_countries:
        print("Missing countries (will geocode):", missing_countries)

        for country in missing_countries:
            try:
                loc = Nominatim(user_agent="GetLoc")
                getLoc = loc.geocode(country)

                geo_data.loc[geo_data["country"] == country, "latitude"] = getLoc.latitude
                geo_data.loc[geo_data["country"] == country, "longitude"] = getLoc.longitude
            except Exception as e:
                print(f"Could not geocode {country}: {e}")

    geo_data = geo_data.dropna(subset=["latitude", "longitude"])

    fig = px.scatter_mapbox(
        geo_data,
        lat="latitude",
        lon="longitude",
        size="count",
        hover_name="country",
        color="count",
        color_continuous_scale=px.colors.sequential.Blues,
        size_max=40,
        zoom=1
    )

    fig.update_layout(mapbox_style="carto-positron")
    fig.update_layout(margin={"r":0,"t":0,"l":0,"b":0})

    st.plotly_chart(fig, use_container_width=True)


def create_distribution_plots(df, virus):
    "creates temporal and locational distributions"
    st.title("📈 Distribution Plots")

    st.subheader("🕗 Temporal Distribution")
    create_temporal_plot(df)

    st.subheader("🌏 Locational Distribution")
    create_geographic_map(df, virus)

def sidebar(virus_list):
    """sidebar navigation for the streamlit"""
    with st.sidebar:

        menu_options = ["Home"] + sorted([visualize(v) for v in virus_list])
        menu_icons = ["🏠"] + ["🦠"] * len(virus_list)

        selected = option_menu(
            menu_title="🧬 OpenRecombinHunt",
            options=menu_options,
            icons=menu_icons,
            menu_icon="cast",
            default_index=0,
            orientation="vertical",
            styles={
                "container": {"padding": "0!important", "background-color": "#fafafa"},
                "icon": {"color": "#4A90E2", "font-size": "18px"},
                "nav-link": {"font-size": "16px", "text-align": "left", "margin": "0px", "--hover-color": "#eee"},
                "nav-link-selected": {"background-color": "#4A90E2"},
            }
        )

        return selected

def show_home_page():
    """display welcome page"""
    st.title("🧬 OpenRecombinHunt Dashboard")
    st.markdown("---")

    st.markdown("""
    ## Welcome to the OpenRecombinHunt Bioinformatics Dashboard
    
    This sophisticated multi-page dashboard provides comprehensive analysis and visualization 
    capabilities for viral recombination detection across multiple virus species.
    
    ### Features:
    
    📊 **Summary Dashboard**
    - Dynamic time-based filtering
    - Key performance metrics
    - Interactive temporal visualizations
    - Geographic distribution analysis
    - Comprehensive summary tables
    
    🔬 **Recombinant Explorer**
    - Advanced filtering capabilities
    - Interactive case selection
    - Detailed report exploration
    - On-demand data loading
    
    ### Getting Started:
    1. Select a virus from the sidebar menu
    2. Choose between Summary Dashboard or Recombinant Explorer tabs
    3. Use filters to customize your analysis
    4. Explore detailed cases in the Recombinant Explorer
    
    ### Data Sources:
    - **SARS-CoV-2**: Nextstrain reformatted data
    - **Other Viruses**: HaploCov reformatted data
    - **Recombination Analysis**: RecombinHunt output files
    
    ---
    *This dashboard is part of a master's thesis project focused on advancing 
    bioinformatics analysis capabilities for viral recombination detection.*
    """)

def show_virus_page(virus):
    """display virus-specific analysis and visualizations"""
    st.title(f"🦠 {visualize(virus)} Dashboard")

    # load master data
    with st.spinner(f"Loading data for {virus}..."):
        master_df = load_master_data(virus)

    if master_df is None or master_df.empty:
        st.error(f"No data found for {virus}.")
        return

    # tabs
    tab1, tab2 = st.tabs(["📊 Summary Dashboard", "🔬 Recombinant Explorer"])

    with tab1:
        # time-based filtering
        summary_df = apply_time_filter(master_df)

        st.markdown("---")

        # create key metrics and display
        create_key_metrics(summary_df)

        st.markdown("---")

        # create summary tables and display
        create_summary_tables(summary_df)

        st.markdown("---")

        create_distribution_plots(summary_df, virus)

def main():

    # discover available viruses
    viruses = discover_viruses()

    viruses_visualized = [visualize(v) for v in viruses]

    # sidebar navigation
    selected = sidebar(viruses)

    if selected == "Home":
        show_home_page()
    elif selected in viruses_visualized:
        show_virus_page(viruses[viruses_visualized.index(selected)])

if __name__ == "__main__":
    main()