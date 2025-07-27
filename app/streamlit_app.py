# app/streamlit_app.py

import streamlit as st

def main():
    """
    Main function to run the Streamlit application.
    """
    # st.set_page_config is used to configure the page's title, icon, etc.
    st.set_page_config(
        page_title="OpenRecombinHunt Results",
        page_icon="🧬",
        layout="wide"
    )

    # st.title() creates the main title of the web page.
    st.title("Welcome to the OpenRecombinHunt Pipeline Results")

    # st.write() can display text, DataFrames, plots, and more.
    st.write("This web application will display the results from the automated recombination analysis pipeline.")
    st.info("We will build this out step-by-step to show the analysis for each virus.")

if __name__ == "__main__":
    main()
