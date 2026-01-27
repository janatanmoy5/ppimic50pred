# run.py - Shiny for Python App for chemical input and MFI/prediction
# Requires: shiny>=0.11, pandas>=2.1.1, requests>=2.32.0

from shiny import App, ui, reactive, render
import pandas as pd
import requests

# ------------------------------
# Local fallback compounds
# ------------------------------
LOCAL_COMPOUNDS = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
}

# ------------------------------
# Shiny UI
# ------------------------------
app_ui = ui.page_fluid(
    ui.h2("Chemical MFI / Prediction Panel"),
    ui.input_text("chemical_name", "Enter Chemical Name:", value="aspirin"),
    ui.input_action_button("run_btn", "Run"),
    ui.output_text_verbatim("status_text"),
    ui.output_table("result_table")
)

# ------------------------------
# Shiny Server
# ------------------------------
def server(input, output, session):

    # Reactive value to store results
    results_df = reactive.Value(pd.DataFrame())

    @input.run_btn
    def fetch_results():
        chem_name = input.chemical_name()
        if not chem_name:
            results_df.set(pd.DataFrame())
            return

        chem_name_lower = chem_name.lower()
        # Lookup local compound
        smiles = LOCAL_COMPOUNDS.get(chem_name_lower)
        if smiles:
            # Dummy MFI prediction values for demonstration
            df = pd.DataFrame({
                "Chemical": [chem_name],
                "SMILES": [smiles],
                "Predicted_MFI": [round(len(smiles) * 1234)]  # just an example formula
            })
        else:
            # If chemical not found, fallback
            df = pd.DataFrame({
                "Chemical": [chem_name],
                "SMILES": ["Not Found"],
                "Predicted_MFI": [None]
            })
        results_df.set(df)

    # Status / messages
    @output
    @render.text
    def status_text():
        df = results_df.get()
        if df.empty:
            return "Enter a chemical name and click Run."
        return f"Showing results for {input.chemical_name()}"

    # Result table
    @output
    @render.table
    def result_table():
        return results_df.get()


# ------------------------------
# Create Shiny App
# ------------------------------
app = App(app_ui, server)
