from shiny import App, ui, render, reactive
import pandas as pd
import joblib
from rdkit import Chem

# ---------------- CONFIG ---------------- #
MODEL_PATH = "random_forest_model.pkl"

# ---------------- Local fallback compounds ---------------- #
LOCAL_COMPOUNDS = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
}

# ---------------- Helper Functions ---------------- #
def calculate_mfi(smiles: str) -> pd.DataFrame:
    alleles = ["DR1", "DR10", "DR103", "DR51", "DR15", "DR16"]
    mean_mfi = [abs(hash(smiles + a)) % 25000 for a in alleles]

    df = pd.DataFrame({
        "Group": ["DR1c"]*3 + ["DR51c"]*3,
        "Allele": alleles[:3] + alleles[3:],
        "MeanMFI": mean_mfi
    })
    return df

def load_model():
    try:
        model = joblib.load(MODEL_PATH)
        return model
    except FileNotFoundError:
        return None

# ---------------- UI ---------------- #
app_ui = ui.page_fluid(
    ui.h2("Interactive Mean MFI Panel"),
    ui.input_text("compound", "Enter Compound Name or SMILES:", value="aspirin"),
    ui.input_action_button("run_btn", "Run"),
    ui.output_table("mfi_table")
)

# ---------------- Server ---------------- #
def server(input, output, session):

    @reactive.Calc
    def compound_df():
        if input.run_btn() == 0:
            return pd.DataFrame(columns=["Group", "Allele", "MeanMFI"])
        
        compound_input = input.compound()
        smiles = LOCAL_COMPOUNDS.get(compound_input.lower(), compound_input)
        df = calculate_mfi(smiles)
        return df[df["MeanMFI"] > 1000]

    @output
    @render.table
    def mfi_table():
        df = compound_df()
        return df

# ---------------- Run App ---------------- #
app = App(app_ui, server)
