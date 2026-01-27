import streamlit as st
import pandas as pd
import subprocess
import sys

# ---------------- Install missing packages dynamically ---------------- #
def install(package):
    """Install package via pip if not already installed."""
    try:
        __import__(package)
    except ModuleNotFoundError:
        st.info(f"Installing {package}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# ---------------- Install required packages ---------------- #
for pkg in ["rdkit-pypi", "chembl_webresource_client", "pandas", "scikit-learn"]:
    install(pkg)

# ---------------- Now import after install ---------------- #
from rdkit import Chem
from rdkit.Chem import Descriptors
from chembl_webresource_client.new_client import new_client
import joblib
import os

# ---------------- Streamlit App ---------------- #
st.set_page_config(page_title="Chemical Activity Predictor", layout="centered")
st.title("💊 Chemical Activity Predictor")

MODEL_PATH = "ml_model.pkl"  # Path to your ML model

# Load ML model if exists
def load_model():
    if os.path.exists(MODEL_PATH):
        model = joblib.load(MODEL_PATH)
        st.success("✅ ML model loaded")
        return model
    else:
        st.warning("⚠️ ML model not found, predictions will be skipped")
        return None

model = load_model()

# ---------------- Helper Functions ---------------- #
def fetch_chembl_info(chemical_name):
    molecule = new_client.molecule
    res = molecule.filter(pref_name__iexact=chemical_name).only(
        ["molecule_chembl_id", "pref_name", "molecule_structures"]
    )
    if res:
        chembl_id = res[0]["molecule_chembl_id"]
        smiles = res[0]["molecule_structures"]["canonical_smiles"]
        return chembl_id, smiles
    return None, None

def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MolecularWeight": Descriptors.MolWt(mol),
        "AlogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
    }

# ---------------- Streamlit UI ---------------- #
chemical_name = st.text_input("Enter a chemical name", "")

if st.button("Search"):
    if not chemical_name.strip():
        st.warning("Please enter a chemical name")
    else:
        chembl_id, smiles = fetch_chembl_info(chemical_name.strip())
        if chembl_id is None:
            st.error("❌ Chemical not found in ChEMBL")
        else:
            st.success(f"Found: {chemical_name.upper()}")
            st.write(f"**ChEMBL ID:** {chembl_id}")
            st.write(f"**SMILES:** {smiles}")

            descriptors = compute_descriptors(smiles)
            if descriptors:
                st.subheader("Molecular Properties (RDKit)")
                st.dataframe(pd.DataFrame(descriptors, index=[0]))
            else:
                st.warning("Could not compute molecular descriptors")

            if model and descriptors:
                try:
                    df_input = pd.DataFrame(descriptors, index=[0])
                    prediction = model.predict(df_input)[0]
                    st.subheader("Predicted Activity")
                    st.write(f"🔹 Activity Score: {prediction:.3f}")
                except Exception as e:
                    st.error(f"Prediction failed: {e}")
            else:
                st.info("Prediction skipped: no ML model or descriptors available")
