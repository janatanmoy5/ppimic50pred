import subprocess
import sys
import os

# ---------------- Auto-install packages ---------------- #
def install_package(pkg, import_name=None):
    import_name = import_name or pkg
    try:
        __import__(import_name)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])

# Install required packages
for pkg, imp_name in [
    ("streamlit", "streamlit"),
    ("pandas", "pandas"),
    ("rdkit-pypi", "rdkit"),
    ("scikit-learn", "sklearn"),
    ("joblib", "joblib"),
    ("chembl_webresource_client", "chembl_webresource_client")
]:
    install_package(pkg, imp_name)

# ---------------- Imports ---------------- #
import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
import joblib
from chembl_webresource_client.new_client import new_client

# ---------------- Page config ---------------- #
st.set_page_config(
    page_title="🔬 PPI MIC50 Predictor",
    layout="wide"
)

st.title("🔬 PPI MIC50 Prediction App")
st.markdown("Enter a chemical name to fetch its structure and predict MIC50.")

# ---------------- Sidebar ---------------- #
st.sidebar.header("Chemical Input")
chemical_name = st.sidebar.text_input("Enter chemical name")
search_button = st.sidebar.button("Search")

# ---------------- Load model ---------------- #
MODEL_PATH = "random_forest_model.pkl"
if os.path.exists(MODEL_PATH):
    model = joblib.load(MODEL_PATH)
else:
    model = None
    st.warning("⚠️ Random Forest model not found. Predictions will not run.")

# ---------------- Descriptor function ---------------- #
def calculate_descriptors(mol):
    """Simple RDKit descriptors for demonstration"""
    return {
        "MolWt": Descriptors.MolWt(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumRotatableBonds": Descriptors.NumRotatableBonds(mol)
    }

# ---------------- Prediction function ---------------- #
def predict_activity(descriptors):
    """Predict MIC50 using the loaded model"""
    if model is None:
        return None
    df = pd.DataFrame([descriptors])
    pred = model.predict(df)
    return pred[0]

# ---------------- Main App ---------------- #
if search_button:
    if chemical_name.strip() == "":
        st.error("Please enter a chemical name!")
    else:
        st.info(f"Searching for **{chemical_name}** in ChEMBL...")
        molecule = new_client.molecule
        res = molecule.filter(pref_name__iexact=chemical_name)
        if res:
            mol_data = res[0]  # Take first match
            smiles = mol_data.get("molecule_structures", {}).get("canonical_smiles")
            if smiles:
                st.success(f"✅ Found SMILES: `{smiles}`")
                mol = Chem.MolFromSmiles(smiles)
                st.image(Draw.MolToImage(mol, size=(300, 300)))

                # Calculate descriptors
                desc = calculate_descriptors(mol)
                st.write("**Descriptors:**")
                st.json(desc)

                # Predict MIC50
                pred = predict_activity(desc)
                if pred is not None:
                    st.metric("Predicted MIC50", f"{pred:.4f}")
                else:
                    st.info("Model not loaded. Cannot predict MIC50.")
            else:
                st.error("SMILES not found for this chemical.")
        else:
            st.error("Chemical not found in ChEMBL database.")
else:
    st.info("Enter a chemical name and click Search to start.")
