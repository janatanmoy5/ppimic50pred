import streamlit as st
import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors
import joblib
import os

st.set_page_config(page_title="Chemical Activity Predictor", layout="centered")
st.title("💊 Chemical Activity Predictor")

# ---------------- ML Model ---------------- #
MODEL_PATH = "ml_model.pkl"  # Replace with your trained model path

def load_model():
    if os.path.exists(MODEL_PATH):
        model = joblib.load(MODEL_PATH)
        st.success("✅ ML model loaded")
        return model
    else:
        st.warning("⚠️ ML model not found, predictions will be random")
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
    else:
        return None, None

def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    descriptors = {
        "MolecularWeight": Descriptors.MolWt(mol),
        "AlogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "NumHDonors": Descriptors.NumHDonors(mol)
    }
    return descriptors

# ---------------- Streamlit UI ---------------- #
chemical_name = st.text_input("Enter a chemical name", "")

if st.button("Search"):
    if chemical_name.strip() == "":
        st.warning("Please enter a chemical name")
    else:
        chembl_id, smiles = fetch_chembl_info(chemical_name.strip())
        if chembl_id is None:
            st.error("❌ Chemical not found in ChEMBL")
        else:
            st.success(f"Found: {chemical_name.upper()}")
            st.write(f"**ChEMBL ID:** {chembl_id}")
            st.write(f"**SMILES:** {smiles}")

            # Compute molecular descriptors
            desc = compute_descriptors(smiles)
            if desc:
                st.subheader("Molecular Properties (RDKit)")
                df = pd.DataFrame(desc, index=[0])
                st.dataframe(df)
            else:
                st.warning("Could not compute molecular descriptors")

            # Predict activity
            if model:
                try:
                    input_df = pd.DataFrame(desc, index=[0])
                    prediction = model.predict(input_df)[0]
                    st.subheader("Predicted Activity")
                    st.write(f"🔹 Activity Score: {prediction:.3f}")
                except Exception as e:
                    st.error(f"Prediction failed: {e}")
            else:
                st.info("Prediction skipped: no ML model loaded")
