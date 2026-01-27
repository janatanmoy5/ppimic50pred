# run.py
import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from chembl_webresource_client.new_client import new_client
from sklearn.ensemble import RandomForestRegressor
import joblib
import os

st.set_page_config(
    page_title="PPI MIC50 Predictor",
    layout="wide"
)

st.title("PPI MIC50 Prediction")

# ---------------- CONFIG ---------------- #
MODEL_PATH = os.path.join(os.path.dirname(__file__), "random_forest_model.pkl")

# ---------------- Local fallback compounds ---------------- #
LOCAL_COMPOUNDS = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
}

# ---------------- Functions ---------------- #
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    desc = {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol)
    }
    return pd.DataFrame([desc])

def load_model():
    if os.path.exists(MODEL_PATH):
        return joblib.load(MODEL_PATH)
    return None

def predict_mic50(model, descriptors_df):
    if model is None or descriptors_df is None:
        return None
    return model.predict(descriptors_df)[0]

# ---------------- UI ---------------- #
compound_name = st.selectbox(
    "Select a compound:",
    options=list(LOCAL_COMPOUNDS.keys())
)

smiles_input = st.text_input(
    "Or enter a SMILES string manually:",
    value=LOCAL_COMPOUNDS.get(compound_name, "")
)

if st.button("Predict MIC50"):
    if not smiles_input:
        st.error("Please provide a SMILES string.")
    else:
        descriptors_df = compute_descriptors(smiles_input)
        if descriptors_df is None:
            st.error("Invalid SMILES string.")
        else:
            model = load_model()
            mic50_pred = predict_mic50(model, descriptors_df)
            if mic50_pred is not None:
                st.success(f"Predicted MIC50: {mic50_pred:.3f} μM")
                st.write("Molecular Descriptors:")
                st.dataframe(descriptors_df)
            else:
                st.warning("Model not found. Please upload your trained model as 'random_forest_model.pkl'.")
