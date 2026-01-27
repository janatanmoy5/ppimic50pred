import streamlit as st
import math
import requests
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Draw
from chembl_webresource_client.new_client import new_client
from sklearn.metrics import r2_score
import joblib
import warnings
import os
import time

warnings.filterwarnings("ignore")
RDLogger.DisableLog('rdApp.*')

# ---------------- CONFIG ---------------- #
MODEL_PATH = os.path.join(os.path.dirname(__file__), "random_forest_model.pkl")

# ---------------- Local Fallback Compounds ---------------- #
LOCAL_COMPOUNDS = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
}

# ---------------- Utility Functions ---------------- #
def is_chembl_id(s):
    return s.upper().startswith("CHEMBL")

def is_smiles(s):
    return Chem.MolFromSmiles(s) is not None

def safe_get(url, retries=3, timeout=10):
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=timeout)
            r.raise_for_status()
            return r
        except requests.RequestException:
            if attempt < retries - 1:
                time.sleep(2)
            else:
                return None

def search_chembl_by_name(name):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={name}&format=json"
    r = safe_get(url)
    if not r:
        return None, None
    molecules = r.json().get("molecules", [])
    if molecules:
        mol = molecules[0]
        cid = mol.get("molecule_chembl_id")
        smi = mol.get("molecule_structures", {}).get("canonical_smiles")
        return cid, smi
    return None, None

def get_smiles_from_chembl(chembl_id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    r = safe_get(url)
    if not r:
        return None
    return r.json().get("molecule_structures", {}).get("canonical_smiles")

def get_chembl_id_from_smiles_similarity(smiles, threshold=0.95):
    url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{smiles}/{threshold}?format=json"
    r = safe_get(url)
    if not r:
        return None
    mols = r.json().get('molecules', [])
    if mols:
        return mols[0].get('molecule_chembl_id')
    return None

def compute_rdkit_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {name: func(mol) for name, func in Descriptors.descList}

def process_input(user_input):
    user_input = user_input.strip().replace(" ", "")
    chembl_id, smiles, compound_name = None, None, None

    # Local fallback
    if user_input.lower() in LOCAL_COMPOUNDS:
        smiles = LOCAL_COMPOUNDS[user_input.lower()]
        compound_name = user_input
    elif is_chembl_id(user_input):
        chembl_id = user_input.upper()
        smiles = get_smiles_from_chembl(chembl_id)
    elif is_smiles(user_input):
        smiles = user_input
        chembl_id = get_chembl_id_from_smiles_similarity(smiles)
    else:
        chembl_id, smiles = search_chembl_by_name(user_input)
        compound_name = user_input

    if not smiles:
        st.error(f"Could not resolve '{user_input}' to a valid SMILES.")
        return None, None, None

    descriptors = compute_rdkit_descriptors(smiles)
    if descriptors is None:
        st.error("Descriptor calculation failed.")
        return None, None, None

    descriptors["chembl_id"] = chembl_id
    descriptors["smiles"] = smiles
    return descriptors, chembl_id, compound_name

# ---------------- Prediction ---------------- #
def predict_ic50(descriptor_dict, model_path):
    df = pd.DataFrame([descriptor_dict])
    drop_cols = ['chembl_id', 'smiles']
    df.drop(columns=[c for c in drop_cols if c in df.columns], inplace=True)

    saved = joblib.load(model_path)
    model, scaler = saved["model"], saved["scaler"]
    X_scaled = scaler.transform(df)
    log_pred = model.predict(X_scaled)[0]

    confidence = None
    if "r2" in saved:
        confidence = saved["r2"]
    elif "X_train" in saved and "y_train" in saved:
        X_train_scaled = scaler.transform(saved["X_train"])
        y_train_pred = model.predict(X_train_scaled)
        confidence = r2_score(saved["y_train"], y_train_pred)

    return log_pred, confidence

# ---------------- Streamlit Interface ---------------- #
st.set_page_config(page_title="PPIM-IC50Pred", layout="wide")
st.title("⚗️ PPIM-IC50Pred Webserver")

user_input = st.text_input("Enter chemical name, ChEMBL ID, or SMILES:")

if user_input:
    descriptors, chembl_id, compound_name = process_input(user_input)

    if descriptors:
        st.subheader("Compound Details")
        st.markdown(f"**Compound Name:** {compound_name or 'Unknown'}")
        st.markdown(f"**ChEMBL ID:** {chembl_id or 'N/A'}")
        st.markdown(f"**SMILES:** {descriptors['smiles']}")

        # 2D Structure
        st.subheader("2D Structure")
        mol = Chem.MolFromSmiles(descriptors['smiles'])
        if mol:
            img = Draw.MolToImage(mol, size=(300, 300))
            st.image(img, caption=compound_name or chembl_id)
        else:
            st.warning("Could not render molecular structure.")

        # Prediction
        st.subheader("Predicted log(IC50) [nM]")
        log_val, conf = predict_ic50(descriptors, MODEL_PATH)
        st.markdown(f"**Predicted:** {log_val:.4f}")
        if conf is not None:
            st.markdown(f"**Model Confidence (R²):** {conf:.4f}")

