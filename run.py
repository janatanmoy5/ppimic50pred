import streamlit as st
import math
import requests
import pandas as pd
import joblib
import warnings
import os
import time

from sklearn.metrics import r2_score
from sklearn.exceptions import InconsistentVersionWarning

warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

# ---------------- RDKit ---------------- #
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, Draw
    RDLogger.DisableLog("rdApp.*")
    RDKit_AVAILABLE = True
except ImportError:
    RDKit_AVAILABLE = False

# ---------------- CONFIG ---------------- #
BASE_DIR = os.path.dirname(__file__)
MODEL_PATH = os.path.join(BASE_DIR, "random_forest_model.pkl")

# ---------------- Local Fallback Compounds ---------------- #
LOCAL_COMPOUNDS = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
}

# ---------------- Utility Functions ---------------- #
def is_chembl_id(x: str) -> bool:
    return x.upper().startswith("CHEMBL")

def is_valid_smiles(smiles: str) -> bool:
    if not RDKit_AVAILABLE:
        return False
    return Chem.MolFromSmiles(smiles) is not None

def safe_get(url, retries=3, timeout=10):
    for _ in range(retries):
        try:
            r = requests.get(url, timeout=timeout)
            r.raise_for_status()
            return r
        except requests.RequestException:
            time.sleep(1)
    return None

# ---------------- ChEMBL Resolution ---------------- #
def chembl_name_to_smiles(name):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={name}&format=json"
    r = safe_get(url)
    if not r:
        return None, None
    mols = r.json().get("molecules", [])
    if not mols:
        return None, None
    m = mols[0]
    return m.get("molecule_chembl_id"), m.get("molecule_structures", {}).get("canonical_smiles")

def chembl_id_to_smiles(chembl_id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    r = safe_get(url)
    if not r:
        return None
    return r.json().get("molecule_structures", {}).get("canonical_smiles")

# ---------------- Descriptors ---------------- #
def compute_descriptors(smiles):
    if not RDKit_AVAILABLE:
        st.error("RDKit not available!")
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    desc = {}
    for name, fn in Descriptors.descList:
        try:
            desc[name] = fn(mol)
        except Exception:
            desc[name] = 0.0
    return desc

# ---------------- Input Processing ---------------- #
def process_input(user_input):
    user_input = user_input.strip()
    smiles = None
    chembl_id = None
    cname = None

    if user_input.lower() in LOCAL_COMPOUNDS:
        smiles = LOCAL_COMPOUNDS[user_input.lower()]
        cname = user_input
    elif is_chembl_id(user_input):
        chembl_id = user_input.upper()
        smiles = chembl_id_to_smiles(chembl_id)
    elif is_valid_smiles(user_input):
        smiles = user_input
    else:
        chembl_id, smiles = chembl_name_to_smiles(user_input)
        cname = user_input

    if not smiles:
        st.error("❌ Could not resolve compound to valid SMILES.")
        return None, None, None

    desc = compute_descriptors(smiles)
    if not desc:
        st.error("❌ RDKit failed to compute descriptors.")
        return None, None, None

    desc["smiles"] = smiles
    desc["chembl_id"] = chembl_id
    return desc, chembl_id, cname

# ---------------- Experimental Data ---------------- #
from chembl_webresource_client.new_client import new_client

def get_top_ic50(chembl_id, top_n=3):
    if not chembl_id:
        return []
    try:
        acts = new_client.activity.filter(molecule_chembl_id=chembl_id, standard_type="IC50")
    except Exception:
        return []
    rows = []
    for a in acts:
        try:
            val = float(a.get("standard_value"))
            log_val = math.log10(val) if val > 0 else None
        except Exception:
            continue
        rows.append({
            "IC50": val,
            "Units": a.get("standard_units"),
            "pChEMBL": a.get("pchembl_value"),
            "log(IC50)": log_val,
            "Target": a.get("target_chembl_id")
        })
    rows.sort(key=lambda x: x["pChEMBL"] or 0, reverse=True)
    return rows[:top_n]

# ---------------- Prediction ---------------- #
def predict_ic50(desc):
    df = pd.DataFrame([desc])
    drop_cols = ["chembl_id", "smiles", "NumRadicalElectrons", "fr_azide", "fr_diazo",
                 "fr_nitroso", "fr_quatN"]
    df.drop(columns=[c for c in drop_cols if c in df.columns], inplace=True)
    saved = joblib.load(MODEL_PATH)
    model, scaler = saved["model"], saved["scaler"]

    # Only keep columns that were in training
    training_cols = saved.get("feature_names")
    if training_cols:
        df = df.reindex(columns=training_cols, fill_value=0.0)

    X = scaler.transform(df)
    log_pred = model.predict(X)[0]

    r2 = saved.get("r2")
    if r2 is None and "X_train" in saved:
        Xtr = scaler.transform(saved["X_train"])
        r2 = r2_score(saved["y_train"], model.predict(Xtr))

    return log_pred, r2

# ---------------- Streamlit UI ---------------- #
st.set_page_config(page_title="PPIM-IC50Pred", layout="wide")
st.title("⚗️ PPIM-IC50Pred")

if not RDKit_AVAILABLE:
    st.error("RDKit is not available. This app cannot run.")
    st.stop()

user_input = st.text_input("Enter compound name, ChEMBL ID, or SMILES:")

if user_input:
    desc, chembl_id, cname = process_input(user_input)
    if desc:
        st.subheader("Compound Details")
        st.write("Name:", cname or "Unknown")
        st.write("ChEMBL ID:", chembl_id or "N/A")
        st.write("SMILES:", desc["smiles"])

        mol = Chem.MolFromSmiles(desc["smiles"])
        if mol:
            st.subheader("2D Structure")
            img = Draw.MolToImage(mol, size=(300, 300))
            st.image(img, use_column_width=False)

        log_ic50, r2 = predict_ic50(desc)
        st.subheader("Prediction")
        st.success(f"Predicted log(IC50): {log_ic50:.4f}")
        if r2:
            st.caption(f"Model R²: {r2:.3f}")

        exp = get_top_ic50(chembl_id)
        if exp:
            st.subheader("Experimental IC50 (ChEMBL)")
            st.dataframe(pd.DataFrame(exp))
        else:
            st.info("No experimental IC50 data found.")
