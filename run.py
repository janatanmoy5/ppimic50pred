# run.py
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

# ----------------- RDKit Imports (NO DRAWING) -----------------
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
    RDLogger.DisableLog("rdApp.*")
    RDKit_AVAILABLE = True
except Exception:
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
def is_chembl_id(s):
    return s.upper().startswith("CHEMBL")

def is_smiles(s):
    if not RDKit_AVAILABLE:
        return False
    return Chem.MolFromSmiles(s) is not None

def safe_get(url, retries=3, timeout=10):
    for i in range(retries):
        try:
            r = requests.get(url, timeout=timeout)
            r.raise_for_status()
            return r
        except requests.RequestException:
            time.sleep(2)
    return None

def search_chembl_by_name(name):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={name}&format=json"
    r = safe_get(url)
    if not r:
        return None, None
    mols = r.json().get("molecules", [])
    if not mols:
        return None, None
    m = mols[0]
    return (
        m.get("molecule_chembl_id"),
        m.get("molecule_structures", {}).get("canonical_smiles"),
    )

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
    mols = r.json().get("molecules", [])
    return mols[0].get("molecule_chembl_id") if mols else None

def compute_rdkit_descriptors(smiles):
    if not RDKit_AVAILABLE:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {name: func(mol) for name, func in Descriptors.descList}

def process_input(user_input):
    user_input = user_input.strip()
    chembl_id = None
    compound_name = None

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
        return None, None, None

    desc = compute_rdkit_descriptors(smiles)
    if desc is None:
        return None, None, None

    desc["chembl_id"] = chembl_id
    desc["smiles"] = smiles
    return desc, chembl_id, compound_name

# ---------------- ChEMBL Experimental IC50 ---------------- #
from chembl_webresource_client.new_client import new_client

def get_top_ic50_values(chembl_id, top_n=3):
    if not chembl_id:
        return []

    activity = new_client.activity
    target_client = new_client.target

    try:
        res = activity.filter(
            molecule_chembl_id=chembl_id,
            standard_type="IC50"
        )
    except Exception:
        return []

    rows = []
    for r in res:
        try:
            val = float(r.get("standard_value"))
            log_val = math.log10(val) if val > 0 else None
        except Exception:
            continue

        tid = r.get("target_chembl_id")
        target_name = "Unknown"
        if tid:
            try:
                target_name = target_client.get(tid).get("pref_name") or "Unknown"
            except Exception:
                pass

        rows.append({
            "ic50": val,
            "units": r.get("standard_units"),
            "pchembl": r.get("pchembl_value"),
            "log_ic50": log_val,
            "target": target_name,
            "target_id": tid
        })

    rows.sort(key=lambda x: x["pchembl"] or 0, reverse=True)
    return rows[:top_n]

# ---------------- Prediction ---------------- #
def predict_ic50(desc, model_path):
    df = pd.DataFrame([desc])

    drop_cols = [
        "chembl_id", "smiles",
        "NumRadicalElectrons",
        "fr_azide", "fr_diazo", "fr_nitroso",
        "fr_prisulfonamd", "fr_quatN"
    ]
    df.drop(columns=[c for c in drop_cols if c in df.columns], inplace=True)

    saved = joblib.load(model_path)
    model = saved["model"]
    scaler = saved["scaler"]

    X = scaler.transform(df)
    log_pred = model.predict(X)[0]

    conf = saved.get("r2", None)
    if conf is None and "X_train" in saved:
        Xtr = scaler.transform(saved["X_train"])
        conf = r2_score(saved["y_train"], model.predict(Xtr))

    return log_pred, conf

# ---------------- Streamlit UI ---------------- #
st.set_page_config(page_title="PPIM-IC50Pred", layout="wide")
st.title("⚗️ PPIM-IC50Pred")

user_input = st.text_input(
    "Enter compound name, ChEMBL ID, or SMILES",
    placeholder="e.g. aspirin | CHEMBL25 | CC(=O)OC1=CC=CC=C1C(=O)O"
)

if user_input:
    desc, chembl_id, cname = process_input(user_input)

    if desc:
        st.subheader("Compound Information")
        st.write("**Name:**", cname or "Unknown")
        st.write("**ChEMBL ID:**", chembl_id or "N/A")
        st.write("**SMILES:**", desc["smiles"])

        st.subheader("Prediction")
        log_ic50, r2 = predict_ic50(desc, MODEL_PATH)
        st.success(f"Predicted log(IC50): **{log_ic50:.4f}**")

        if r2 is not None:
            st.caption(f"Model R²: {r2:.3f}")

        exp = get_top_ic50_values(chembl_id)
        if exp:
            st.subheader("Experimental IC50 (ChEMBL)")
            st.dataframe(pd.DataFrame(exp))
        else:
            st.info("No experimental IC50 data found.")

    else:
        st.error("Could not process this compound.")
