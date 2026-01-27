import os
import time
import math
import joblib
import warnings
import requests
import numpy as np
import streamlit as st

from chembl_webresource_client.new_client import new_client
from sklearn.metrics import r2_score
from sklearn.exceptions import InconsistentVersionWarning

warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

# ---------------- CONFIG ---------------- #
MODEL_PATH = os.path.join(os.path.dirname(__file__), "random_forest_model.pkl")

LOCAL_COMPOUNDS = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
}

# Remote RDKit descriptor API (you can swap this to your own service)
RDKIT_DESCRIPTOR_API = "https://rdkit-api.johnsnowlabs.com/descriptors"

# ---------------- Utility Functions ---------------- #

def is_chembl_id(s: str) -> bool:
    return s.upper().startswith("CHEMBL")

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

def compute_rdkit_descriptors_remote(smiles: str):
    """
    Call a remote RDKit descriptor service that returns a dict of descriptors.
    Assumes the API returns keys compatible with the model training.
    """
    try:
        resp = requests.post(
            RDKIT_DESCRIPTOR_API,
            json={"smiles": smiles},
            timeout=20,
        )
        resp.raise_for_status()
        data = resp.json()
        if not isinstance(data, dict) or not data:
            st.error("Descriptor API returned an unexpected response.")
            return None
        return data
    except Exception as e:
        st.error(f"Descriptor API error: {e}")
        return None

def process_input(user_input):
    user_input = user_input.strip()
    chembl_id, smiles, compound_name = None, None, None

    # Local fallback
    key = user_input.lower().replace(" ", "")
    if key in LOCAL_COMPOUNDS:
        smiles = LOCAL_COMPOUNDS[key]
        compound_name = user_input

    elif is_chembl_id(user_input):
        chembl_id = user_input.upper()
        smiles = get_smiles_from_chembl(chembl_id)

    else:
        # Try name search first
        chembl_id, smiles = search_chembl_by_name(user_input)
        compound_name = user_input
        # If that fails, treat as SMILES and try similarity
        if not smiles:
            smiles = user_input
            chembl_id = get_chembl_id_from_smiles_similarity(smiles)

    if not smiles:
        st.error(f"Could not resolve '{user_input}' to a valid SMILES.")
        return None, None, None

    descriptors = compute_rdkit_descriptors_remote(smiles)
    if descriptors is None:
        st.error("Descriptor calculation failed.")
        return None, None, None

    descriptors["chembl_id"] = chembl_id
    descriptors["smiles"] = smiles
    return descriptors, chembl_id, compound_name

# ---------------- Experimental Data ---------------- #

def get_top_ic50_values(chembl_id, top_n=3):
    if not chembl_id:
        return []
    activity = new_client.activity
    target_client = new_client.target
    try:
        res = activity.filter(molecule_chembl_id=chembl_id, standard_type="IC50")
    except Exception:
        return []
    valid = []
    for entry in res:
        pchembl = entry.get('pchembl_value')
        val = entry.get('standard_value')
        units = entry.get('standard_units')
        if pchembl is None or val is None or units is None:
            continue
        try:
            val_num = float(val)
            log10_val = math.log10(val_num) if val_num > 0 else None
        except Exception:
            log10_val = None
        target_name = "Unknown"
        tid = entry.get('target_chembl_id')
        if tid:
            try:
                target_data = target_client.get(tid)
                target_name = target_data.get('pref_name') or "Unknown"
            except Exception:
                pass
        valid.append({
            'chembl_id': chembl_id,
            'ic50_value': val,
            'units': units,
            'pchembl_value': pchembl,
            'log10_ic50': log10_val,
            'target_name': target_name,
            'target_id': tid
        })
    valid.sort(key=lambda x: x['pchembl_value'], reverse=True)
    return valid[:top_n]

# ---------------- Prediction ---------------- #

def build_feature_vector(descriptor_dict, saved_obj):
    """
    Build a NumPy feature vector in the exact order expected by the scaler/model.
    We try, in order:
      1) saved_obj["feature_names"] if present
      2) saved_obj["feature_order"] if present
      3) infer from descriptor_dict keys minus dropped columns, sorted
    """
    cols_to_drop = {
        'NumRadicalElectrons', 'SMR_VSA8', 'SlogP_VSA9', 'fr_aldehyde', 'fr_azide', 'fr_barbitur',
        'fr_benzodiazepine', 'fr_diazo', 'fr_epoxide', 'fr_isocyan', 'fr_lactam', 'fr_nitroso',
        'fr_prisulfonamd', 'fr_quatN', 'fr_thiocyan', 'fr_term_acetylene', 'fr_phos_ester',
        'fr_oxime', 'fr_dihydropyridine', 'fr_phos_acid', 'fr_hdrzine', 'fr_N_O',
        'chembl_id', 'smiles'
    }

    feature_names = None
    if isinstance(saved_obj, dict):
        if "feature_names" in saved_obj:
            feature_names = saved_obj["feature_names"]
        elif "feature_order" in saved_obj:
            feature_names = saved_obj["feature_order"]

    if feature_names is None:
        # Infer from descriptor_dict keys
        keys = [k for k in descriptor_dict.keys() if k not in cols_to_drop]
        keys.sort()
        feature_names = keys

    # Build vector
    vec = []
    for name in feature_names:
        val = descriptor_dict.get(name, 0.0)
        try:
            val = float(val)
        except Exception:
            val = 0.0
        vec.append(val)

    X = np.array(vec, dtype=float).reshape(1, -1)
    return X, feature_names

def predict_ic50(descriptor_dict, model_path):
    saved = joblib.load(model_path)
    model = saved["model"]
    scaler = saved["scaler"]

    X, feature_names = build_feature_vector(descriptor_dict, saved)
    X_scaled = scaler.transform(X)
    log_pred = model.predict(X_scaled)[0]

    confidence = None
    if "r2" in saved:
        confidence = saved["r2"]
    elif "X_train" in saved and "y_train" in saved:
        X_train = saved["X_train"]
        y_train = saved["y_train"]
        # If X_train was stored as pandas originally, it might be array-like now
        X_train = np.array(X_train)
        X_train_scaled = scaler.transform(X_train)
        y_train_pred = model.predict(X_train_scaled)
        confidence = r2_score(y_train, y_train_pred)

    return log_pred, confidence

# ---------------- Streamlit Interface ---------------- #

st.set_page_config(page_title="PPIM-IC50Pred", layout="wide")
st.title("⚗️ PPIM-IC50Pred Webserver")

user_input = st.text_input("Enter chemical name, ChEMBL ID, or SMILES:")

if user_input:
    descriptors, chembl_id, compound_name = process_input(user_input)

    if descriptors:
        st.subheader("Compound Details")
        st.markdown(f"**Compound Name:** {compound_name if compound_name else 'Unknown'}")
        st.markdown(f"**ChEMBL ID:** {chembl_id if chembl_id else 'N/A'}")
        st.markdown(f"**SMILES:** {descriptors['smiles']}")

        # --- 2D Structure via external image services ---
        st.subheader("2D Structure Visualization")
        img_url = None
        if chembl_id:
            img_url = f"https://www.ebi.ac.uk/chembl/api/data/image/{chembl_id}"
        else:
            img_url = f"https://cactus.nci.nih.gov/chemical/structure/{descriptors['smiles']}/image"
        if img_url:
            st.image(img_url, caption=f"{compound_name if compound_name else chembl_id}", use_column_width=False)
        else:
            st.warning("Could not render molecular structure.")

        # --- Prediction ---
        st.subheader("Prediction Details")
        log_val, conf = predict_ic50(descriptors, MODEL_PATH)
        st.markdown(f"**Predicted log(IC50) [nM]:** {log_val:.4f}")
        if conf is not None:
            st.markdown(f"**Model Confidence (R²):** {conf:.4f} ({conf*100:.2f}%)")

        # --- Experimental Data ---
        exp_entries = get_top_ic50_values(chembl_id, top_n=3)
        if exp_entries:
            st.subheader("Experimental IC50 Values from ChEMBL")
            for i, e in enumerate(exp_entries, 1):
                exp_log_ic50 = float(e['log10_ic50'])
                diff_val = log_val - exp_log_ic50
                perc_diff = (diff_val / exp_log_ic50) * 100 if exp_log_ic50 != 0 else 0
                st.markdown(f"**#{i} Target:** {e['target_name']} ({e['target_id']})")
                st.markdown(f"- IC50: {float(e['ic50_value']):.4f} {e['units']}")
                st.markdown(f"- pChEMBL: {float(e['pchembl_value']):.4f}")
                st.markdown(f"- log(IC50): {exp_log_ic50:.4f}")
                st.markdown(f"- 🔄 Difference (Predicted – Experimental): {diff_val:.4f} nM ({perc_diff:.2f}%)")
        else:
            if chembl_id:
                st.warning("No experimental IC50 values found in ChEMBL.")
            else:
                st.info("Experimental IC50 data not available for local compounds.")
    else:
        st.error("Could not make prediction.")
