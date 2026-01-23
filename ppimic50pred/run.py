import streamlit as st
import math
import requests
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from chembl_webresource_client.new_client import new_client
from sklearn.metrics import r2_score
import joblib
import warnings
from sklearn.exceptions import InconsistentVersionWarning
import os
import time

warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
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

    if user_input.lower() in LOCAL_COMPOUNDS:
        smiles = LOCAL_COMPOUNDS[user_input.lower()]
        chembl_id = None
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

def predict_ic50(descriptor_dict, model_path):
    df = pd.DataFrame([descriptor_dict])
    cols_to_drop = [
        'NumRadicalElectrons', 'SMR_VSA8', 'SlogP_VSA9', 'fr_aldehyde', 'fr_azide', 'fr_barbitur',
        'fr_benzodiazepine', 'fr_diazo', 'fr_epoxide', 'fr_isocyan', 'fr_lactam', 'fr_nitroso',
        'fr_prisulfonamd', 'fr_quatN', 'fr_thiocyan', 'fr_term_acetylene', 'fr_phos_ester',
        'fr_oxime', 'fr_dihydropyridine', 'fr_phos_acid', 'fr_hdrzine', 'fr_N_O',
        'chembl_id', 'smiles'
    ]
    df.drop(columns=[c for c in cols_to_drop if c in df.columns], inplace=True)
    
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
        st.markdown(f"**Compound Name:** {compound_name if compound_name else 'Unknown'}")
        st.markdown(f"**ChEMBL ID:** {chembl_id if chembl_id else 'N/A'}")
        st.markdown(f"**SMILES:** {descriptors['smiles']}")

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
