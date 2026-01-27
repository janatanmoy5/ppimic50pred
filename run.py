import os
import time
import math
import joblib
import warnings
import requests
import pandas as pd
import streamlit as st

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Draw
from chembl_webresource_client.new_client import new_client
from sklearn.metrics import r2_score
from sklearn.exceptions import InconsistentVersionWarning

warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
RDLogger.DisableLog('rdApp.*')

# ---------------- CONFIG ---------------- #
MODEL_PATH = os.path.join(os.path.dirname(__file__), "random_forest_model.pkl")

LOCAL_COMPOUNDS = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
}

# ---------------- PAGE CONFIG ---------------- #

st.set_page_config(
    page_title="PPIM‑IC50Pred",
    page_icon="🧬",
    layout="wide"
)

# ---------------- RESPONSIVE PUBCHEM‑STYLE CSS ---------------- #

st.markdown("""
<style>

.stApp {
    background-color: #ffffff;
    font-family: "Segoe UI", sans-serif;
}

/* Header */
.header-box {
    background-color: #1f77b4;
    padding: 2rem;
    border-radius: 6px;
    color: white;
    margin-bottom: 2rem;
    text-align: center;
}

/* Section title */
.section-title {
    font-size: 1.7rem;
    font-weight: 600;
    color: #1f4e79;
    margin-top: 1.5rem;
    margin-bottom: 0.8rem;
}

/* Card */
.card {
    background: #f8f9fb;
    padding: 1.5rem;
    border-radius: 6px;
    border: 1px solid #dce3eb;
    margin-bottom: 1.5rem;
}

/* Responsive image */
img {
    max-width: 100%;
    height: auto;
}

/* Footer */
.footer {
    text-align: center;
    color: #777;
    margin-top: 3rem;
    padding: 1rem;
    font-size: 0.9rem;
}

</style>
""", unsafe_allow_html=True)

# ---------------- HEADER ---------------- #

st.markdown("""
<div class="header-box">
    <h1 style="margin-bottom:0;">PPIM‑IC50Pred</h1>
    <h3 style="margin-top:5px;">IC50 Prediction Server</h3>
</div>
""", unsafe_allow_html=True)

# ---------------- FUNCTIONS ---------------- #

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

# ---------------- MAIN UI ---------------- #

st.markdown('<div class="section-title">Search Compound</div>', unsafe_allow_html=True)

user_input = st.text_input("Enter chemical name, ChEMBL ID, or SMILES:")

if user_input:
    descriptors, chembl_id, compound_name = process_input(user_input)

    if descriptors:

        st.markdown('<div class="section-title">Compound Summary</div>', unsafe_allow_html=True)

        col1, col2 = st.columns([1, 2])

        with col1:
            st.markdown('<div class="card">', unsafe_allow_html=True)
            mol = Chem.MolFromSmiles(descriptors['smiles'])
            if mol:
                img = Draw.MolToImage(mol, size=(250, 250))
                st.image(img, caption="Structure", use_column_width=True)
            st.markdown('</div>', unsafe_allow_html=True)

        with col2:
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.write(f"**Name:** {compound_name if compound_name else 'Unknown'}")
            st.write(f"**ChEMBL ID:** {chembl_id if chembl_id else 'N/A'}")
            st.write(f"**SMILES:** {descriptors['smiles']}")
            st.markdown('</div>', unsafe_allow_html=True)

        st.markdown('<div class="section-title">Prediction</div>', unsafe_allow_html=True)

        st.markdown('<div class="card">', unsafe_allow_html=True)
        log_val, conf = predict_ic50(descriptors, MODEL_PATH)
        st.write(f"**Predicted log(IC50) [nM]:** {log_val:.4f}")
        if conf is not None:
            st.write(f"**Model Confidence (R²):** {conf:.4f}")
        st.markdown('</div>', unsafe_allow_html=True)

        st.markdown('<div class="section-title">Experimental IC50 Data</div>', unsafe_allow_html=True)

        st.markdown('<div class="card">', unsafe_allow_html=True)
        exp_entries = get_top_ic50_values(chembl_id, top_n=3)
        if exp_entries:
            for i, e in enumerate(exp_entries, 1):
                st.write(f"**#{i} Target:** {e['target_name']} ({e['target_id']})")
                st.write(f"- IC50: {float(e['ic50_value']):.4f} {e['units']}")
                st.write(f"- pChEMBL: {float(e['pchembl_value']):.4f}")
                st.write(f"- log(IC50): {float(e['log10_ic50']):.4f}")
                st.write("---")
        else:
            st.info("No experimental IC50 values found.")
        st.markdown('</div>', unsafe_allow_html=True)

# ---------------- FOOTER ---------------- #

st.markdown("""
<div class="footer">
Designed by <b>Tanmoy Jana</b> — 2026
</div>
""", unsafe_allow_html=True)
