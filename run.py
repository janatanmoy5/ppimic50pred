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

# ----------------- RDKit Imports -----------------
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, Draw
    RDLogger.DisableLog('rdApp.*')
    RDKit_AVAILABLE = True
except ModuleNotFoundError:
    RDKit_AVAILABLE = False

# ---------------- CONFIG ---------------- #
MODEL_PATH = os.path.join(os.path.dirname(__file__), "random_forest_model.pkl")

# ---------------- Local Fallback Compounds ---------------- #
LOCAL_COMPOUNDS = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
}

# ---------------- Utility Functions ---------------- #
def is_chembl_id(s: str) -> bool:
    return s.upper().startswith("CHEMBL")

def is_smiles(s: str) -> bool:
    if not RDKit_AVAILABLE:
        return False
    return Chem.MolFromSmiles(s) is not None

def safe_get(url, retries=3, timeout=10, method="GET", json_data=None, headers=None):
    for attempt in range(retries):
        try:
            if method == "GET":
                r = requests.get(url, timeout=timeout, headers=headers)
            else:
                r = requests.post(url, json=json_data, timeout=timeout, headers=headers)
            r.raise_for_status()
            return r
        except requests.RequestException:
            if attempt < retries - 1:
                time.sleep(1.5)
            else:
                return None

# ---------------- PubChem: name → CID → SMILES ---------------- #
def pubchem_name_to_cid(name: str):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
    r = safe_get(url)
    if not r:
        return None
    cids = r.json().get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None

def pubchem_synonym_to_cid(name: str):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/synonyms/JSON"
    r = safe_get(url)
    if not r:
        return None
    info = r.json().get("InformationList", {}).get("Information", [])
    if not info:
        return None
    return info[0].get("CID")

def pubchem_cid_to_smiles(cid: int):
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid}/property/CanonicalSMILES/JSON"
    )
    r = safe_get(url)
    if not r:
        return None
    props = r.json().get("PropertyTable", {}).get("Properties", [])
    if not props:
        return None
    return props[0].get("CanonicalSMILES")

def resolve_name_to_smiles_pubchem(name: str):
    # 1) direct name → CID
    cid = pubchem_name_to_cid(name)
    if cid:
        smi = pubchem_cid_to_smiles(cid)
        if smi:
            return smi

    # 2) synonym search → CID
    cid = pubchem_synonym_to_cid(name)
    if cid:
        smi = pubchem_cid_to_smiles(cid)
        if smi:
            return smi

    # 3) variants
    variants = {
        name.strip(),
        name.lower(),
        name.upper(),
        name.replace("-", ""),
        name.replace(" ", ""),
        name.replace("‑", ""),
    }
    for v in variants:
        if v == name:
            continue
        cid = pubchem_name_to_cid(v)
        if cid:
            smi = pubchem_cid_to_smiles(cid)
            if smi:
                return smi
        cid = pubchem_synonym_to_cid(v)
        if cid:
            smi = pubchem_cid_to_smiles(cid)
            if smi:
                return smi

    return None

# ---------------- ChEMBL REST helpers ---------------- #
def get_smiles_from_chembl(chembl_id: str):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    r = safe_get(url)
    if not r:
        return None
    return r.json().get("molecule_structures", {}).get("canonical_smiles")

def get_chembl_id_from_smiles_similarity(smiles: str, threshold=0.95):
    url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{smiles}/{threshold}?format=json"
    r = safe_get(url)
    if not r:
        return None
    mols = r.json().get('molecules', [])
    if mols:
        return mols[0].get('molecule_chembl_id')
    return None

def search_chembl_by_name(name: str):
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

def get_top_ic50_values_rest(chembl_id: str, top_n=3):
    if not chembl_id:
        return []
    url = (
        "https://www.ebi.ac.uk/chembl/api/data/activity.json"
        f"?molecule_chembl_id={chembl_id}&standard_type=IC50"
    )
    r = safe_get(url)
    if not r:
        return []
    activities = r.json().get("activities", [])
    valid = []
    for entry in activities:
        pchembl = entry.get('pchembl_value')
        val = entry.get('standard_value')
        units = entry.get('standard_units')
        tid = entry.get('target_chembl_id')
        if pchembl is None or val is None or units is None:
            continue
        try:
            val_num = float(val)
            log10_val = math.log10(val_num) if val_num > 0 else None
        except Exception:
            log10_val = None

        target_name = "Unknown"
        if tid:
            t_url = f"https://www.ebi.ac.uk/chembl/api/data/target/{tid}.json"
            t_res = safe_get(t_url)
            if t_res:
                target_name = t_res.json().get("pref_name") or "Unknown"

        valid.append({
            'chembl_id': chembl_id,
            'ic50_value': val_num,
            'units': units,
            'pchembl_value': pchembl,
            'log10_ic50': log10_val,
            'target_name': target_name,
            'target_id': tid
        })
    valid.sort(key=lambda x: x['pchembl_value'], reverse=True)
    return valid[:top_n]

# ---------------- RDKit descriptors ---------------- #
def compute_rdkit_descriptors(smiles: str):
    if not RDKit_AVAILABLE:
        return {}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {name: func(mol) for name, func in Descriptors.descList}

# ---------------- Input processing (SMILES-first) ---------------- #
def process_input(user_input: str):
    user_input_raw = user_input.strip()
    user_input = user_input_raw.replace(" ", "")

    chembl_id, smiles, compound_name = None, None, None

    # 1) Local fallback
    if user_input_raw.lower() in LOCAL_COMPOUNDS:
        smiles = LOCAL_COMPOUNDS[user_input_raw.lower()]
        compound_name = user_input_raw

    # 2) ChEMBL ID
    elif is_chembl_id(user_input_raw):
        chembl_id = user_input_raw.upper()
        smiles = get_smiles_from_chembl(chembl_id)

    # 3) SMILES directly
    elif is_smiles(user_input_raw):
        smiles = user_input_raw
        chembl_id = get_chembl_id_from_smiles_similarity(smiles)

    # 4) Name → PubChem → SMILES, then optional ChEMBL
    else:
        smiles = resolve_name_to_smiles_pubchem(user_input_raw)
        compound_name = user_input_raw
        if smiles:
            chembl_id = get_chembl_id_from_smiles_similarity(smiles)
        else:
            # fallback: ChEMBL name search
            chembl_id, smiles = search_chembl_by_name(user_input_raw)

    if not smiles:
        st.error(f"Could not resolve '{user_input_raw}' to a valid SMILES using PubChem/ChEMBL.")
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

# ---------------- Potency comparison ---------------- #
def compare_potency(pred_nM, exp_nM, units: str):
    if exp_nM <= 0:
        return "Experimental IC50 value is non-positive; comparison not meaningful."

    ratio = pred_nM / exp_nM

    # Within 2-fold → moderate / good agreement
    if 0.5 <= ratio <= 2.0:
        return (
            f"Model prediction is in good agreement with experimental IC50 "
            f"({pred_nM:.2f} {units} vs {exp_nM:.2f} {units}); "
            f"moderate difference."
        )

    # Predicted IC50 much higher → model underestimates potency (weaker)
    if ratio > 2.0:
        return (
            f"Model underestimates potency (predicts weaker activity): "
            f"predicted {pred_nM:.2f} {units} vs experimental {exp_nM:.2f} {units}."
        )

    # Predicted IC50 much lower → model overestimates potency (stronger)
    if ratio < 0.5:
        return (
            f"Model overestimates potency (predicts stronger activity): "
            f"predicted {pred_nM:.2f} {units} vs experimental {exp_nM:.2f} {units}."
        )

    return (
        f"Model and experimental IC50 differ (predicted {pred_nM:.2f} {units} "
        f"vs experimental {exp_nM:.2f} {units})."
    )

# ---------------- Streamlit Interface ---------------- #
st.set_page_config(page_title="PPIM-IC50Pred", layout="wide")
st.title("⚗️ PPIM-IC50Pred Webserver (SMILES-first)")

user_input = st.text_input("Enter chemical name, ChEMBL ID, or SMILES:")

if user_input:
    descriptors, chembl_id, compound_name = process_input(user_input)

    if descriptors:
        st.subheader("Compound Details")
        st.markdown(f"**Compound Name:** {compound_name if compound_name else 'Unknown'}")
        st.markdown(f"**ChEMBL ID:** {chembl_id if chembl_id else 'N/A'}")
        st.markdown(f"**SMILES:** `{descriptors['smiles']}`")

        # --- 2D Structure ---
        st.subheader("2D Structure Visualization")
        if RDKit_AVAILABLE:
            mol = Chem.MolFromSmiles(descriptors['smiles'])
            if mol:
                img = Draw.MolToImage(mol, size=(300, 300))
                st.image(img, caption=f"{compound_name if compound_name else chembl_id}")
            else:
                st.warning("Could not render molecular structure.")
        else:
            st.warning("2D structure visualization not available in this environment.")

        # --- Prediction ---
        st.subheader("Prediction Details")
        log_val, conf = predict_ic50(descriptors, MODEL_PATH)
        pred_nM = 10 ** log_val
        st.markdown(f"**Predicted log(IC50) [nM]:** `{log_val:.4f}`")
        st.markdown(f"**Predicted IC50:** `{pred_nM:.2f} nM`")
        if conf is not None:
            st.markdown(f"**Model Confidence (R²):** `{conf:.4f}`")

        # --- Experimental Data via ChEMBL REST ---
        st.subheader("Experimental IC50 Values from ChEMBL")
        if chembl_id:
            exp_entries = get_top_ic50_values_rest(chembl_id, top_n=3)
            if exp_entries:
                # Use best (lowest IC50) for comparison
                best = sorted(exp_entries, key=lambda x: x["ic50_value"])[0]
                best_ic50 = best["ic50_value"]
                best_units = best["units"]

                for i, e in enumerate(exp_entries, 1):
                    st.markdown(f"**#{i} Target:** {e['target_name']} ({e['target_id']})")
                    st.markdown(f"- IC50: `{e['ic50_value']:.4f} {e['units']}`")
                    st.markdown(f"- pChEMBL: `{e['pchembl_value']:.4f}`")
                    if e["log10_ic50"] is not None:
                        st.markdown(f"- log(IC50): `{e['log10_ic50']:.4f}`")

                st.markdown("---")
                st.subheader("Predicted vs Experimental Potency")
                comment = compare_potency(pred_nM, best_ic50, best_units)
                st.markdown(comment)
            else:
                st.warning("No experimental IC50 values found in ChEMBL.")
        else:
            st.info("No ChEMBL ID available → experimental IC50 data not retrieved.")
    else:
        st.error("Could not make prediction.")
