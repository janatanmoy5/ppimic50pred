import streamlit as st
import math
import pandas as pd
import pubchempy as pcp
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Draw
from sklearn.metrics import r2_score
import joblib
import warnings
import os
import requests
import urllib

# ---------------- Setup ---------------- #
warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")
MODEL_PATH = os.path.join(os.path.dirname(__file__), "random_forest_model.pkl")

# ---------------- Utility ---------------- #

def normalize_name(name: str):
    """Normalize input chemical name (replace fancy dashes, strip)."""
    return name.replace("‑", "-").replace("–", "-").strip()

def is_smiles(s: str) -> bool:
    return Chem.MolFromSmiles(s) is not None

def is_chembl_id(s: str) -> bool:
    return s.upper().startswith("CHEMBL")

# ---------------- Robust PubChem + ChEMBL Resolver ---------------- #

def resolve_name_to_smiles(name: str):
    """
    Resolve chemical name to SMILES using PubChemPy and ChEMBL fallback.
    """
    name = normalize_name(name)

    # 1) PubChem exact name search
    try:
        compounds = pcp.get_compounds(name, 'name')
        if compounds:
            return compounds[0].canonical_smiles
    except:
        pass

    # 2) PubChem synonym search
    try:
        compounds = pcp.get_compounds(name, 'synonym')
        if compounds:
            return compounds[0].canonical_smiles
    except:
        pass

    # 3) PubChem formula / similarity search
    try:
        compounds = pcp.get_compounds(name, 'formula')
        if compounds:
            return compounds[0].canonical_smiles
    except:
        pass

    # 4) ChEMBL pref_name fallback
    try:
        chembl_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name={urllib.parse.quote(name)}"
        r = requests.get(chembl_url, timeout=10)
        if r.status_code == 200:
            mols = r.json().get("molecules", [])
            if mols:
                smi = mols[0].get("molecule_structures", {}).get("canonical_smiles")
                if smi:
                    return smi
    except:
        pass

    return None

# ---------------- ChEMBL Functions ---------------- #

def chembl_get_smiles_from_id(chembl_id: str):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    r = requests.get(url, timeout=10)
    if r.status_code != 200:
        return None
    return r.json().get("molecule_structures", {}).get("canonical_smiles")

def chembl_substructure_search(smiles: str, max_hits: int = 20):
    smiles_enc = urllib.parse.quote(smiles)
    url = f"https://www.ebi.ac.uk/chembl/api/data/substructure?smiles={smiles_enc}&limit={max_hits}"
    r = requests.get(url, headers={"Accept": "application/json"}, timeout=10)
    if r.status_code != 200:
        return []
    molecules = r.json().get("molecules", [])
    return [m.get("molecule_chembl_id") for m in molecules if m.get("molecule_chembl_id")]

def chembl_get_ic50_for_molecule(chembl_id: str, top_n: int = 3):
    url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id={chembl_id}&standard_type=IC50"
    r = requests.get(url, timeout=10)
    if r.status_code != 200:
        return []
    activities = r.json().get("activities", [])
    out = []
    for a in activities:
        val = a.get("standard_value")
        units = a.get("standard_units")
        pchembl = a.get("pchembl_value")
        tid = a.get("target_chembl_id")
        if val is None or units is None:
            continue
        try:
            val_num = float(val)
            log10_val = math.log10(val_num) if val_num > 0 else None
        except Exception:
            log10_val = None
        target_name = "Unknown"
        if tid:
            t_url = f"https://www.ebi.ac.uk/chembl/api/data/target/{tid}.json"
            t_res = requests.get(t_url, timeout=10)
            if t_res.status_code == 200:
                target_name = t_res.json().get("pref_name") or "Unknown"
        out.append({
            "chembl_id": chembl_id,
            "ic50_value": val_num,
            "units": units,
            "pchembl_value": pchembl,
            "log10_ic50": log10_val,
            "target_name": target_name,
            "target_id": tid
        })
    out.sort(key=lambda x: x["ic50_value"])
    return out[:top_n]

def chembl_get_targets_from_ids(chembl_ids):
    rows = []
    seen = set()
    for mid in chembl_ids:
        url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism.json?molecule_chembl_id={mid}"
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            continue
        mechs = r.json().get("mechanisms", [])
        for m in mechs:
            key = (mid, m.get("target_chembl_id"), m.get("mechanism_of_action"), m.get("action_type"))
            if key in seen:
                continue
            seen.add(key)
            t_id = m.get("target_chembl_id")
            target_type = organism = "NA"
            if t_id:
                t_url = f"https://www.ebi.ac.uk/chembl/api/data/target/{t_id}.json"
                t_res = requests.get(t_url, timeout=10)
                if t_res.status_code == 200:
                    tj = t_res.json()
                    target_type = tj.get("target_type", "NA")
                    organism = tj.get("organism", "NA")
            rows.append({
                "Molecule ChEMBL ID": mid,
                "Target ChEMBL ID": t_id or "NA",
                "Target Type": target_type,
                "Organism": organism,
                "Mechanism": m.get("mechanism_of_action") or "NA",
                "Action Type": m.get("action_type") or "NA"
            })
    return rows

# ---------------- RDKit descriptors ---------------- #

def compute_rdkit_descriptors(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {name: func(mol) for name, func in Descriptors.descList}

# ---------------- Input Processing ---------------- #

def process_input(user_input: str):
    user_input = normalize_name(user_input)
    smiles = None
    chembl_ids = []

    if is_smiles(user_input):
        smiles = user_input
    elif is_chembl_id(user_input):
        smiles = chembl_get_smiles_from_id(user_input.upper())
        if smiles:
            chembl_ids.append(user_input.upper())
    else:
        smiles = resolve_name_to_smiles(user_input)

    if not smiles:
        return None, []

    ids = chembl_substructure_search(smiles)
    chembl_ids = list(set(chembl_ids + ids))
    return smiles, chembl_ids

# ---------------- Prediction ---------------- #

def predict_ic50(descriptor_dict, model_path):
    df = pd.DataFrame([descriptor_dict])
    df.drop(columns=[c for c in ["smiles","chembl_id"] if c in df.columns], inplace=True)
    saved = joblib.load(model_path)
    model, scaler = saved["model"], saved["scaler"]
    X_scaled = scaler.transform(df)
    log_pred = model.predict(X_scaled)[0]
    conf = saved.get("r2", None)
    return log_pred, conf

def compare_potency(pred_nM, exp_nM, units):
    if exp_nM <= 0:
        return "Experimental IC50 non-positive; comparison not meaningful."
    ratio = pred_nM / exp_nM
    if 0.5 <= ratio <= 2.0:
        return f"Prediction is in good agreement with experimental IC50 ({pred_nM:.2f} vs {exp_nM:.2f} {units})"
    if ratio > 2.0:
        return f"Prediction underestimates potency: predicted {pred_nM:.2f} {units} vs experimental {exp_nM:.2f} {units}"
    return f"Prediction overestimates potency: predicted {pred_nM:.2f} {units} vs experimental {exp_nM:.2f} {units}"

# ---------------- Streamlit UI ---------------- #

st.set_page_config(page_title="PPIM-IC50Pred", layout="centered")
st.markdown("<h2 style='text-align:center;'>PPIM-IC50Pred</h2><p style='text-align:center;color:#555;'>IC50 Prediction Server</p>", unsafe_allow_html=True)
st.markdown("---")

user_input = st.text_input(
    "Enter SMILES, ChEMBL ID, or chemical name:",
    placeholder="Nutlin, Nutlin-3a, CHEMBL1201733, CC(=O)OC1=CC=CC=C1C(=O)O"
)

if user_input:
    smiles, chembl_ids = process_input(user_input)
    if not smiles:
        st.error("Could not resolve input to SMILES via PubChemPy or ChEMBL.")
        st.stop()

    descriptors = compute_rdkit_descriptors(smiles)
    if descriptors is None:
        st.error("RDKit could not compute descriptors for this SMILES.")
        st.stop()
    descriptors["smiles"] = smiles

    st.markdown("### SMILES")
    st.code(smiles, language="text")

    st.markdown("### 2D Structure")
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(260, 260))
        st.image(img, use_column_width=False)

    st.markdown("### Predicted IC50")
    log_val, conf = predict_ic50(descriptors, MODEL_PATH)
    pred_nM = 10 ** log_val
    st.markdown(f"- **Predicted log(IC50) [nM]:** `{log_val:.4f}`")
    st.markdown(f"- **Predicted IC50:** `{pred_nM:.2f} nM`")
    if conf:
        st.markdown(f"- **Model confidence (R²):** `{conf:.4f}`")

    st.markdown("### Experimental IC50 (ChEMBL)")
    if chembl_ids:
        all_exp = []
        for cid in chembl_ids:
            all_exp.extend(chembl_get_ic50_for_molecule(cid))
        if all_exp:
            df_exp = pd.DataFrame(all_exp).sort_values("ic50_value")
            best = df_exp.iloc[0]
            st.markdown(f"**Top IC50:** {best['ic50_value']:.2f} {best['units']} for target {best['target_name']}")
            st.markdown(compare_potency(pred_nM, best["ic50_value"], best["units"]))
        else:
            st.info("No IC50 values found in ChEMBL for these molecules.")
    else:
        st.info("No ChEMBL IDs found from input.")

    st.markdown("### Target information (ChEMBL)")
    if chembl_ids:
        targets = chembl_get_targets_from_ids(chembl_ids)
        if targets:
            st.table(pd.DataFrame(targets))
        else:
            st.info("No target information available.")
    else:
        st.info("No ChEMBL IDs → target information not retrieved.")
