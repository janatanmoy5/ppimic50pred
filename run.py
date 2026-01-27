import streamlit as st
import math
import requests
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Draw
from sklearn.metrics import r2_score
import joblib
import warnings
from sklearn.exceptions import InconsistentVersionWarning
import os
import time

warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
RDLogger.DisableLog("rdApp.*")

# ---------------- CONFIG ---------------- #
MODEL_PATH = os.path.join(os.path.dirname(__file__), "random_forest_model.pkl")

# ---------------- Utility ---------------- #

def is_smiles(s: str) -> bool:
    return Chem.MolFromSmiles(s) is not None

def is_chembl_id(s: str) -> bool:
    return s.upper().startswith("CHEMBL")

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
                time.sleep(0.8)
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

def pubchem_synonym_to_cid(name: str):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/synonyms/JSON"
    r = safe_get(url)
    if not r:
        return None
    info = r.json().get("InformationList", {}).get("Information", [])
    if not info:
        return None
    return info[0].get("CID")

def resolve_name_to_smiles(name: str):
    # Try direct name → CID
    cid = pubchem_name_to_cid(name)
    if cid:
        smi = pubchem_cid_to_smiles(cid)
        if smi:
            return smi

    # Try synonym search → CID
    cid = pubchem_synonym_to_cid(name)
    if cid:
        smi = pubchem_cid_to_smiles(cid)
        if smi:
            return smi

    # Try variants
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

# ---------------- ChEMBL: SMILES → ChEMBL IDs + IC50 ---------------- #

def chembl_get_smiles_from_id(chembl_id: str):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    r = safe_get(url)
    if not r:
        return None
    return r.json().get("molecule_structures", {}).get("canonical_smiles")

def chembl_substructure_search(smiles: str, max_hits: int = 20):
    url = "https://www.ebi.ac.uk/chembl/api/data/substructure.json"
    headers = {
        "X-HTTP-Method-Override": "GET",
        "Content-Type": "application/json",
    }
    payload = {"smiles": smiles}
    r = safe_get(url, method="POST", json_data=payload, headers=headers)
    if not r:
        return []
    molecules = r.json().get("molecules", [])
    ids = [m.get("molecule_chembl_id") for m in molecules if m.get("molecule_chembl_id")]
    return ids[:max_hits]

def chembl_get_ic50_for_molecule(chembl_id: str, top_n: int = 3):
    url = (
        "https://www.ebi.ac.uk/chembl/api/data/activity.json"
        f"?molecule_chembl_id={chembl_id}&standard_type=IC50"
    )
    r = safe_get(url)
    if not r:
        return []
    activities = r.json().get("activities", [])
    out = []
    for a in activities:
        pchembl = a.get("pchembl_value")
        val = a.get("standard_value")
        units = a.get("standard_units")
        tid = a.get("target_chembl_id")
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

        out.append(
            {
                "chembl_id": chembl_id,
                "ic50_value": val_num,
                "units": units,
                "pchembl_value": pchembl,
                "log10_ic50": log10_val,
                "target_name": target_name,
                "target_id": tid,
            }
        )
    out.sort(key=lambda x: x["pchembl_value"], reverse=True)
    return out[:top_n]

def chembl_get_targets_from_ids(chembl_ids):
    rows = []
    seen = set()
    for mid in chembl_ids:
        mech_url = (
            "https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
            f"?molecule_chembl_id={mid}"
        )
        r = safe_get(mech_url)
        if not r:
            continue
        mechs = r.json().get("mechanisms", [])
        for m in mechs:
            t_id = m.get("target_chembl_id")
            mechanism = m.get("mechanism_of_action")
            action_type = m.get("action_type")
            refs = m.get("mechanism_refs", [])

            key = (mid, t_id, mechanism, action_type)
            if key in seen:
                continue
            seen.add(key)

            target_type = None
            organism = None
            if t_id:
                t_url = f"https://www.ebi.ac.uk/chembl/api/data/target/{t_id}.json"
                t_res = safe_get(t_url)
                if t_res:
                    tj = t_res.json()
                    target_type = tj.get("target_type")
                    organism = tj.get("organism")

            ref_source = None
            if isinstance(refs, list) and refs:
                ref_source = refs[0].get("ref_type")

            rows.append(
                {
                    "Molecule ChEMBL ID": mid,
                    "Target ChEMBL ID": t_id or "NA",
                    "Target Type": target_type or "NA",
                    "Organism": organism or "NA",
                    "Mechanism": mechanism or "NA",
                    "Action Type": action_type or "NA",
                    "Ref Source": ref_source or "NA",
                }
            )
    return rows

# ---------------- RDKit descriptors ---------------- #

def compute_rdkit_descriptors(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {name: func(mol) for name, func in Descriptors.descList}

# ---------------- Input processing (SMILES-first) ---------------- #

def process_input(user_input: str):
    user_input = user_input.strip()
    smiles = None
    chembl_ids = []

    # 1) If SMILES directly
    if is_smiles(user_input):
        smiles = user_input

    # 2) If ChEMBL ID
    elif is_chembl_id(user_input):
        smiles = chembl_get_smiles_from_id(user_input.upper())
        if smiles:
            chembl_ids.append(user_input.upper())

    # 3) Otherwise: treat as name → PubChem → SMILES
    else:
        smiles = resolve_name_to_smiles(user_input)

    if not smiles:
        return None, []

    # 4) From SMILES → ChEMBL IDs (substructure)
    ids = chembl_substructure_search(smiles)
    chembl_ids = list(set(chembl_ids + ids))

    return smiles, chembl_ids

# ---------------- Prediction ---------------- #

def predict_ic50(descriptor_dict, model_path):
    df = pd.DataFrame([descriptor_dict])

    cols_to_drop = [
        "NumRadicalElectrons","SMR_VSA8","SlogP_VSA9","fr_aldehyde","fr_azide",
        "fr_barbitur","fr_benzodiazepine","fr_diazo","fr_epoxide","fr_isocyan",
        "fr_lactam","fr_nitroso","fr_prisulfonamd","fr_quatN","fr_thiocyan",
        "fr_term_acetylene","fr_phos_ester","fr_oxime","fr_dihydropyridine",
        "fr_phos_acid","fr_hdrzine","fr_N_O","chembl_id","smiles",
    ]
    df.drop(columns=[c for c in cols_to_drop if c in df.columns], inplace=True)

    saved = joblib.load(model_path)
    model, scaler = saved["model"], saved["scaler"]

    X_scaled = scaler.transform(df)
    log_pred = model.predict(X_scaled)[0]

    confidence = saved.get("r2", None)
    return log_pred, confidence

# ---------------- Potency comparison ---------------- #

def compare_potency(pred_nM, exp_nM, units):
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

# ---------------- Streamlit UI ---------------- #

st.set_page_config(page_title="PPIM‑IC50Pred", layout="centered")

st.markdown(
    "<h2 style='text-align:center;'>PPIM‑IC50Pred</h2>"
    "<p style='text-align:center; color:#555;'>IC50 Prediction Server (SMILES‑first)</p>",
    unsafe_allow_html=True,
)

st.markdown("---")

user_input = st.text_input(
    "Enter SMILES, ChEMBL ID, or chemical name:",
    placeholder="e.g., Nutlin-3a, CHEMBL1201733, CC(=O)OC1=CC=CC=C1C(=O)O",
)

if user_input:
    smiles, chembl_ids = process_input(user_input)

    if not smiles:
        st.error("Could not resolve input to a valid SMILES using PubChem/ChEMBL.")
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
    else:
        st.info("Could not render molecular structure.")

    st.markdown("### Predicted IC50")
    log_val, conf = predict_ic50(descriptors, MODEL_PATH)
    pred_nM = 10 ** log_val
    st.markdown(f"- **Predicted log(IC50) [nM]:** `{log_val:.4f}`")
    st.markdown(f"- **Predicted IC50:** `{pred_nM:.2f} nM`")
    if conf is not None:
        st.markdown(f"- **Model confidence (R²):** `{conf:.4f}`")

    st.markdown("### Experimental IC50 (ChEMBL)")
    if not chembl_ids:
        st.info("No matching ChEMBL molecules found from SMILES; experimental IC50 not available.")
    else:
        all_exp = []
        for cid in chembl_ids:
            all_exp.extend(chembl_get_ic50_for_molecule(cid, top_n=3))

        if not all_exp:
            st.info("No experimental IC50 values found in ChEMBL for these structures.")
        else:
            # Show top few entries and compare to best (lowest IC50)
            df_exp = pd.DataFrame(all_exp)
            df_exp_sorted = df_exp.sort_values("ic50_value")
            best = df_exp_sorted.iloc[0]
            best_ic50 = best["ic50_value"]
            best_units = best["units"]

            st.markdown("**Top experimental IC50 entries:**")
            for i, row in df_exp_sorted.head(5).iterrows():
                st.markdown(
                    f"- **Molecule:** {row['chembl_id']} | "
                    f"**Target:** {row['target_name']} ({row['target_id']}) | "
                    f"IC50: `{row['ic50_value']:.2f} {row['units']}` | "
                    f"pChEMBL: `{row['pchembl_value']}`"
                )

            st.markdown("---")
            st.markdown("### Predicted vs Experimental Potency")
            comment = compare_potency(pred_nM, best_ic50, best_units)
            st.markdown(comment)

    st.markdown("### Target information (ChEMBL)")
    if chembl_ids:
        target_rows = chembl_get_targets_from_ids(chembl_ids)
        if target_rows:
            st.table(pd.DataFrame(target_rows))
        else:
            st.info("No target mechanism information available in ChEMBL for these molecules.")
    else:
        st.info("No ChEMBL IDs available → target information not retrieved.")
