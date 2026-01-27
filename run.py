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
from sklearn.exceptions import InconsistentVersionWarning
import os
import time

warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
RDLogger.DisableLog("rdApp.*")

# ---------------- CONFIG ---------------- #
MODEL_PATH = os.path.join(os.path.dirname(__file__), "random_forest_model.pkl")

# ---------------- Local Fallback Compounds ---------------- #
LOCAL_COMPOUNDS = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
}

# ---------------- Utility Functions ---------------- #

def is_chembl_id(s: str) -> bool:
    return s.upper().startswith("CHEMBL")

def is_smiles(s: str) -> bool:
    return Chem.MolFromSmiles(s) is not None

def safe_get(url, retries=3, timeout=10):
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=timeout)
            r.raise_for_status()
            return r
        except requests.RequestException:
            if attempt < retries - 1:
                time.sleep(1.5)
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

def get_similar_names(query):
    """Return up to 5 similar chemical names from ChEMBL."""
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={query}&format=json"
    r = safe_get(url)
    if not r:
        return []
    molecules = r.json().get("molecules", [])
    suggestions = []
    for m in molecules[:5]:
        name = m.get("pref_name")
        cid = m.get("molecule_chembl_id")
        if name:
            suggestions.append((name, cid))
    return suggestions

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
    if mols:
        return mols[0].get("molecule_chembl_id")
    return None

def compute_rdkit_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {name: func(mol) for name, func in Descriptors.descList}

def process_input(user_input):
    user_input = user_input.strip()
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
        return None, None, None

    descriptors = compute_rdkit_descriptors(smiles)
    if descriptors is None:
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
        pchembl = entry.get("pchembl_value")
        val = entry.get("standard_value")
        units = entry.get("standard_units")
        if pchembl is None or val is None or units is None:
            continue

        try:
            val_num = float(val)
            log10_val = math.log10(val_num) if val_num > 0 else None
        except Exception:
            log10_val = None

        target_name = "Unknown"
        tid = entry.get("target_chembl_id")
        if tid:
            try:
                target_data = target_client.get(tid)
                target_name = target_data.get("pref_name") or "Unknown"
            except Exception:
                pass

        valid.append(
            {
                "chembl_id": chembl_id,
                "ic50_value": val,
                "units": units,
                "pchembl_value": pchembl,
                "log10_ic50": log10_val,
                "target_name": target_name,
                "target_id": tid,
            }
        )

    valid.sort(key=lambda x: x["pchembl_value"], reverse=True)
    return valid[:top_n]

# ---------------- Prediction ---------------- #

def predict_ic50(descriptor_dict, model_path):
    df = pd.DataFrame([descriptor_dict])

    cols_to_drop = [
        "NumRadicalElectrons",
        "SMR_VSA8",
        "SlogP_VSA9",
        "fr_aldehyde",
        "fr_azide",
        "fr_barbitur",
        "fr_benzodiazepine",
        "fr_diazo",
        "fr_epoxide",
        "fr_isocyan",
        "fr_lactam",
        "fr_nitroso",
        "fr_prisulfonamd",
        "fr_quatN",
        "fr_thiocyan",
        "fr_term_acetylene",
        "fr_phos_ester",
        "fr_oxime",
        "fr_dihydropyridine",
        "fr_phos_acid",
        "fr_hdrzine",
        "fr_N_O",
        "chembl_id",
        "smiles",
    ]

    df.drop(columns=[c for c in cols_to_drop if c in df.columns], inplace=True)

    saved = joblib.load(model_path)
    model, scaler = saved["model"], saved["scaler"]

    X_scaled = scaler.transform(df)
    log_pred = model.predict(X_scaled)[0]

    confidence = None
    if "r2" in saved:
        confidence = saved["r2"]
    else:
        try:
            X_train_scaled = scaler.transform(saved["X_train"])
            y_train_pred = model.predict(X_train_scaled)
            confidence = r2_score(saved["y_train"], y_train_pred)
        except Exception:
            confidence = None

    return log_pred, confidence

# ---------------- Streamlit Layout (DyeMind‑style, Responsive) ---------------- #

st.set_page_config(page_title="PPIM‑IC50Pred", layout="wide")

st.markdown(
    """
    <style>
    body {
        background-color: #f5f6fa;
    }
    .main-container {
        max-width: 900px;
        margin: 0 auto;
        padding: 1.5rem 1rem 3rem 1rem;
    }
    .app-title {
        text-align: center;
        margin-bottom: 0.25rem;
        font-size: 1.9rem;
        font-weight: 600;
    }
    .app-subtitle {
        text-align: center;
        color: #555;
        margin-bottom: 1.5rem;
        font-size: 0.95rem;
    }
    .card {
        background-color: #ffffff;
        border-radius: 12px;
        padding: 1.2rem 1.3rem;
        margin-bottom: 1rem;
        box-shadow: 0 2px 6px rgba(15, 23, 42, 0.06);
    }
    .card-title {
        font-weight: 600;
        margin-bottom: 0.6rem;
        font-size: 1.05rem;
    }
    .small-label {
        font-size: 0.9rem;
        color: #555;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

with st.container():
    st.markdown("<div class='main-container'>", unsafe_allow_html=True)

    st.markdown("<div class='app-title'>PPIM‑IC50Pred</div>", unsafe_allow_html=True)
    st.markdown(
        "<div class='app-subtitle'>IC50 prediction server for small molecules</div>",
        unsafe_allow_html=True,
    )

    # -------- Input Card -------- #
    st.markdown("<div class='card'>", unsafe_allow_html=True)
    st.markdown("<div class='card-title'>Search compound</div>", unsafe_allow_html=True)
    user_input = st.text_input(
        "Enter chemical name, ChEMBL ID, or SMILES:",
        placeholder="e.g., Nutlin-3a, CHEMBL1201733, CC(=O)OC1=CC=CC=C1C(=O)O",
        label_visibility="collapsed",
    )
    st.markdown("</div>", unsafe_allow_html=True)

    if user_input:
        descriptors, chembl_id, compound_name = process_input(user_input)

        # If not found → suggestions + PubMed
        if descriptors is None:
            st.markdown("<div class='card'>", unsafe_allow_html=True)
            st.markdown(
                "<div class='card-title'>Compound not found in ChEMBL</div>",
                unsafe_allow_html=True,
            )
            st.write(
                "The name or identifier you entered could not be resolved to a compound in ChEMBL."
            )

            suggestions = get_similar_names(user_input)
            if suggestions:
                st.markdown("**Did you mean:**")
                for name, cid in suggestions:
                    st.markdown(f"- **{name}** — `{cid}`")
            else:
                st.info("No similar compounds found in ChEMBL.")

            st.markdown("**PubMed reference search:**")
            st.markdown(
                f"[Search PubMed for **{user_input}**](https://pubmed.ncbi.nlm.nih.gov/?term={user_input})"
            )
            st.markdown("</div>", unsafe_allow_html=True)
            st.markdown("</div>", unsafe_allow_html=True)
            st.stop()

        # -------- Compound Summary Card -------- #
        st.markdown("<div class='card'>", unsafe_allow_html=True)
        st.markdown(
            "<div class='card-title'>Compound summary</div>", unsafe_allow_html=True
        )
        st.markdown(
            f"**Name:** {compound_name if compound_name else 'Unknown'}",
        )
        st.markdown(f"**ChEMBL ID:** {chembl_id if chembl_id else 'N/A'}")
        st.markdown(f"**SMILES:** `{descriptors['smiles']}`")
        st.markdown("</div>", unsafe_allow_html=True)

        # -------- Structure Card -------- #
        st.markdown("<div class='card'>", unsafe_allow_html=True)
        st.markdown(
            "<div class='card-title'>2D structure</div>", unsafe_allow_html=True
        )
        mol = Chem.MolFromSmiles(descriptors["smiles"])
        if mol:
            img = Draw.MolToImage(mol, size=(260, 260))
            st.image(img, use_column_width=True)
        else:
            st.info("Could not render molecular structure.")
        st.markdown("</div>", unsafe_allow_html=True)

        # -------- Prediction Card -------- #
        st.markdown("<div class='card'>", unsafe_allow_html=True)
        st.markdown(
            "<div class='card-title'>Predicted IC50</div>", unsafe_allow_html=True
        )
        log_val, conf = predict_ic50(descriptors, MODEL_PATH)
        st.markdown(f"**Predicted log(IC50) [nM]:** `{log_val:.4f}`")
        if conf is not None:
            st.markdown(f"**Model confidence (R²):** `{conf:.4f}`")
        st.markdown("</div>", unsafe_allow_html=True)

        # -------- Experimental Data Card -------- #
        st.markdown("<div class='card'>", unsafe_allow_html=True)
        st.markdown(
            "<div class='card-title'>Experimental IC50 data (ChEMBL)</div>",
            unsafe_allow_html=True,
        )
        exp_entries = get_top_ic50_values(chembl_id, top_n=3)

        if exp_entries:
            for i, e in enumerate(exp_entries, 1):
                st.markdown(f"**#{i} Target:** {e['target_name']} ({e['target_id']})")
                st.markdown(
                    f"- IC50: `{float(e['ic50_value']):.4f} {e['units']}`",
                )
                st.markdown(
                    f"- pChEMBL: `{float(e['pchembl_value']):.4f}`",
                )
                if e["log10_ic50"] is not None:
                    st.markdown(
                        f"- log(IC50): `{float(e['log10_ic50']):.4f}`",
                    )
                st.markdown("---")
        else:
            st.info("No experimental IC50 values found in ChEMBL for this compound.")
        st.markdown("</div>", unsafe_allow_html=True)

    st.markdown("</div>", unsafe_allow_html=True)
