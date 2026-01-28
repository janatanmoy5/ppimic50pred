import streamlit as st
import math
import requests
import pandas as pd
import pubchempy as pcp
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Draw
import joblib
import os
import time

RDLogger.DisableLog("rdApp.*")
MODEL_PATH = os.path.join(os.path.dirname(__file__), "random_forest_model.pkl")

# ---------------- NAME NORMALIZATION ---------------- #
NAME_CORRECTIONS = {"nultin": "nutlin", "nutline": "nutlin"}

def normalize_name(name: str):
    name = name.strip().replace("-", "-").replace("–", "-").replace("—", "-")
    return NAME_CORRECTIONS.get(name.lower(), name)

# ---------------- NETWORK UTILITY ---------------- #
def safe_get(url, retries=3, timeout=10):
    for _ in range(retries):
        try:
            r = requests.get(url, timeout=timeout)
            r.raise_for_status()
            return r
        except:
            time.sleep(1)
    return None

# ---------------- PUBCHEM RESOLUTION ---------------- #
def resolve_pubchem(name):
    try:
        c = pcp.get_compounds(name, "name")
        if c:
            return c[0]
    except:
        pass
    try:
        c = pcp.get_compounds(name, "synonym")
        if c:
            return c[0]
    except:
        pass
    return None

# ---------------- ChEMBL LOOKUP ---------------- #
def chembl_name_to_smiles(name):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__icontains={name}"
    r = safe_get(url)
    if r:
        mols = r.json().get("molecules", [])
        if mols and mols[0].get("molecule_structures"):
            return mols[0]["molecule_structures"]["canonical_smiles"]
    return None

def chembl_get_ids_from_smiles(smiles):
    url = "https://www.ebi.ac.uk/chembl/api/data/substructure.json"
    try:
        r = requests.post(url, json={"smiles": smiles}, timeout=15)
        r.raise_for_status()
        mols = r.json().get("molecules", [])
        return [m["molecule_chembl_id"] for m in mols if m.get("molecule_chembl_id")]
    except:
        return []

def chembl_get_ic50_for_molecule(chembl_id, max_records=20):
    url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id={chembl_id}&standard_type=IC50&relation==&limit=1000"
    results = []

    while url and len(results) < max_records:
        r = safe_get(url)
        if not r:
            break

        data = r.json()
        for a in data.get("activities", []):
            val, units = a.get("standard_value"), a.get("standard_units")
            target = a.get("target_chembl_id")
            if not val or not units:
                continue
            try:
                val = float(val)
            except:
                continue

            units = units.lower()
            if units == "nm":
                ic50 = val
            elif units == "um":
                ic50 = val * 1000
            elif units == "mm":
                ic50 = val * 1_000_000
            else:
                continue

            results.append({"Molecule": chembl_id, "Target": target, "IC50_nM": ic50})

        url = data.get("page_meta", {}).get("next")

    results.sort(key=lambda x: x["IC50_nM"])
    return results[:max_records]

# ---------------- RDKit DESCRIPTORS ---------------- #
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    return {name: func(mol) for name, func in Descriptors.descList}

# ---------------- ML PREDICTION ---------------- #
def predict_ic50(descriptor_dict):
    saved = joblib.load(MODEL_PATH)
    model, scaler = saved["model"], saved["scaler"]
    train_features = saved["feature_names"]

    df = pd.DataFrame([descriptor_dict])

    for col in train_features:
        if col not in df.columns:
            df[col] = 0

    df = df[train_features]
    X_scaled = scaler.transform(df)
    log_pred = model.predict(X_scaled)[0]
    return log_pred, saved.get("r2")

# ---------------- STREAMLIT UI ---------------- #
st.set_page_config(page_title="Smart IC50 Predictor", layout="centered")
st.title("🧪 Smart IC50 Prediction Server")

user_input = st.text_input("Enter SMILES, ChEMBL ID, or chemical name:")

if user_input:
    name = normalize_name(user_input)
    smiles = None
    compound = None

    if Chem.MolFromSmiles(name):
        smiles = name
    elif name.upper().startswith("CHEMBL"):
        r = safe_get(f"https://www.ebi.ac.uk/chembl/api/data/molecule/{name.upper()}.json")
        if r:
            smiles = r.json()["molecule_structures"]["canonical_smiles"]
    else:
        compound = resolve_pubchem(name)
        smiles = compound.canonical_smiles if compound else chembl_name_to_smiles(name)

    if not smiles:
        st.error("Could not resolve chemical name to structure.")
        st.stop()

    st.subheader("SMILES")
    st.code(smiles)
    st.image(Draw.MolToImage(Chem.MolFromSmiles(smiles), size=(250, 250)))

    if compound:
        st.subheader("PubChem Info")
        st.write("Formula:", compound.molecular_formula)
        st.write("Molecular Weight:", compound.molecular_weight)
        st.write("XLogP:", compound.xlogp)
        st.write("IUPAC:", compound.iupac_name)

    desc = compute_descriptors(smiles)
    log_ic50, r2 = predict_ic50(desc)
    pred_nM = 10 ** log_ic50

    st.subheader("Predicted IC50")
    st.write(f"log(IC50 nM): {log_ic50:.3f}")
    st.write(f"IC50: {pred_nM:.2f} nM")
    if r2:
        st.write(f"Model R²: {r2:.3f}")

    st.subheader("Experimental IC50 (ChEMBL)")
    chembl_ids = chembl_get_ids_from_smiles(smiles)

    if not chembl_ids:
        st.info("No matching ChEMBL molecules found.")
    else:
        all_data = []
        for cid in chembl_ids[:3]:
            all_data.extend(chembl_get_ic50_for_molecule(cid))

        if not all_data:
            st.info("No experimental IC50 data available.")
        else:
            df_exp = pd.DataFrame(all_data)
            st.dataframe(df_exp)
            best = df_exp.sort_values("IC50_nM").iloc[0]
            st.success(f"Most potent IC50: {best['IC50_nM']:.2f} nM (Target: {best['Target']})")
