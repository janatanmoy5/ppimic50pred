import streamlit as st
import math
import requests
import pandas as pd
import joblib
import warnings
import os
import time
from sklearn.exceptions import InconsistentVersionWarning

warnings.filterwarnings("ignore", category=InconsistentVersionWarning)

# ----------------- RDKit Imports -----------------
try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, Draw
    RDLogger.DisableLog("rdApp.*")
    RDKit_AVAILABLE = True
except ModuleNotFoundError:
    RDKit_AVAILABLE = False

# ----------------- PubChemPy (optional) -----------------
try:
    import pubchempy as pcp
    PUBCHEMPY_AVAILABLE = True
except ModuleNotFoundError:
    PUBCHEMPY_AVAILABLE = False

# ----------------- py3Dmol (optional, for 3D structures) -----------------
try:
    import py3Dmol
    PY3DMOL_AVAILABLE = True
except ModuleNotFoundError:
    PY3DMOL_AVAILABLE = False

# ---------------- CONFIG ---------------- #
MODEL_PATH = os.path.join(os.path.dirname(__file__), "random_forest_model.pkl")

# ---------------- Local Fallback Compounds ---------------- #
LOCAL_COMPOUNDS = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
}

# ---------------- Utility ---------------- #
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
                time.sleep(0.8)
            else:
                return None

# ---------------- PubChem (REST + PubChemPy) ---------------- #
def pubchem_name_to_cid_rest(name: str):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
    r = safe_get(url)
    if not r:
        return None
    cids = r.json().get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None

def pubchem_synonym_to_cid_rest(name: str):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/synonyms/JSON"
    r = safe_get(url)
    if not r:
        return None
    info = r.json().get("InformationList", {}).get("Information", [])
    if not info:
        return None
    return info[0].get("CID")

def pubchem_cid_to_props_rest(cid: int):
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid}/property/CanonicalSMILES,MolecularFormula,InChIKey/JSON"
    )
    r = safe_get(url)
    if not r:
        return None
    props = r.json().get("PropertyTable", {}).get("Properties", [])
    return props[0] if props else None

def pubchem_autocomplete_name(name: str):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/{name}/JSON"
    r = safe_get(url)
    if not r:
        return []
    return r.json().get("dictionary_terms", {}).get("compound", [])

def resolve_name_to_pubchem(name: str):
    # Try PubChemPy first if available
    if PUBCHEMPY_AVAILABLE:
        try:
            comps = pcp.get_compounds(name, "name")
            if comps:
                c = comps[0]
                return {
                    "smiles": c.canonical_smiles,
                    "cid": c.cid,
                    "formula": c.molecular_formula,
                    "inchikey": c.inchikey,
                    "iupac": c.iupac_name,
                }
        except Exception:
            pass

    # REST: direct name → CID
    cid = pubchem_name_to_cid_rest(name)
    if cid:
        props = pubchem_cid_to_props_rest(cid)
        if props:
            return {
                "smiles": props.get("CanonicalSMILES"),
                "cid": cid,
                "formula": props.get("MolecularFormula"),
                "inchikey": props.get("InChIKey"),
                "iupac": None,
            }

    # REST: synonym search → CID
    cid = pubchem_synonym_to_cid_rest(name)
    if cid:
        props = pubchem_cid_to_props_rest(cid)
        if props:
            return {
                "smiles": props.get("CanonicalSMILES"),
                "cid": cid,
                "formula": props.get("MolecularFormula"),
                "inchikey": props.get("InChIKey"),
                "iupac": None,
            }

    # Variants
    variants = {
        name.strip(),
        name.lower(),
        name.upper(),
        name.replace("-", ""),
        name.replace(" ", ""),
        name.replace("‑", ""),
    }
    for v in variants:
        cid = pubchem_name_to_cid_rest(v)
        if cid:
            props = pubchem_cid_to_props_rest(cid)
            if props:
                return {
                    "smiles": props.get("CanonicalSMILES"),
                    "cid": cid,
                    "formula": props.get("MolecularFormula"),
                    "inchikey": props.get("InChIKey"),
                    "iupac": None,
                }
        cid = pubchem_synonym_to_cid_rest(v)
        if cid:
            props = pubchem_cid_to_props_rest(cid)
            if props:
                return {
                    "smiles": props.get("CanonicalSMILES"),
                    "cid": cid,
                    "formula": props.get("MolecularFormula"),
                    "inchikey": props.get("InChIKey"),
                    "iupac": None,
                }

    # Autocomplete fallback
    suggestions = pubchem_autocomplete_name(name)
    for s in suggestions:
        cid = pubchem_name_to_cid_rest(s)
        if cid:
            props = pubchem_cid_to_props_rest(cid)
            if props:
                return {
                    "smiles": props.get("CanonicalSMILES"),
                    "cid": cid,
                    "formula": props.get("MolecularFormula"),
                    "inchikey": props.get("InChIKey"),
                    "iupac": None,
                }

    return None

# ---------------- ChEMBL REST ---------------- #
def chembl_get_smiles_from_id(chembl_id: str):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    r = safe_get(url)
    if not r:
        return None
    return r.json().get("molecule_structures", {}).get("canonical_smiles")

def chembl_similarity_search(smiles: str, threshold: int = 90, max_hits: int = 20):
    url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{smiles}/{threshold}?format=json"
    r = safe_get(url)
    if not r:
        return []
    mols = r.json().get("molecules", [])
    ids = [m.get("molecule_chembl_id") for m in mols if m.get("molecule_chembl_id")]
    return ids[:max_hits]

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
    mols = r.json().get("molecules", [])
    ids = [m.get("molecule_chembl_id") for m in mols if m.get("molecule_chembl_id")]
    return ids[:max_hits]

def chembl_search_name(name: str):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={name}&format=json"
    r = safe_get(url)
    if not r:
        return None, None
    mols = r.json().get("molecules", [])
    if not mols:
        return None, None
    mol = mols[0]
    return mol.get("molecule_chembl_id"), mol.get("molecule_structures", {}).get("canonical_smiles")

def chembl_get_ic50_for_molecule(chembl_id: str, top_n: int = 5):
    url = (
        "https://www.ebi.ac.uk/chembl/api/data/activity.json"
        f"?molecule_chembl_id={chembl_id}&standard_type=IC50"
    )
    r = safe_get(url)
    if not r:
        return []
    acts = r.json().get("activities", [])
    out = []
    for a in acts:
        pchembl = a.get("pchembl_value")
        val = a.get("standard_value")
        units = a.get("standard_units")
        if pchembl is None or val is None or units is None:
            continue
        try:
            val_num = float(val)
            log10_val = math.log10(val_num) if val_num > 0 else None
            pchembl_val = float(pchembl)
        except Exception:
            continue

        tid = a.get("target_chembl_id")
        target_name = "Unknown"
        target_type = None
        organism = None
        assay_type = a.get("assay_type")

        if tid:
            t_url = f"https://www.ebi.ac.uk/chembl/api/data/target/{tid}.json"
            t_res = safe_get(t_url)
            if t_res:
                tj = t_res.json()
                target_name = tj.get("pref_name") or "Unknown"
                target_type = tj.get("target_type")
                organism = tj.get("organism")

        out.append(
            {
                "chembl_id": chembl_id,
                "ic50_value": val_num,
                "units": units,
                "pchembl_value": pchembl_val,
                "log10_ic50": log10_val,
                "target_name": target_name,
                "target_id": tid,
                "target_type": target_type,
                "organism": organism,
                "assay_type": assay_type,
            }
        )
    out.sort(key=lambda x: x["pchembl_value"], reverse=True)
    return out[:top_n]

def chembl_get_mechanisms_for_ids(chembl_ids):
    rows = []
    seen = set()
    for mid in chembl_ids:
        url = (
            "https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
            f"?molecule_chembl_id={mid}"
        )
        r = safe_get(url)
        if not r:
            continue
        mechs = r.json().get("mechanisms", [])
        for m in mechs:
            tid = m.get("target_chembl_id")
            mech = m.get("mechanism_of_action")
            action_type = m.get("action_type")
            refs = m.get("mechanism_refs", [])
            key = (mid, tid, mech, action_type)
            if key in seen:
                continue
            seen.add(key)

            target_type = None
            organism = None
            if tid:
                t_url = f"https://www.ebi.ac.uk/chembl/api/data/target/{tid}.json"
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
                    "Target ChEMBL ID": tid or "NA",
                    "Target Name": m.get("target_name") or "NA",
                    "Target Type": target_type or "NA",
                    "Organism": organism or "NA",
                    "Mechanism": mech or "NA",
                    "Action Type": action_type or "NA",
                    "Ref Source": ref_source or "NA",
                }
            )
    return rows

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
    raw = user_input.strip()
    chembl_ids = set()
    pubchem_meta = None
    smiles = None
    compound_name = raw

    # 1) Local fallback
    if raw.lower() in LOCAL_COMPOUNDS:
        smiles = LOCAL_COMPOUNDS[raw.lower()]
        return smiles, list(chembl_ids), compound_name, pubchem_meta

    # 2) ChEMBL ID
    if is_chembl_id(raw):
        cid = raw.upper()
        smi = chembl_get_smiles_from_id(cid)
        if smi:
            smiles = smi
            chembl_ids.add(cid)
            return smiles, list(chembl_ids), compound_name, pubchem_meta

    # 3) SMILES directly
    if is_smiles(raw):
        smiles = raw
        sim_ids = chembl_similarity_search(smiles, threshold=90)
        sub_ids = chembl_substructure_search(smiles)
        chembl_ids.update(sim_ids)
        chembl_ids.update(sub_ids)
        return smiles, list(chembl_ids), compound_name, pubchem_meta

    # 4) Name → PubChem → SMILES
    meta = resolve_name_to_pubchem(raw)
    if meta and meta.get("smiles"):
        smiles = meta["smiles"]
        pubchem_meta = meta
        sim_ids = chembl_similarity_search(smiles, threshold=90)
        sub_ids = chembl_substructure_search(smiles)
        chembl_ids.update(sim_ids)
        chembl_ids.update(sub_ids)
        return smiles, list(chembl_ids), compound_name, pubchem_meta

    # 5) Fallback: ChEMBL name search
    cid, smi = chembl_search_name(raw)
    if smi:
        smiles = smi
        if cid:
            chembl_ids.add(cid)
        return smiles, list(chembl_ids), compound_name, pubchem_meta

    return None, [], compound_name, pubchem_meta

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

# ---------------- Per-row impression (log-scale) ---------------- #
def row_impression(pred_log, exp_log):
    if exp_log is None:
        return "N/A"
    delta_log = pred_log - exp_log
    fold = 10 ** abs(delta_log)
    if abs(delta_log) < 0.3:
        label = "close (≈2‑fold)"
    elif delta_log > 0:
        label = "predicted weaker"
    else:
        label = "predicted stronger"
    return f"{label}: Δlog={delta_log:.4f}, ~{fold:.1f}×"

# ---------------- Log-scale potency interpretation ---------------- #
def interpret_potency_logscale(pred_log, exp_log, pred_nM, exp_nM, units):
    if exp_log is None:
        return "Experimental log(IC50) unavailable."

    delta_log = pred_log - exp_log
    fold = 10 ** abs(delta_log)

    if abs(delta_log) < 0.3:
        comment = "Model prediction is in close agreement (within ~2-fold)."
    elif delta_log > 0:
        comment = (
            f"Model underestimates potency by {delta_log:.2f} log units "
            f"(~{fold:.1f}-fold weaker than experiment)."
        )
    else:
        comment = (
            f"Model overestimates potency by {abs(delta_log):.2f} log units "
            f"(~{fold:.1f}-fold stronger than experiment)."
        )

    detail = f"Predicted {pred_nM:.2f} {units} vs experimental {exp_nM:.2f} {units}."
    return f"{comment} {detail}"

# ---------------- Streamlit UI ---------------- #
st.set_page_config(page_title="PPIM‑IC50Pred", layout="wide")
st.markdown(
    "<h2 style='text-align:center;'>PPIM‑IC50Pred</h2>"
    "<p style='text-align:center; color:#555;'>IC50 prediction for small molecules "
    "(PubChem + ChEMBL, SMILES‑first, multi‑target)</p>",
    unsafe_allow_html=True,
)
st.markdown("---")

user_input = st.text_input(
    "Enter chemical name, ChEMBL ID, or SMILES:",
    placeholder="e.g., Nutlin-3a, CHEMBL211045, CC(=O)OC1=CC=CC=C1C(=O)O",
)

df_ic50 = None  # keep in scope for docking & references
pubchem_meta = None

if user_input:
    smiles, chembl_ids, compound_name, pubchem_meta = process_input(user_input)

    if not smiles:
        st.error("Could not resolve input to a valid SMILES using PubChem/ChEMBL.")
        st.stop()

    descriptors = compute_rdkit_descriptors(smiles)
    if descriptors is None:
        st.error("RDKit could not compute descriptors for this SMILES.")
        st.stop()

    descriptors["smiles"] = smiles
    if chembl_ids:
        descriptors["chembl_id"] = chembl_ids[0]

    col1, col2 = st.columns([1, 1])

    with col1:
        st.subheader("Compound & Structure")
        st.markdown(f"**Input Name:** {compound_name}")
        if pubchem_meta:
            st.markdown(f"**PubChem CID:** {pubchem_meta.get('cid', 'N/A')}")
            st.markdown(f"**Formula:** {pubchem_meta.get('formula', 'N/A')}")
            st.markdown(f"**InChIKey:** {pubchem_meta.get('inchikey', 'N/A')}")
        st.markdown(
            f"**Primary ChEMBL IDs (similar/substructure):** "
            f"{', '.join(chembl_ids) if chembl_ids else 'N/A'}"
        )
        st.markdown(f"**SMILES:** `{smiles}`")

        if RDKit_AVAILABLE:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(280, 280))
                st.image(img, caption=compound_name or "Structure", use_column_width=False)
            else:
                st.info("Could not render molecular structure.")
        else:
            st.info("RDKit not available; structure rendering disabled.")

    with col2:
        st.subheader("Predicted IC50")
        log_val, conf = predict_ic50(descriptors, MODEL_PATH)
        pred_nM = 10 ** log_val
        st.markdown(f"- **Predicted log(IC50) [nM]:** `{log_val:.4f}`")
        st.markdown(f"- **Predicted IC50:** `{pred_nM:.2f} nM`")
        if conf is not None:
            st.markdown(f"- **Model R²:** `{conf:.4f}`")

    st.markdown("---")
    st.subheader("Experimental IC50 (ChEMBL) & Potency Interpretation")

    all_ic50_rows = []
    for cid in chembl_ids:
        all_ic50_rows.extend(chembl_get_ic50_for_molecule(cid, top_n=5))

    if not all_ic50_rows:
        st.info("No experimental IC50 values found in ChEMBL for these structures.")
    else:
        df_ic50 = pd.DataFrame(all_ic50_rows)
        df_ic50_display = df_ic50.copy()
        df_ic50_display["IC50 (nM)"] = df_ic50_display["ic50_value"]
        df_ic50_display["log10(IC50)"] = df_ic50_display["log10_ic50"]

        # Impression column only (predicted vs experimental log(IC50))
        impressions = []
        for _, row in df_ic50_display.iterrows():
            exp_log = row["log10_ic50"]
            impressions.append(row_impression(log_val, exp_log))
        df_ic50_display["Impression"] = impressions

        df_ic50_display = df_ic50_display[
            [
                "chembl_id",
                "target_name",
                "target_id",
                "target_type",
                "organism",
                "IC50 (nM)",
                "units",
                "pchembl_value",
                "log10(IC50)",
                "Impression",
                "assay_type",
            ]
        ]
        st.dataframe(df_ic50_display, use_container_width=True)

        best = df_ic50.sort_values("ic50_value").iloc[0]
        best_ic50 = best["ic50_value"]
        best_units = best["units"]
        best_log = best["log10_ic50"]

        st.markdown("**Best experimental IC50 (most potent):**")
        st.markdown(
            f"- Molecule: `{best['chembl_id']}` | Target: `{best['target_name']}` "
            f"({best['target_id']})"
        )
        st.markdown(f"- IC50: `{best_ic50:.2f} {best_units}`")
        if best_log is not None:
            st.markdown(f"- log10(IC50): `{best_log:.4f}`")

        st.markdown("**Predicted vs Experimental (nM + log scale):**")
        exp_log = best_log if best_log is not None else None
        pred_log = log_val

        if exp_log is not None:
            delta_log = pred_log - exp_log
            st.markdown(f"- Predicted log(IC50): `{pred_log:.4f}`")
            st.markdown(f"- Experimental log(IC50): `{exp_log:.4f}`")
            st.markdown(f"- Δlog(IC50) (pred − exp): `{delta_log:.4f}`")
        else:
            st.markdown("- Experimental log(IC50): `N/A`")

        interpretation = interpret_potency_logscale(
            pred_log, exp_log, pred_nM, best_ic50, best_units
        )
        st.markdown(f"- **Interpretation:** {interpretation}")

    st.markdown("---")
    st.subheader("Target Landscape (protein, PPI, cell line, etc.)")

    mech_rows = chembl_get_mechanisms_for_ids(chembl_ids)
    if mech_rows:
        df_mech = pd.DataFrame(mech_rows)
        st.dataframe(df_mech, use_container_width=True)
        st.markdown(
            "_Target types may include single proteins, protein complexes, "
            "protein–protein interactions, and cell line–related targets depending on ChEMBL annotations._"
        )
    else:
        st.info("No mechanism/target annotations found in ChEMBL for these molecules.")

    # ---------------- Docking & 3D Structure Section ---------------- #
    st.markdown("---")
    st.subheader("Docking Study & 3D Structures")

    st.markdown(
        "This section connects predicted/experimental IC50 data with structural "
        "information that can be used for molecular docking studies."
    )

    # Identify most potent target
    best_target_id = None
    best_target_name = None

    if df_ic50 is not None and not df_ic50.empty:
        best_row = df_ic50.sort_values("ic50_value").iloc[0]
        best_target_id = best_row.get("target_id")
        best_target_name = best_row.get("target_name")

    # TARGET STRUCTURE LINKS
    if best_target_id:
        st.markdown(
            f"**Most potent target (from ChEMBL IC50):** "
            f"{best_target_name} ({best_target_id})"
        )

        chembl_target_url = f"https://www.ebi.ac.uk/chembl/target_report_card/{best_target_id}/"
        st.markdown(f"- **ChEMBL Target Page:** [{best_target_id}]({chembl_target_url})")

        rcsb_query = f"https://www.rcsb.org/search?query={best_target_name}"
        st.markdown(
            f"- **Search 3D structures (RCSB PDB):** "
            f"[Open RCSB search]({rcsb_query})"
        )

        uniprot_query = f"https://www.uniprot.org/uniprotkb?query={best_target_name}"
        st.markdown(f"- **UniProt Search:** [Open UniProt]({uniprot_query})")

        pdbe_query = f"https://www.ebi.ac.uk/pdbe/search/pdb?q={best_target_name}"
        st.markdown(f"- **PDBe Search:** [Open PDBe]({pdbe_query})")
    else:
        st.info("No target with experimental IC50 found.")

    # LIGAND STRUCTURE LINKS (PubChem)
    if pubchem_meta and pubchem_meta.get("cid"):
        cid = pubchem_meta["cid"]
        st.markdown("### Ligand 3D Structures (PubChem)")

        sdf_3d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
        st.markdown(f"- **3D SDF Download:** [SDF 3D]({sdf_3d})")

        pdb_3d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/PDB/?record_type=3d"
        st.markdown(f"- **3D PDB Download:** [PDB 3D]({pdb_3d})")

        viewer = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=3D-Conformer"
        st.markdown(f"- **PubChem 3D Viewer:** [Open Viewer]({viewer})")
    else:
        st.info("No PubChem CID available for ligand 3D structure links.")

    st.markdown(
        """
### How to Perform a Docking Study

1. **Obtain a protein 3D structure (PDB)** for the target.  
2. **Prepare the ligand 3D structure** from the SMILES (PubChem links above).  
3. Use a **docking engine** such as AutoDock Vina, Smina, Glide, or GOLD.  
4. Analyze **docking scores and poses** alongside predicted/experimental IC50 values.
"""
    )

    # Optional PDB viewer
    pdb_id = st.text_input(
        "Optional: Enter a PDB ID to visualize the protein structure (e.g., 4HG7):",
        value="",
    )

    if pdb_id:
        if PY3DMOL_AVAILABLE:
            try:
                pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                r = safe_get(pdb_url)
                if r:
                    pdb_block = r.text
                    view = py3Dmol.view(width=500, height=400)
                    view.addModel(pdb_block, "pdb")
                    view.setStyle({"cartoon": {"color": "spectrum"}})
                    view.zoomTo()
                    st.components.v1.html(view._make_html(), height=420, scrolling=False)
                else:
                    st.info("Could not download PDB file.")
            except Exception:
                st.info("3D visualization failed.")
        else:
            st.info(
                "py3Dmol is not installed; 3D visualization is disabled. "
                "Install py3Dmol to enable in-app structure viewing."
            )

    # ---------------- Chemical References & Notes ---------------- #
    st.markdown("---")
    st.subheader("Chemical References & Notes")

    st.markdown(
        "**Yes, ChEMBL provides extensive data linking ligands (small molecules) to "
        "target proteins (receptors), and it is widely used to obtain the necessary "
        "3D structures and 2D data for molecular docking studies.**"
    )

    if df_ic50 is not None and not df_ic50.empty:
        pubmed_links = []
        for _, row in df_ic50.iterrows():
            tid = row.get("target_id")
            if tid:
                query = f"https://pubmed.ncbi.nlm.nih.gov/?term={tid}"
                pubmed_links.append(f"- [{tid} PubMed Search]({query})")

        if pubmed_links:
            st.markdown("**PubMed References (Target-related):**")
            st.markdown("\n".join(pubmed_links))
        else:
            st.markdown("No PubMed references available for these targets.")
    else:
        st.markdown("No experimental IC50/target data available to generate PubMed references.")
