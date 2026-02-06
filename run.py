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
    from rdkit.Chem import AllChem
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

# ---------------- RDKit 3D conformer generation ---------------- #
def generate_3d_molblock(smiles: str):
    """
    Generate a 3D conformer from SMILES using RDKit (ETKDG + UFF),
    and return an SDF/Mol block string for py3Dmol.
    """
    if not RDKit_AVAILABLE:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xf00d
    try:
        AllChem.EmbedMolecule(mol, params)
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        return None
    return Chem.MolToMolBlock(mol)

def show_ligand_3d(smiles: str, width=400, height=350):
    """
    Show a spinning 3D ligand from SMILES using py3Dmol.
    """
    if not PY3DMOL_AVAILABLE:
        st.info("py3Dmol is not installed; 3D ligand visualization is disabled.")
        return
    molblock = generate_3d_molblock(smiles)
    if molblock is None:
        st.info("Could not generate 3D conformer for this ligand.")
        return
    view = py3Dmol.view(width=width, height=height)
    view.addModel(molblock, "sdf")
    view.setStyle({"stick": {"radius": 0.2}, "sphere": {"scale": 0.25}})
    view.zoomTo()
    view.spin(True)
    st.components.v1.html(view._make_html(), height=height + 20, scrolling=False)

# ---------------- Input processing (SMILES-first) ---------------- #
def process_input(user_input: str):
    raw = user_input.strip()
    chembl_ids = set()
    pubchem_meta = None
    smiles = None
    compound_name = raw

    if raw.lower() in LOCAL_COMPOUNDS:
        smiles = LOCAL_COMPOUNDS[raw.lower()]
        return smiles, list(chembl_ids), compound_name, pubchem_meta

    if is_chembl_id(raw):
        cid = raw.upper()
        smi = chembl_get_smiles_from_id(cid)
        if smi:
            smiles = smi
            chembl_ids.add(cid)
            return smiles, list(chembl_ids), compound_name, pubchem_meta

    if is_smiles(raw):
        smiles = raw
        sim_ids = chembl_similarity_search(smiles, threshold=90)
        sub_ids = chembl_substructure_search(smiles)
        chembl_ids.update(sim_ids)
        chembl_ids.update(sub_ids)
        return smiles, list(chembl_ids), compound_name, pubchem_meta

    meta = resolve_name_to_pubchem(raw)
    if meta and meta.get("smiles"):
        smiles = meta["smiles"]
        pubchem_meta = meta
        sim_ids = chembl_similarity_search(smiles, threshold=90)
        sub_ids = chembl_substructure_search(smiles)
        chembl_ids.update(sim_ids)
        chembl_ids.update(sub_ids)
        return smiles, list(chembl_ids), compound_name, pubchem_meta

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

# ---------------- PubChem synonyms for PDB search ---------------- #
def get_pubchem_synonyms(name):
    """Return PubChem depositor-supplied synonyms for a chemical name."""
    if PUBCHEMPY_AVAILABLE:
        try:
            comps = pcp.get_compounds(name, "name")
            if comps:
                syns = comps[0].synonyms
                if syns:
                    return syns
        except Exception:
            pass

    cid = pubchem_name_to_cid_rest(name)
    if not cid:
        return []

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
    r = safe_get(url)
    if not r:
        return []

    try:
        return r.json()["InformationList"]["Information"][0]["Synonym"]
    except Exception:
        return []

# ---------------- RCSB / PDB helpers ---------------- #
def rcsb_search_ligand_safe(chemical_name, max_rows=100):
    """
    Search PDB using PubChem synonyms to avoid JSONDecodeError when
    the original brand name is not present in PDB.
    """
    synonyms = get_pubchem_synonyms(chemical_name)
    search_terms = [chemical_name] + synonyms[:20]

    url = "https://search.rcsb.org/rcsbsearch/v2/query?json="

    for term in search_terms:
        query = {
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {"value": term}
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": max_rows},
                "results_content_type": ["experimental"]
            }
        }

        r = requests.post(url, json=query)
        try:
            data = r.json()
        except Exception:
            continue

        hits = [item["identifier"] for item in data.get("result_set", [])]
        if hits:
            return term, hits

    return None, []

def rcsb_get_assembly_info(pdb_id, assembly_id="1"):
    """Fetch assembly metadata (symmetry, oligomeric state, etc.)."""
    url = f"https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}/{assembly_id}"
    r = safe_get(url)
    if not r:
        return None
    try:
        return r.json()
    except Exception:
        return None

def rcsb_download_url(pdb_id, filetype="pdb"):
    if filetype == "pdb":
        return f"https://files.rcsb.org/download/{pdb_id}.pdb"
    if filetype == "cif":
        return f"https://files.rcsb.org/download/{pdb_id}.cif"
    return None

def show_3d_structure(pdb_id):
    """Display PDB structure using py3Dmol."""
    if not PY3DMOL_AVAILABLE:
        st.info("py3Dmol is not installed; 3D visualization is disabled.")
        return
    pdb_url = rcsb_download_url(pdb_id, "pdb")
    r = safe_get(pdb_url)
    if not r:
        st.info("Could not download PDB file.")
        return
    pdb_data = r.text
    view = py3Dmol.view(width=600, height=500)
    view.addModel(pdb_data, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.zoomTo()
    st.components.v1.html(view._make_html(), height=520, scrolling=False)

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

df_ic50 = None

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

        # 2D structure
        if RDKit_AVAILABLE:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(280, 280))
                st.image(img, caption="2D Structure", use_column_width=False)
            else:
                st.info("Could not render 2D molecular structure.")
        else:
            st.info("RDKit not available; 2D structure rendering disabled.")

        # 3D ligand structure (spinning)
        st.markdown("**3D Ligand Structure (py3Dmol)**")
        show_ligand_3d(smiles, width=380, height=320)

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

    # ---------------- PDB Structures & Assemblies (Chemical Name Search) ---------------- #
    st.markdown("---")
    st.subheader("PDB Structures & Assemblies (Chemical Name Search)")

    st.markdown(
        "Enter any **chemical name** (e.g., Nutlin, Imatinib, Aspirin, Calpol). "
        "The app will resolve PubChem synonyms, search the RCSB PDB, fetch matching "
        "structures, display them in 3D, and extract assembly metadata."
    )

    chem_query = st.text_input(
        "Chemical name for PDB search:",
        value=compound_name if user_input else "",
    )

    if chem_query:
        st.info(f"Searching PDB for ligand (using PubChem synonyms): **{chem_query}** ...")
        matched_term, pdb_hits = rcsb_search_ligand_safe(chem_query)

        if not matched_term or not pdb_hits:
            st.warning("No PDB structures found for this chemical or its synonyms.")
        else:
            st.success(
                f"Found {len(pdb_hits)} structures using synonym: **{matched_term}**"
            )

            st.markdown("### Matching PDB IDs")
            for pid in pdb_hits:
                st.write(
                    f"- **{pid}** — "
                    f"[PDB]({rcsb_download_url(pid,'pdb')}) | "
                    f"[mmCIF]({rcsb_download_url(pid,'cif')})"
                )

            selected_pdb = st.selectbox("Select a PDB ID to visualize:", pdb_hits)

            if selected_pdb:
                st.markdown("### 3D Protein Structure Viewer (PDB)")
                show_3d_structure(selected_pdb)

                st.markdown("### Assembly & Symmetry Information")
                asm = rcsb_get_assembly_info(selected_pdb, assembly_id="1")

                if asm is None:
                    st.warning("No assembly metadata available.")
                else:
                    info = asm.get("rcsb_assembly_info", {})
                    sym = asm.get("rcsb_struct_symmetry", [])

                    st.write(f"**Entry ID:** {info.get('entry_id','N/A').strip()}")
                    st.write(
                        f"**Polymer Composition:** "
                        f"{info.get('polymer_composition','N/A').strip()}"
                    )
                    st.write(
                        f"**Selected Polymer Types:** "
                        f"{info.get('selected_polymer_entity_types','N/A').strip()}"
                    )
                    st.write(
                        f"**Buried Surface Area:** "
                        f"{info.get('total_assembly_buried_surface_area','N/A')}"
                    )

                    st.markdown("### Symmetry / Oligomeric State")
                    if sym:
                        for s in sym:
                            st.markdown(
                                f"- **Kind:** {s.get('kind','N/A').strip()}  \n"
                                f"  **Type:** {s.get('type','N/A').strip()}  \n"
                                f"  **Oligomeric State:** {s.get('oligomeric_state','N/A').strip()}  \n"
                                f"  **Symbol:** {s.get('symbol','N/A').strip()}"
                            )
                    else:
                        st.write("No symmetry annotations available.")
