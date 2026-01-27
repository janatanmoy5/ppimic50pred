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
                time.sleep(1.0)
            else:
                return None

# ---------------- ChEMBL helpers (REST only) ---------------- #

def chembl_search_name_first(name):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={name}&format=json"
    r = safe_get(url)
    if not r:
        return None, None
    molecules = r.json().get("molecules", [])
    if not molecules:
        return None, None
    mol = molecules[0]
    cid = mol.get("molecule_chembl_id")
    smi = mol.get("molecule_structures", {}).get("canonical_smiles")
    return cid, smi

def chembl_get_smiles_from_id(chembl_id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    r = safe_get(url)
    if not r:
        return None
    return r.json().get("molecule_structures", {}).get("canonical_smiles")

def chembl_similarity_smiles_to_id(smiles, threshold=0.95):
    url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{smiles}/{threshold}?format=json"
    r = safe_get(url)
    if not r:
        return None
    mols = r.json().get("molecules", [])
    if mols:
        return mols[0].get("molecule_chembl_id")
    return None

def chembl_get_ic50(chembl_id, top_n=3):
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
                "ic50_value": val,
                "units": units,
                "pchembl_value": pchembl,
                "log10_ic50": log10_val,
                "target_name": target_name,
                "target_id": tid,
            }
        )
    out.sort(key=lambda x: x["pchembl_value"], reverse=True)
    return out[:top_n]

def chembl_get_targets_all_molecules(query, chembl_id):
    molecule_ids = set()
    if chembl_id:
        molecule_ids.add(chembl_id)
    if query and (not is_chembl_id(query)) and (not is_smiles(query)):
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={query}&format=json"
        r = safe_get(url)
        if r:
            for m in r.json().get("molecules", [])[:25]:
                mid = m.get("molecule_chembl_id")
                if mid:
                    molecule_ids.add(mid)
    if not molecule_ids:
        return []

    rows = []
    seen = set()
    for mid in molecule_ids:
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

def chembl_similar_names(query):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={query}&format=json"
    r = safe_get(url)
    if not r:
        return []
    suggestions = []
    for m in r.json().get("molecules", [])[:5]:
        name = m.get("pref_name")
        cid = m.get("molecule_chembl_id")
        if name:
            suggestions.append((name, cid))
    return suggestions

# ---------------- PubChem helpers ---------------- #

def pubchem_name_to_cid(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
    r = safe_get(url)
    if not r:
        return None
    cids = r.json().get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None

def pubchem_cid_to_smiles(cid):
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

def pubchem_smiles_normalize(smiles):
    # optional: just validate via PubChem; here we just return original if valid
    return smiles if Chem.MolFromSmiles(smiles) is not None else None

# ---------------- RDKit descriptors ---------------- #

def compute_rdkit_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {name: func(mol) for name, func in Descriptors.descList}

# ---------------- Input processing (ChEMBL + PubChem fallback) ---------------- #

def process_input(user_input):
    user_input = user_input.strip()
    chembl_id, smiles, compound_name, pubchem_cid, source = None, None, None, None, None

    # Local fallback
    if user_input.lower() in LOCAL_COMPOUNDS:
        smiles = LOCAL_COMPOUNDS[user_input.lower()]
        compound_name = user_input
        source = "local"

    # ChEMBL ID
    elif is_chembl_id(user_input):
        chembl_id = user_input.upper()
        smiles = chembl_get_smiles_from_id(chembl_id)
        if smiles:
            source = "chembl_id"

    # SMILES
    elif is_smiles(user_input):
        smiles = user_input
        chembl_id = chembl_similarity_smiles_to_id(smiles)
        source = "smiles"

    # Name: try ChEMBL first, then PubChem
    else:
        chembl_id, smiles = chembl_search_name_first(user_input)
        compound_name = user_input
        if smiles:
            source = "chembl_name"
        else:
            cid = pubchem_name_to_cid(user_input)
            if cid:
                pubchem_cid = cid
                smiles = pubchem_cid_to_smiles(cid)
                if smiles:
                    source = "pubchem_name"

    if not smiles:
        return None, None, None, None, None

    smiles = pubchem_smiles_normalize(smiles) or smiles
    descriptors = compute_rdkit_descriptors(smiles)
    if descriptors is None:
        return None, None, None, None, None

    descriptors["chembl_id"] = chembl_id
    descriptors["smiles"] = smiles
    return descriptors, chembl_id, compound_name, pubchem_cid, source

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

# ---------------- Streamlit Layout ---------------- #

st.set_page_config(page_title="PPIM‑IC50Pred", layout="wide")

st.markdown("""
<style>
body { background-color: #f5f6fa; }
.pc-section {
    background-color: #ffffff;
    border-radius: 8px;
    padding: 1.1rem 1.2rem;
    margin-bottom: 1rem;
    box-shadow: 0 1px 4px rgba(15, 23, 42, 0.06);
}
.pc-section-title {
    font-weight: 600;
    margin-bottom: 0.6rem;
    font-size: 1.05rem;
}
.pc-structure-img img {
    max-width: 260px;
    width: 100%;
    height: auto;
}
</style>
""", unsafe_allow_html=True)

st.markdown("&nbsp;", unsafe_allow_html=True)
st.markdown("## PPIM‑IC50Pred")
st.markdown("IC50 prediction server for small molecules (ChEMBL + PubChem)")

st.markdown("<div class='pc-section'>", unsafe_allow_html=True)
st.markdown("<div class='pc-section-title'>Search compound</div>", unsafe_allow_html=True)
user_input = st.text_input(
    "Enter chemical name, ChEMBL ID, or SMILES:",
    placeholder="e.g., Nutlin, Nutlin-3a, CHEMBL1201733, CC(=O)OC1=CC=CC=C1C(=O)O",
)
st.markdown("</div>", unsafe_allow_html=True)

if user_input:
    descriptors, chembl_id, compound_name, pubchem_cid, source = process_input(user_input)

    if descriptors is None:
        st.error("No matching compound found in ChEMBL or PubChem.")

        st.markdown("""
        ### Why this happens
        - The name may be a trade/brand name or abbreviation  
        - The compound may contain metals or salts (Na⁺, K⁺, Mg²⁺, Ca²⁺)  
        - Naming may follow a different convention (PubChem vs ChEMBL)  
        - The compound may be stored under a different synonym  
        - The name may be misspelled  

        ### What to try
        - Try the IUPAC name  
        - Try the SMILES string  
        - Try the InChIKey  
        - Remove metals/salts (e.g., “Sodium XYZ” → “XYZ”)  
        - Check synonyms in PubChem  
        """)

        st.markdown("### External search")
        st.markdown(f"- PubChem: https://pubchem.ncbi.nlm.nih.gov/#query={user_input}")
        st.markdown(f"- PubMed: https://pubmed.ncbi.nlm.nih.gov/?term={user_input}")

        st.markdown("### Possible ChEMBL suggestions")
        suggestions = chembl_similar_names(user_input)
        if suggestions:
            for name, cid in suggestions:
                st.markdown(f"- **{name}** — `{cid}`")
        else:
            st.info("No similar names found in ChEMBL.")
        st.stop()

    st.markdown("""
    ### 🧪 Important Notes
    - Predicted IC50 comes from the ML model using RDKit descriptors.  
    - Experimental IC50 values (when shown) are retrieved from ChEMBL.  
    - Many chemicals bind to multiple diverse targets, so IC50 varies across:
      - Different proteins  
      - Different cell lines  
      - Different assay conditions  
    - If ChEMBL data are missing, PubChem is still used to obtain SMILES and run the model.
    """)

    st.markdown("<div class='pc-section'>", unsafe_allow_html=True)
    st.markdown("<div class='pc-section-title'>Compound summary</div>", unsafe_allow_html=True)
    st.markdown(f"**Name:** {compound_name or 'Unknown'}")
    st.markdown(f"**ChEMBL ID:** {chembl_id or 'N/A'}")
    st.markdown(f"**PubChem CID:** {pubchem_cid or 'N/A'}")
    st.markdown(f"**Source used for SMILES:** {source or 'Unknown'}")
    st.markdown(f"**SMILES:** `{descriptors['smiles']}`")
    st.markdown("</div>", unsafe_allow_html=True)

    st.markdown("<div class='pc-section pc-structure-img'>", unsafe_allow_html=True)
    st.markdown("<div class='pc-section-title'>2D structure</div>", unsafe_allow_html=True)
    mol = Chem.MolFromSmiles(descriptors["smiles"])
    if mol:
        img = Draw.MolToImage(mol, size=(260, 260))
        st.image(img)
    else:
        st.info("Could not render molecular structure.")
    st.markdown("</div>", unsafe_allow_html=True)

    st.markdown("<div class='pc-section'>", unsafe_allow_html=True)
    st.markdown("<div class='pc-section-title'>Predicted IC50</div>", unsafe_allow_html=True)
    log_val, conf = predict_ic50(descriptors, MODEL_PATH)
    st.markdown(f"**Predicted log(IC50) [nM]:** `{log_val:.4f}`")
    ic50_nM = 10 ** log_val
    st.markdown(f"**Predicted IC50 (nM):** `{ic50_nM:.2f}`")
    if ic50_nM <= 100:
        st.success("Strong predicted potency (low IC50).")
    elif ic50_nM <= 1000:
        st.info("Moderate predicted potency.")
    else:
        st.warning("Weak predicted potency (high IC50).")
    if conf is not None:
        st.markdown(f"**Model confidence (R²):** `{conf:.4f}`")
    st.markdown("</div>", unsafe_allow_html=True)

    st.markdown("<div class='pc-section'>", unsafe_allow_html=True)
    st.markdown("<div class='pc-section-title'>Experimental IC50 (ChEMBL)</div>", unsafe_allow_html=True)
    if chembl_id:
        exp_entries = chembl_get_ic50(chembl_id, top_n=3)
        if exp_entries:
            for i, e in enumerate(exp_entries, 1):
                st.markdown(f"**#{i} Target:** {e['target_name']} ({e['target_id']})")
                st.markdown(f"- IC50: `{e['ic50_value']} {e['units']}`")
                st.markdown(f"- pChEMBL: `{e['pchembl_value']}`")
                if e["log10_ic50"] is not None:
                    st.markdown(f"- log(IC50): `{e['log10_ic50']}`")
                st.markdown("---")
        else:
            st.info("No experimental IC50 values found in ChEMBL for this compound.")
    else:
        st.info("No ChEMBL ID available → experimental IC50 cannot be retrieved.")
    st.markdown("</div>", unsafe_allow_html=True)

    st.markdown("<div class='pc-section'>", unsafe_allow_html=True)
    st.markdown("<div class='pc-section-title'>Target information (ChEMBL, all matching molecules)</div>", unsafe_allow_html=True)
    if chembl_id:
        target_rows = chembl_get_targets_all_molecules(user_input, chembl_id)
        if target_rows:
            st.table(pd.DataFrame(target_rows))
        else:
            st.info("No target information available in ChEMBL for this query.")
    else:
        st.info("No ChEMBL ID available → target information cannot be retrieved.")
    st.markdown("</div>", unsafe_allow_html=True)

    st.markdown("<div class='pc-section'>", unsafe_allow_html=True)
    st.markdown("<div class='pc-section-title'>External resources</div>", unsafe_allow_html=True)
    query_name = compound_name if compound_name else descriptors["smiles"]
    st.markdown(f"- **PubChem compound search:** https://pubchem.ncbi.nlm.nih.gov/#query={query_name}")
    st.markdown(f"- **PubMed literature search:** https://pubmed.ncbi.nlm.nih.gov/?term={query_name}")
    if chembl_id:
        st.markdown(f"- **ChEMBL compound page:** https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}")
    st.markdown("</div>", unsafe_allow_html=True)
