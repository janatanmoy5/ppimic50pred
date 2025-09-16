import sys
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
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
import os




# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')

# ---------------- CONFIG ---------------- #

MODEL_PATH = os.path.join(os.path.dirname(__file__), "random_forest_model.pkl")

# ---------------- Utility Functions ---------------- #

def is_chembl_id(s):
    return s.upper().startswith("CHEMBL")

def is_smiles(s):
    return Chem.MolFromSmiles(s) is not None

def search_chembl_by_name(name):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={name}&format=json"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        molecules = r.json().get("molecules", [])
        if molecules:
            mol = molecules[0]
            cid = mol.get("molecule_chembl_id")
            smi = mol.get("molecule_structures", {}).get("canonical_smiles")
            return cid, smi
    except Exception as e:
        print(f"[ERROR] ChEMBL search failed for '{name}': {e}")
    return None, None

def get_smiles_from_chembl(chembl_id):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}.json"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        return r.json().get("molecule_structures", {}).get("canonical_smiles")
    except Exception as e:
        print(f"[ERROR] SMILES retrieval failed for {chembl_id}: {e}")
        return None

def get_chembl_id_from_smiles_similarity(smiles, threshold=0.95):
    url = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{smiles}/{threshold}?format=json"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        mols = r.json().get('molecules', [])
        if mols:
            return mols[0].get('molecule_chembl_id')
    except Exception as e:
        print(f"[ERROR] Similarity search failed: {e}")
    return None

def print_main_page(title):
    """Prints a big main page heading like a server home page."""
    border = "=" * (len(title) + 12)
    print("\n" + border)
    print(f"||   {title.upper()}   ||")
    print(border + "\n")

def print_section(title, subtitle=None):
    """Prints a section heading with optional subtitle."""
    section_border = "-" * (len(title) + 8)
    print(section_border)
    print(f"--  {title}  --")
    print(section_border)
    
    if subtitle:
        print(f"   {subtitle}")
        print("   " + "-" * len(subtitle))
    print()


def compute_rdkit_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {name: func(mol) for name, func in Descriptors.descList}

# ---------------- Data Retrieval ---------------- #

def process_input(user_input):
    user_input = user_input.strip().replace(" ", "")  # Remove all spaces
    chembl_id, smiles, compound_name = None, None, None

    if is_chembl_id(user_input):
        chembl_id = user_input.upper()
        smiles = get_smiles_from_chembl(chembl_id)
    elif is_smiles(user_input):
        smiles = user_input
        chembl_id = get_chembl_id_from_smiles_similarity(smiles)
    else:
        chembl_id, smiles = search_chembl_by_name(user_input)
        compound_name = user_input

    if not smiles:
        print(f"[ERROR] Could not resolve '{user_input}' to a valid SMILES.")
        return None, None, None

    descriptors = compute_rdkit_descriptors(smiles)
    if descriptors is None:
        print("[ERROR] Descriptor calculation failed.")
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
    except Exception as e:
        print(f"[ERROR] Failed to fetch activities: {e}")
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
    # Sort by pchembl_value descending (higher is better affinity)
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

# ---------------- Main ---------------- #

def main():
    print_main_page("PPIM-IC50Pred Webserver")
    CHEM_LOGO = "⚗️ "

    # Keep asking until user provides a valid non-empty input
    user_input = ""
    while not user_input.strip():
        user_input = input(f"{CHEM_LOGO}Enter chemical name, ChEMBL ID, SMILES, or molecular formula: ").strip()

    descriptors, chembl_id, compound_name = process_input(user_input)

    if descriptors:
        print_section("Compound Details")
        print(f"🔎 Compound Name: {compound_name if compound_name else 'Unknown'}")
        print(f"🆔 ChEMBL ID: {chembl_id}")
        print(f"🧪 SMILES: {descriptors['smiles']}")

        print_section("Prediction Details")
        log_val, conf = predict_ic50(descriptors, MODEL_PATH)
        print(f"✅ Predicted Half maximal inhibitory concentration log(IC50): {log_val:.4f} nM")
        if conf is not None:
            print(f"📈 Model Confidence (R²): {conf:.4f} ({conf*100:.2f}%)")

        exp_entries = get_top_ic50_values(chembl_id, top_n=3)
        if exp_entries:
            print_section("Experimental Details")
            print("📊 Top Experimental IC50 values from ChEMBL:")
            for i, e in enumerate(exp_entries, 1):
                exp_log_ic50 = float(e['log10_ic50'])
                diff_val = log_val - exp_log_ic50
                perc_diff = (diff_val / exp_log_ic50) * 100 if exp_log_ic50 != 0 else 0
                print(f"\n#{i} 🎯 Target: {e['target_name']} ({e['target_id']})")
                print(f"   IC50: {float(e['ic50_value']):.4f} {e['units']}")
                print(f"   pChEMBL: {float(e['pchembl_value']):.4f}")
                print(f"   log(IC50): {exp_log_ic50:.4f}")
                print(f"   🔄 Difference (Predicted – Experimental): {diff_val:.4f} nM ({perc_diff:.2f}%)")
        else:
            print("⚠ No experimental IC50 values found.")
    else:
        print("❌ Could not make prediction.")

if __name__ == "__main__":
    main()
