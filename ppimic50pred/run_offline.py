#!/usr/bin/env python3
"""
PPIM-IC50Pred Webserver (CLI version, offline)
Predicts log(IC50) values from chemical SMILES using a Random Forest model.
"""

import os
import sys
import warnings
import pandas as pd
import joblib
from sklearn.metrics import r2_score
from sklearn.exceptions import InconsistentVersionWarning
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

# ---------------- CONFIG ---------------- #

MODEL_PATH = os.path.join(os.path.dirname(__file__), "random_forest_model.pkl")

# Suppress warnings
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
RDLogger.DisableLog("rdApp.*")


# ---------------- Utility Functions ---------------- #

def is_smiles(s: str) -> bool:
    """Check if string is valid SMILES."""
    return Chem.MolFromSmiles(s) is not None


def compute_rdkit_descriptors(smiles: str):
    """Generate all RDKit descriptors for a given SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {name: func(mol) for name, func in Descriptors.descList}


def process_input(smiles: str):
    """Resolve SMILES string into descriptors."""
    if not is_smiles(smiles):
        print(f"[ERROR] '{smiles}' is not a valid SMILES.")
        return None

    descriptors = compute_rdkit_descriptors(smiles)
    if descriptors is None:
        print("[ERROR] Descriptor calculation failed.")
        return None

    descriptors["smiles"] = smiles
    return descriptors


# ---------------- Prediction ---------------- #

def predict_ic50(descriptor_dict: dict, model_path: str):
    """Predict log(IC50) using pretrained model."""
    df = pd.DataFrame([descriptor_dict])

    # drop unused descriptors
    cols_to_drop = [
        "NumRadicalElectrons", "SMR_VSA8", "SlogP_VSA9", "fr_aldehyde", "fr_azide", "fr_barbitur",
        "fr_benzodiazepine", "fr_diazo", "fr_epoxide", "fr_isocyan", "fr_lactam", "fr_nitroso",
        "fr_prisulfonamd", "fr_quatN", "fr_thiocyan", "fr_term_acetylene", "fr_phos_ester",
        "fr_oxime", "fr_dihydropyridine", "fr_phos_acid", "fr_hdrzine", "fr_N_O",
        "smiles"
    ]
    df.drop(columns=[c for c in cols_to_drop if c in df.columns], inplace=True)

    # Load model + scaler
    saved = joblib.load(model_path)
    model, scaler = saved["model"], saved["scaler"]

    # Scale and predict
    X_scaled = scaler.transform(df)
    log_pred = model.predict(X_scaled)[0]

    # Confidence from training R²
    confidence = None
    if "r2" in saved:
        confidence = saved["r2"]
    elif "X_train" in saved and "y_train" in saved:
        X_train_scaled = scaler.transform(saved["X_train"])
        y_train_pred = model.predict(X_train_scaled)
        confidence = r2_score(saved["y_train"], y_train_pred)

    return log_pred, confidence


# ---------------- UI Helpers ---------------- #

def print_main_page(title: str):
    border = "=" * (len(title) + 12)
    print("\n" + border)
    print(f"||   {title.upper()}   ||")
    print(border + "\n")


def print_section(title: str):
    section_border = "-" * (len(title) + 8)
    print(section_border)
    print(f"--  {title}  --")
    print(section_border)
    print()


# ---------------- Main ---------------- #

def main():
    print_main_page("PPIM-IC50Pred Webserver (Offline)")
    CHEM_LOGO = "⚗️ "

    # Keep asking until valid SMILES is provided
    while True:
        user_input = input(f"{CHEM_LOGO}Enter chemical SMILES: ").strip()
        if user_input and is_smiles(user_input):
            break
        print("⚠️ Invalid SMILES. Please try again.\n")

    descriptors = process_input(user_input)
    if not descriptors:
        print("❌ Could not make prediction.")
        return

    # Compound details
    print_section("Compound Details")
    print(f"🧪 SMILES: {descriptors['smiles']}")

    # Prediction
    print_section("Prediction Details")
    log_val, conf = predict_ic50(descriptors, MODEL_PATH)
    print(f"✅ Predicted log(IC50): {log_val:.4f} nM")
    if conf is not None:
        print(f"📈 Model Confidence (R²): {conf:.4f} ({conf*100:.2f}%)")


if __name__ == "__main__":
    main()

