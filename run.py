# run.py
import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from chembl_webresource_client.new_client import new_client
import os
import joblib

# ---------------- CONFIG ---------------- #
MODEL_PATH = "random_forest_model.pkl"

st.set_page_config(page_title="Chemical Activity Predictor", layout="wide")
st.title("💊 Chemical Activity Predictor")

# ---------------- LOAD MODEL ---------------- #
model = None
if os.path.exists(MODEL_PATH):
    try:
        model = joblib.load(MODEL_PATH)
        st.success(f"✅ ML model loaded (type: {type(model)})")
    except Exception as e:
        st.warning(f"⚠ Model could not be loaded: {e}")
else:
    st.info("⚠ Model file not found. Prediction disabled.")

# ---------------- USER INPUT ---------------- #
chem_name = st.text_input("Enter a chemical name to fetch details and predict activity.")
search_button = st.button("Search")

if search_button and chem_name:
    chembl = new_client.molecule
    res = chembl.filter(pref_name__iexact=chem_name).only("molecule_chembl_id", "pref_name", "molecule_structures").first()
    
    if res:
        chembl_id = res.get("molecule_chembl_id", "N/A")
        smiles = res.get("molecule_structures", {}).get("canonical_smiles", "N/A")
        st.subheader(f"Found: {chem_name.upper()}")
        st.write(f"**ChEMBL ID:** {chembl_id}")
        st.write(f"**SMILES:** {smiles}")

        # ---------------- MOLECULAR PROPERTIES ---------------- #
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mw = Descriptors.MolWt(mol)
            alogp = Descriptors.MolLogP(mol)
            psa = Descriptors.TPSA(mol)

            df_props = pd.DataFrame({
                "Property": ["Molecular Weight", "AlogP", "Polar Surface Area"],
                "Value": [mw, alogp, psa]
            })
            st.subheader("Molecular Properties (RDKit)")
            st.table(df_props)

            # ---------------- PREDICTION ---------------- #
            if model and hasattr(model, "predict"):
                try:
                    features = [[mw, alogp, psa]]
                    prediction = model.predict(features)[0]
                    st.subheader("📈 Predicted Activity (IC50)")
                    st.success(f"{prediction:.4f}")
                except Exception as e:
                    st.warning(f"Prediction failed: {e}")
            else:
                if not model:
                    st.info("Prediction unavailable: model missing.")
                else:
                    st.info(f"Prediction unavailable: loaded object is not a model (type: {type(model)})")
        else:
            st.warning("Invalid SMILES. Cannot compute molecular properties.")
    else:
        st.error(f"No results found for '{chem_name}'.")
