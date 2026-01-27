import streamlit as st
import pandas as pd
import requests
import joblib
import os
from chembl_webresource_client.new_client import new_client

# ---------------- CONFIG ---------------- #
MODEL_PATH = "random_forest_model.pkl"

st.set_page_config(page_title="IC50 Predictor", layout="centered")
st.title("🧪 Compound IC50 Prediction App")

st.markdown("Enter a **chemical name** to fetch details and predict activity.")

# ---------------- Load Model ---------------- #
model = None
if os.path.exists(MODEL_PATH):
    try:
        model = joblib.load(MODEL_PATH)
        st.success("✅ ML model loaded")
    except Exception as e:
        st.warning(f"Model could not be loaded: {e}")
else:
    st.warning("Model file not found. Prediction disabled.")

# ---------------- Helper: Fetch from ChEMBL ---------------- #
def fetch_from_chembl(name):
    try:
        res = new_client.molecule.search(name)
        if res and len(res) > 0:
            mol = res[0]
            return {
                "name": mol.get("pref_name"),
                "smiles": mol.get("molecule_structures", {}).get("canonical_smiles"),
                "chembl_id": mol.get("molecule_chembl_id"),
                "mw": mol.get("molecule_properties", {}).get("full_mwt"),
                "alogp": mol.get("molecule_properties", {}).get("alogp"),
                "psa": mol.get("molecule_properties", {}).get("psa"),
            }
    except:
        pass
    return None

# ---------------- UI ---------------- #
compound_name = st.text_input("Chemical Name")

if st.button("🔍 Search"):

    if not compound_name.strip():
        st.error("Please enter a compound name.")
    else:
        with st.spinner("Fetching compound data..."):
            data = fetch_from_chembl(compound_name)

        if data is None:
            st.error("Compound not found in ChEMBL.")
        else:
            st.success(f"Found: {data['name']}")

            st.markdown(f"**ChEMBL ID:** {data['chembl_id']}")
            st.markdown(f"**SMILES:** `{data['smiles']}`")

            st.subheader("Molecular Properties (ChEMBL)")
            props_df = pd.DataFrame({
                "Property": ["Molecular Weight", "AlogP", "Polar Surface Area"],
                "Value": [data["mw"], data["alogp"], data["psa"]]
            })
            st.table(props_df)

            # ---------------- Prediction ---------------- #
            if model and data["mw"] and data["alogp"] and data["psa"]:
                try:
                    features = [[data["mw"], data["alogp"], data["psa"]]]
                    prediction = model.predict(features)[0]
                    st.subheader("📈 Predicted IC50")
                    st.success(f"{prediction:.4f}")
                except Exception as e:
                    st.warning(f"Prediction failed: {e}")
            else:
                st.info("Prediction unavailable (missing model or properties).")
