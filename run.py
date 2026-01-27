import streamlit as st
import requests
import pandas as pd
import joblib
import os

# ---------------- CONFIG ---------------- #
MODEL_PATH = "model.pkl"   # Make sure your model file is in the repo

st.set_page_config(page_title="PPI IC50 Predictor", layout="centered")
st.title("🧪 PPI IC50 Prediction App")
st.write("Enter a chemical name to fetch molecular features from ChEMBL and predict activity.")

# ---------------- LOAD MODEL ---------------- #
@st.cache_resource
def load_model():
    if os.path.exists(MODEL_PATH):
        return joblib.load(MODEL_PATH)
    return None

model = load_model()

# ---------------- FETCH DATA FROM CHEMBL ---------------- #
def get_chembl_features(chem_name):
    search_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={chem_name}&format=json"
    res = requests.get(search_url)

    if res.status_code != 200:
        return None, "ChEMBL API error"

    data = res.json()
    if not data["molecules"]:
        return None, "No molecule found"

    mol = data["molecules"][0]
    props = mol.get("molecule_properties")

    if not props:
        return None, "No molecular properties available"

    try:
        features = {
            "alogp": float(props.get("alogp", 0)),
            "full_mwt": float(props.get("full_mwt", 0)),
            "hba": float(props.get("hba", 0)),
            "hbd": float(props.get("hbd", 0)),
            "psa": float(props.get("psa", 0)),
            "rtb": float(props.get("rtb", 0)),
            "aromatic_rings": float(props.get("aromatic_rings", 0)),
            "heavy_atoms": float(props.get("heavy_atoms", 0)),
        }
    except:
        return None, "Error parsing molecular properties"

    return pd.DataFrame([features]), None

# ---------------- UI ---------------- #
chem_name = st.text_input("🔎 Chemical Name", "")

if st.button("Predict"):
    if chem_name.strip() == "":
        st.warning("Please enter a chemical name.")
    else:
        with st.spinner("Fetching data from ChEMBL..."):
            df, error = get_chembl_features(chem_name)

        if error:
            st.error(error)
        else:
            st.success("Molecular features retrieved!")

            st.subheader("🧬 Molecular Descriptors")
            st.dataframe(df)

            if model:
                prediction = model.predict(df)[0]
                st.subheader("📈 Predicted pIC50")
                st.success(f"Predicted Value: {prediction:.3f}")
            else:
                st.warning("Model file not found. Upload model.pkl to enable predictions.")

# ---------------- FOOTER ---------------- #
st.markdown("---")
st.caption("Data source: ChEMBL Database | No RDKit required 🚀")
