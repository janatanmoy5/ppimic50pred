import streamlit as st
import pandas as pd
from chembl_webresource_client.new_client import new_client
import joblib
import os

st.set_page_config(page_title="Chemical Activity Predictor", layout="centered")

st.title("Chemical Activity Predictor")
st.write("Enter a chemical name to fetch details and predict activity.")

# --------------------------
# Load ML model
# --------------------------
MODEL_PATH = "ml_model.pkl"

if os.path.exists(MODEL_PATH):
    model = joblib.load(MODEL_PATH)
    st.success("✅ ML model loaded")
else:
    model = None
    st.warning("⚠️ ML model not found. Predictions will not work.")

# --------------------------
# Input chemical
# --------------------------
chem_name = st.text_input("Chemical Name", "")

if st.button("Search"):
    if chem_name.strip() == "":
        st.error("Please enter a chemical name.")
    else:
        # Fetch molecule from ChEMBL
        molecule = new_client.molecule
        res = molecule.filter(pref_name__iexact=chem_name).only(['molecule_chembl_id', 'pref_name', 'molecule_properties'])
        
        if len(res) == 0:
            st.error(f"No molecule found for '{chem_name}'.")
        else:
            mol = res[0]
            st.success(f"Found: {mol.get('pref_name')}")
            st.write(f"**ChEMBL ID:** {mol.get('molecule_chembl_id')}")
            
            # Show molecular properties
            props = mol.get('molecule_properties')
            if props:
                df_props = pd.DataFrame.from_dict(props, orient='index', columns=['Value'])
                st.subheader("Molecular Properties (ChEMBL)")
                st.dataframe(df_props)
            else:
                st.info("No molecular properties available.")
            
            # --------------------------
            # Predict activity
            # --------------------------
            if model:
                try:
                    # Use available numeric properties for prediction
                    feature_cols = [
                        'full_molweight', 'alogp', 'psa', 'hba', 'hbd'
                    ]
                    X = []
                    for col in feature_cols:
                        val = props.get(col)
                        if val is not None:
                            X.append(float(val))
                        else:
                            X.append(0.0)  # missing values = 0

                    pred = model.predict([X])[0]
                    st.subheader("Predicted Activity")
                    st.write(pred)
                except Exception as e:
                    st.error(f"Prediction failed: {e}")
            else:
                st.info("ML model not loaded, cannot predict.")
