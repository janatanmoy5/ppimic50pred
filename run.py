import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import joblib
import os

# ---------------- Page config ---------------- #
st.set_page_config(
    page_title="PPI MIC50 Predictor",
    layout="wide"
)

st.title("🔬 PPI MIC50 Prediction App")
st.markdown("Predict MIC50 values from SMILES input. Enter SMILES below.")

# ---------------- Sidebar ---------------- #
st.sidebar.header("Input SMILES")
smiles_input = st.sidebar.text_area(
    "Enter one SMILES per line:",
    height=150
)
smiles_list = [s.strip() for s in smiles_input.split("\n") if s.strip()]

# ---------------- Load model ---------------- #
MODEL_PATH = "random_forest_model.pkl"
if os.path.exists(MODEL_PATH):
    model = joblib.load(MODEL_PATH)
else:
    model = None
    st.warning("⚠️ Random Forest model not found. Predictions will not run.")

# ---------------- Descriptor function ---------------- #
def calculate_descriptors(mol):
    """
    Simple RDKit descriptors for demonstration
    """
    from rdkit.Chem import Descriptors
    return {
        "MolWt": Descriptors.MolWt(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumRotatableBonds": Descriptors.NumRotatableBonds(mol)
    }

# ---------------- Prediction function ---------------- #
def predict_activity(descriptors):
    """
    Predict MIC50 using the loaded model
    """
    if model is None:
        return None
    df = pd.DataFrame([descriptors])
    pred = model.predict(df)
    return pred[0]

# ---------------- Main App ---------------- #
if smiles_list:
    for idx, smi in enumerate(smiles_list, 1):
        st.markdown(f"### Compound {idx}")
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            st.error(f"❌ Invalid SMILES: {smi}")
            continue

        st.success(f"✅ Valid SMILES: {smi}")

        # Draw molecule
        st.image(Draw.MolToImage(mol, size=(300, 300)))

        # Calculate descriptors
        desc = calculate_descriptors(mol)
        st.write("**Descriptors:**")
        st.json(desc)

        # Predict MIC50
        pred = predict_activity(desc)
        if pred is not None:
            st.metric("Predicted MIC50", f"{pred:.4f}")
        else:
            st.info("Model not loaded. Cannot predict MIC50.")

else:
    st.info("Enter SMILES in the sidebar to start predictions.")
