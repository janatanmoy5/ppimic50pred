import streamlit as st
import subprocess
import sys

# -----------------------------
# Automatic package installer
# -----------------------------
def install(package):
    try:
        __import__(package)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Install required packages if missing
for pkg in ["rdkit", "pandas", "numpy", "scikit_learn", "chembl_webresource_client", "requests"]:
    install(pkg)

# -----------------------------
# Imports (after install)
# -----------------------------
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from chembl_webresource_client.new_client import new_client
from PIL import Image
from io import BytesIO

# -----------------------------
# Streamlit App UI
# -----------------------------
st.set_page_config(page_title="Chemical Info Finder", layout="centered")
st.title("🔬 Chemical Info Finder")

st.write("Enter a chemical name to fetch its SMILES and properties from ChEMBL:")

# Input field
chemical_name = st.text_input("Chemical Name", "")

if st.button("Search") and chemical_name.strip():
    st.info(f"Searching for '{chemical_name}'...")

    # Fetch compound from ChEMBL
    molecule = new_client.molecule
    res = molecule.filter(pref_name__icontains=chemical_name).only(["molecule_chembl_id", "pref_name", "molecule_structures"]).first()

    if not res:
        st.error("❌ Chemical not found in ChEMBL.")
    else:
        st.success("✅ Chemical found!")

        chembl_id = res["molecule_chembl_id"]
        name = res.get("pref_name", "N/A")
        smiles = res.get("molecule_structures", {}).get("canonical_smiles", "")

        st.write(f"**ChEMBL ID:** {chembl_id}")
        st.write(f"**Name:** {name}")
        st.write(f"**SMILES:** `{smiles}`")

        # Show molecule image
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(300, 300))
                st.image(img, caption=name)

        # Compute basic properties
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hba = Descriptors.NumHAcceptors(mol)
            hbd = Descriptors.NumHDonors(mol)

            st.write("### Properties")
            st.write(f"- Molecular Weight: {mw:.2f}")
            st.write(f"- LogP: {logp:.2f}")
            st.write(f"- H-Bond Acceptors: {hba}")
            st.write(f"- H-Bond Donors: {hbd}")
