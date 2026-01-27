import subprocess
import sys
import importlib

# -----------------------------
# Function to install packages
# -----------------------------
def install_package(pkg):
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", pkg])

# -----------------------------
# Required packages
# -----------------------------
required_packages = [
    "streamlit",
    "pandas",
    "rdkit-pypi",
    "chembl-webresource-client",
]

# -----------------------------
# Ensure packages are installed
# -----------------------------
for pkg in required_packages:
    try:
        importlib.import_module(pkg.split("-")[0])
    except ImportError:
        install_package(pkg)

# -----------------------------
# Imports (after ensuring install)
# -----------------------------
import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from chembl_webresource_client.new_client import new_client
from io import BytesIO
from PIL import Image

# -----------------------------
# Streamlit app
# -----------------------------
st.set_page_config(page_title="Chemical Structure Viewer", layout="centered")

st.title("🧪 Chemical Structure Viewer")
st.write("Enter a chemical name to fetch its SMILES and display the structure.")

# User input
chem_name = st.text_input("Chemical Name", placeholder="e.g., Aspirin")

if st.button("Search"):
    if not chem_name.strip():
        st.warning("Please enter a chemical name!")
    else:
        with st.spinner("Fetching chemical information..."):
            try:
                # Search ChEMBL for the molecule
                molecule = new_client.molecule
                res = molecule.filter(pref_name__icontains=chem_name).only(["molecule_chembl_id", "pref_name", "molecule_structures"]).first()
                
                if not res:
                    st.error(f"No molecule found for '{chem_name}'.")
                else:
                    smiles = res.get("molecule_structures", {}).get("canonical_smiles")
                    chembl_id = res.get("molecule_chembl_id")
                    pref_name = res.get("pref_name")
                    
                    if not smiles:
                        st.error(f"SMILES not available for '{pref_name}'.")
                    else:
                        st.success(f"Found molecule: **{pref_name}** (ChEMBL ID: {chembl_id})")
                        st.write(f"**SMILES:** `{smiles}`")

                        # Draw molecule
                        mol = Chem.MolFromSmiles(smiles)
                        img = Draw.MolToImage(mol, size=(400, 400))
                        st.image(img, caption=f"{pref_name} Structure", use_column_width=False)
            except Exception as e:
                st.error(f"An error occurred: {e}")
