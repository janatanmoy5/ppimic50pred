# ppimic50pred

Supervised machine learning framework for predicting small-molecule bioactivity (IC₅₀) against protein–protein interactions (PPIs), single proteins, and cell-based assays.

This repository accompanies the manuscript **"Supervised Learning-Driven Prediction of Small-Molecule Modulator Bioactivity Against Protein–Protein Interactions"** (Jana T, Karmakar J, Banerjee R, Saha S; submitted to *Journal of Computer-Aided Molecular Design*). It provides the trained random forest model, prediction scripts, and analysis code used in that study, benchmarking RDKit, PubChem, and PaDEL molecular descriptors across four supervised regressors (random forest, gradient boosting, support vector regression, and LSTM) on a curated ChEMBL dataset of 3,451 compounds spanning 176 biological targets.

🌐 **Run** [ppimic50pred.onrender.com](https://ppimic50pred.onrender.com/)

## Overview

- **Input:** a small-molecule SMILES string, ChEMBL ID, compound name, or molecular formula.
- **Output:** predicted Log₁₀(IC₅₀), with an optional comparison against experimentally reported ChEMBL bioactivity values when run in online mode.
- **Model:** random forest regressor trained on 217 RDKit descriptors (195 after feature selection), achieving R² = 0.75 / RMSE = 0.78 in cross-validation and R² = 0.74 / RMSE = 0.79 on a held-out blind set. See the manuscript for full performance benchmarks, feature-importance analysis, and known generalization limits on structurally novel targets.

## Repository structure

This repository provides **two separate interfaces** to the same underlying random forest model:

| Path | Description |
|---|---|
| `run.py` (repo root) | **Streamlit web application** powering the live demo; launched via `streamlit run run.py` |
| `random_forest_model.pkl` (repo root) | Trained random forest model + fitted scaler used by the Streamlit app |
| `ppimic50pred/` | Nested Python package providing a **command-line interface**: `run_offline.py`, `run_online.py`, its own `run.py`, `random_forest_model.pkl`, `requirements.txt`, `setup.py`, `__init__.py` |
| `analysis/analysis/PPIM-IC50pred.ipynb` | Jupyter notebook reproducing the study's descriptor generation, feature selection/SHAP analysis, model training and cross-validation, clustering, and multitarget validation — see [Notebook](#notebook) below |
| `analysis/dataset/` | Processed descriptor/bioactivity CSVs: `rdkit_train.csv`, `rdkit_blind.csv`, `pubchem_train.csv`, `pubchem_blind.csv`, `padel_blind.csv`, `padel_train.csv` ,|
| `analysis/dataset/rdkit/` | Downstream RDKit-based analysis outputs: feature importances, chemical-class breakdowns, structure- and feature-based cluster assignments, Tanimoto similarity heatmap — see [Data](#data) below |
| `Dockerfile` | Container definition; launches the Streamlit app via `streamlit run run.py` |
| `render.yaml` | Deployment configuration for the hosted demo (Render) |
| `requirements.txt` (repo root) | Dependencies for the Streamlit app |


## Installation

```bash
git clone https://github.com/janatanmoy5/ppimic50pred.git
cd ppimic50pred
```

## Usage

### Option A: command-line interface (offline / online)

The CLI scripts live in the nested `ppimic50pred/` package folder — install its dependencies from there, not the repo-root `requirements.txt`:

```bash
cd ppimic50pred
pip install -r requirements.txt
```

**Offline mode** — predicts Log₁₀(IC₅₀) from a SMILES string using only the local model; no network access required:

```bash
python run_offline.py
```

```
=============================================
||   PPIM-IC50PRED WEBSERVER (OFFLINE)    ||
=============================================

⚗️ Enter chemical SMILES: CCOC(=O)c1nc(N(C)Cc2ccccc2)c2c(C)noc2n1
------------------------
--  Compound Details  --
------------------------
🧪 SMILES: CCOC(=O)c1nc(N(C)Cc2ccccc2)c2c(C)noc2n1

--------------------------
--  Prediction Details  --
--------------------------
✅ Predicted log(IC50): 3.4318 nM
```

**Online mode** — accepts a chemical name, ChEMBL ID, or SMILES, and additionally reports the closest matching experimental IC₅₀ value(s) from ChEMBL for comparison (requires internet access):

```bash
python run_online.py
```

```
===================================
||   PPIM-IC50PRED WEBSERVER    ||
===================================

⚗️ Enter chemical name, ChEMBL ID, SMILES, or molecular formula: CHEMBL5414626
------------------------
--  Compound Details  --
------------------------
🔎 Compound Name: Unknown
🆔 ChEMBL ID: CHEMBL5414626
🧪 SMILES: CCOC(=O)c1nc(N(C)Cc2ccccc2)c2c(C)noc2n1

--------------------------
--  Prediction Details  --
--------------------------
✅ Predicted Half maximal inhibitory concentration log(IC50): 3.4318 nM

----------------------------
--  Experimental Details  --
----------------------------
📊 Top Experimental IC50 values from ChEMBL:

#1 🎯 Target: Solute carrier family 26 member 6 (CHEMBL5465351)
   IC50: 1000.0000 nM
   pChEMBL: 6.0000
   log(IC50): 3.0000
   🔄 Difference (Predicted – Experimental): 0.4318 nM (14.39%)
```

Example compound: [CHEMBL5414626](https://www.ebi.ac.uk/chembl/explore/compound/CHEMBL5414626)

### Option B: Streamlit web app (same experience as the live demo)

From the repository root:

```bash
pip install -r requirements.txt
streamlit run run.py
```

Opens a browser interface at `http://localhost:8501` accepting a chemical name, ChEMBL ID, or SMILES string, with 2D/3D structure rendering, experimental ChEMBL comparison, and related PDB structure lookup via RCSB.

🌐 Hosted version: [ppimic50pred.onrender.com](https://ppimic50pred.onrender.com/)

## Notebook

`analysis/analysis/` contains four Jupyter notebooks with the full runtime code and results for the study's analysis pipeline (descriptor generation, feature selection, model training/cross-validation, clustering, and multitarget validation):

- `PPIM-IC50pred.ipynb`
- `PPIM-IC50pred_2nd_part.ipynb`
- `PPIM-IC50pred_3rd_part.ipynb`
- `PPIM-IC50pred_4th_part.ipynb`

## Data

`analysis/dataset/` contains the processed descriptor/bioactivity datasets for the three descriptor sets used in the study:

| Descriptor set | Training | Blind set |
|---|---|---|
| RDKit | `rdkit_train.csv` (2,760 molecules, 195 features) | `rdkit_blind.csv` (691 molecules) |
| PubChem | `pubchem_train.csv` (2,755 molecules, 17 features) | `pubchem_blind.csv` (689 molecules) |
| PaDEL | `padel_train.csv` (2,340 molecules, 1,356 features) | `padel_blind.csv` (586 molecules) |

Other files in this folder:

- `get2p2idb.xlsx`, `ALL_PPIS_dataset_new.xlsx`, `ALL_PPIS_ic50_results.csv` — source/raw dataset files underlying the descriptor sets above
- `rdkit/` — downstream RDKit-based analysis outputs (feature importances, structure- and feature-based clustering, Tanimoto similarity heatmap) supporting Section 4.4–4.5, Fig. 4

## Citation

If you use this tool or model in your research, please cite:

> Jana T, Karmakar J, Banerjee R, Saha S. Supervised Learning-Driven Prediction of Small-Molecule Modulator Bioactivity Against Protein–Protein Interactions. *Journal of Computer-Aided Molecular Design* (in press).


## License

No license file is currently specified. Since the manuscript describes this as a "freely available, open-source implementation," 

## Contact

Tanmoy Jana — tanmoy.jana@ucsf.edu
Department of Surgery, University of California San Francisco

## Funding

This work was carried out without external funding.
