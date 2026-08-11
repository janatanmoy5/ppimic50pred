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
| `analysis/dataset/` | Processed descriptor/bioactivity CSVs: `rdkit_train.csv`, `rdkit_blind.csv`, `pubchem_train.csv`, `padel_blind.csv` |
| `Manuscript/` | `Main_Manuscript_PPIM_IC50pred_TJ.docx`, `SM_PPIM_IC50pred_TJ.docx` |
| `Dockerfile` | Container definition; launches the Streamlit app via `streamlit run run.py` |
| `render.yaml` | Deployment configuration for the hosted demo (Render) |
| `requirements.txt` (repo root) | Dependencies for the Streamlit app |

> **Before submission:** `analysis/dataset/` currently has 4 of the 6 expected descriptor-set CSVs — `pubchem_blind.csv` and `padel_train.csv` are not present. Add them (or note why they're withheld) so the datasets underlying Tables 1–2 are fully traceable. The `Manuscript/` files should also be updated to match the JCAMD-formatted submission rather than an earlier draft, and `ppimic50pred/__pycache__/` should be removed from version control via `.gitignore`.

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

All analysis reported in the manuscript — descriptor generation (RDKit/PubChem/PaDEL), feature selection and SHAP analysis (Fig. 2), model training and cross-validation for RF/GB/SVR/LSTM (Table 1, Table 2, Fig. 3), structure- and feature-based clustering (Fig. 4), and the multitarget external validation (Fig. 5) — is contained in a single notebook:

```
analysis/analysis/PPIM-IC50pred.ipynb
```

(Note the repeated `analysis/analysis/` path — this looks like it may be an unintentional nested folder rather than a deliberate structure; consider flattening it to `analysis/PPIM-IC50pred.ipynb` for clarity, though either works as long as the README and repo agree.)

Given this notebook covers the full pipeline in one file, it's worth adding markdown section headers inside it that map directly to manuscript sections/figures (e.g., "## 4.1 Feature selection — Fig. 2", "## 4.6 Multitarget validation — Fig. 5") if not already present, so a reviewer can navigate directly to the code behind any given result.

## Data

The processed descriptor/bioactivity datasets are in `analysis/dataset/`:


## Citation

If you use this tool or model in your research, please cite:

> Jana T, Karmakar J, Banerjee R, Saha S. Supervised Learning-Driven Prediction of Small-Molecule Modulator Bioactivity Against Protein–Protein Interactions. *Journal of Computer-Aided Molecular Design* (in press).

*(Update with the final volume/page/DOI once the article is published.)*

## License

[No license file currently specified. Since the manuscript describes this as a "freely available, open-source implementation," add a `LICENSE` file — MIT or Apache-2.0 are common, permissive choices for academic software — so that reuse terms are unambiguous to readers and reviewers.]

## Contact

Tanmoy Jana — tanmoy.jana@ucsf.edu
Department of Surgery, University of California San Francisco

## Funding

This work was carried out without external funding.
