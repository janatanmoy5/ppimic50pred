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

All analysis reported in the manuscript — descriptor generation (RDKit/PubChem/PaDEL), feature selection and SHAP analysis (Fig. 2), model training and cross-validation for RF/GB/SVR/LSTM (Table 1, Table 2, Fig. 3), structure- and feature-based clustering (Fig. 4), and the multitarget external validation (Fig. 5) — is contained in a single notebook:

```
analysis/analysis/PPIM-IC50pred.ipynb
```

(Note the repeated `analysis/analysis/` path — this looks like it may be an unintentional nested folder rather than a deliberate structure; consider flattening it to `analysis/PPIM-IC50pred.ipynb` for clarity, though either works as long as the README and repo agree.)

Given this notebook covers the full pipeline in one file, it's worth adding markdown section headers inside it that map directly to manuscript sections/figures (e.g., "## 4.1 Feature selection — Fig. 2", "## 4.6 Multitarget validation — Fig. 5") if not already present, so a reviewer can navigate directly to the code behind any given result.

## Data

The processed descriptor/bioactivity datasets are in `analysis/dataset/`:

| File | Descriptor set | Split |
|---|---|---|
| `rdkit_train.csv` | RDKit | Training (2,760 molecules, 195 features) |
| `rdkit_blind.csv` | RDKit | Blind set (691 molecules) |
| `pubchem_train.csv` | PubChem | Training (2,755 molecules, 17 features) |
| `padel_blind.csv` | PaDEL | Blind set (586 molecules) |
| `padel_train.csv` | PaDEL | Training (2,340 molecules, 1,356 features) |

**Missing:** `pubchem_blind.csv` (689 molecules), reported in the manuscript's Tables 1–2, is not yet present in this folder. Add it before submission so all six train/blind combinations across the three descriptor sets are reproducible from this repository, or update the manuscript's Data Availability statement to accurately reflect what is and isn't included.

The external multitarget validation set (1,528 compounds) referenced in Section 4.6/Fig. 5 is not listed above either — confirm whether it's included elsewhere in the repo (e.g., inside the notebook itself) or needs to be added.

`analysis/dataset/rdkit/` contains downstream outputs derived from the RDKit descriptor set, supporting the structure- and feature-based clustering analysis (Section 4.4–4.5, Fig. 4):

| File(s) | Description |
|---|---|
| `feature_importances.csv`, `feature_importances_plot.png` | Feature importance ranking (Fig. 2A) |
| `rdkit_cluster_0.csv` – `rdkit_cluster_4.csv` | Structure-based (Butina) cluster assignments |
| `rdkit_chemical_class_*.csv` (35 files) | Compounds grouped by the six physicochemical-property classes (size, polarity, ring type, flexibility, hydrophobicity) described in Section 4.5 |
| `rdkit_kmeans_cluster_0.csv` – `rdkit_kmeans_cluster_4.csv` | Feature-based (K-means) cluster assignments |
| `rdkit_clusters_and_classes.csv`, `rdkit_clusters_and_classes_reduced.csv`, `rdkit_clusters_with_ids.csv` | Combined cluster/class mapping tables |
| `tanimoto_similarity_clustermap.png` | Tanimoto similarity heatmap (Fig. 4B) |
| `rf_surrogate_linear_coeffs.csv` | Coefficients from a linear surrogate model fit to approximate the random forest — supplementary interpretability output not otherwise described in the manuscript |

## Citation

If you use this tool or model in your research, please cite:

> Jana T, Karmakar J, Banerjee R, Saha S. Supervised Learning-Driven Prediction of Small-Molecule Modulator Bioactivity Against Protein–Protein Interactions. *Journal of Computer-Aided Molecular Design* (in press).

*(Update with the final volume/page/DOI once the article is published.)*

## License

No license file is currently specified. Since the manuscript describes this as a "freely available, open-source implementation," add a `LICENSE` file — MIT or Apache-2.0 are common, permissive choices for academic software — so that reuse terms are unambiguous to readers and reviewers.

## Contact

Tanmoy Jana — tanmoy.jana@ucsf.edu
Department of Surgery, University of California San Francisco

## Funding

This work was carried out without external funding.
