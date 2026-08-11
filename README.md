# ppimic50pred 


Supervised machine learning framework for predicting small-molecule bioactivity (IC50) against protein–protein interactions (PPIs), single proteins, and cell-based assays.

This repository accompanies the manuscript "Supervised Learning-Driven Prediction of Small-Molecule Modulator Bioactivity Against Protein-Protein Interactions" (Jana T, Karmakar J, Banerjee R, Saha S; submitted to Journal of Computer-Aided Molecular Design). It provides the trained random forest model, prediction scripts, and analysis code used in that study, benchmarking RDKit, PubChem, and PaDEL molecular descriptors across four supervised regressors (random forest, gradient boosting, support vector regression, and LSTM) on a curated ChEMBL dataset of 3,451 compounds spanning 176 biological targets.
It is available at https://ppimic50pred.onrender.com/

# Overview
Input: a small-molecule SMILES string, ChEMBL ID, compound name, or molecular formula.
Output: predicted Log₁₀(IC₅₀), with an optional comparison against experimentally reported ChEMBL bioactivity values when run in online mode.
Model: random forest regressor trained on 217 RDKit descriptors (195 after feature selection), achieving R² = 0.75 / RMSE = 0.78 in cross-validation and R² = 0.74 / RMSE = 0.79 on a held-out blind set. See the manuscript for full performance benchmarks, feature-importance analysis, and known generalization limits on structurally novel targets.


# ⚙️ Installation

Clone the repository and install required dependencies:

```bash
git clone https://github.com/janatanmoy5/ppimic50pred.git
cd ppimic50pred
pip install -r requirements.txt

# ⚙️ How to Run

Conside the example of small chemical https://www.ebi.ac.uk/chembl/explore/compound/CHEMBL5414626

1. Offline mode
2. Online mode connected with CHEMBL resource
___________________________
1. Offline mode

(base) janat@SURG-6KH9X4Q-LT ppimic50pred % python run_offline.py 

=============================================
||   PPIM-IC50PRED WEBSERVER (OFFLINE)   ||
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
(base) janat@SURG-6KH9X4Q-LT ppimic50pred %

________________________________________________________
1. Online mode

(base) janat@SURG-6KH9X4Q-LT ppimic50pred % python run_online.py 

===================================
||   PPIM-IC50PRED WEBSERVER   ||
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


