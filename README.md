ppimic50pred 


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


