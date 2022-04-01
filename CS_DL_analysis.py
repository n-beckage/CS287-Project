from rdkit.Chem import QED
from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Crippen
import numpy as np
import pandas as pd

# WHAT'S LEFT TO DO FOR ANALYSIS
# - basically put GPCR activity, smile strings, and all appropriate QED molecular features into a data frame
# - will also need to add the features necessary for geber
# - print out a csv file and you'll be good to go!


############################### Lipinksi Ro5 ##############################
# Returns True if mol passes Linpinski's rule of 5
def Ro5(mol):
    strikes=0
    if Lipinski.NumHDonors(mol)>5:
        strikes+=1
    if Lipinski.NumHAcceptors(mol)>10:
        strikes=+1
    if Descriptors.ExactMolWt(mol)>499:
        strikes+=1
    if Descriptors.MolLogP(mol)>4.99:
        strikes+=1
    if strikes>1:
        result=False
    elif strikes<=1:
        result=True
    return result

########################## Program Body ##################################
# need to load in GPCR predictions to harvest activity data
# Note that the numpy objects on this local machine are from random sample of the pubChem data, n=1000 (smaller for ease of processing)
smiles=np.load("CS_numpy_objs/smiles.npy",allow_pickle=True)
preds=np.load("CS_numpy_objs/GPCR_predictions.npy",allow_pickle=True)

# This for loop will be generating lists of features for each mol object
mols=[]
activity=[]
qed=[]
# Molar refractivity
MR=[]
# Total # of atoms
n_atom=[]
for i in np.arange(len(preds)): # len(preds) is equal to len(mols) FYI
    activity.append(preds[i][0])
    mols.append(Chem.MolFromSmiles(smiles[i]))
    qed.append(QED.default(mols[i]))
    MR.append(Crippen.MolMR(mols[i]))
    n_atom.append(Mol.GetNumAtoms(mols[i]))

# should return true
print(len(activity)==len(mols)==len(qed)==len(preds)==len(MR)==len(n_atom))

################### Calculating QED Properties #############################
### make a dataframe of QED properties for each mol
##### will probaly want to add more properties later for Ghose method

df=pd.DataFrame(columns=["MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS","QED"])

for i in np.arange(len(mols)):
   row=list(QED.properties(mols[i]))
   row.append(QED.default(mols[i]))
   df.loc[len(df.index)]=row

df['MR']=MR
df['ATOM']=n_atom
df['SMILES']=smiles
df['GPCR_act']=activity

column_order=['SMILES','GPCR_act','MR','ATOM',"MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS","QED"]
df=df.reindex(columns=column_order)
print(df)


df.to_csv("data/PubChemData.csv",index=False)
print("done")



















########################## Plotting Activity score by Ro5 pass/fail - Boxplot ######################
lip=np.array(lip) # convert to array for easy boolean indexing
activity=np.array(activity) # dido

# subset of activity scores for Ro5 passes
lip_pass=activity[lip]

# subset of activity scores for Ro5 fails
lip_fail=activity[np.invert(lip)]

# combining the data for the box plot
lip_data=[lip_pass,lip_fail]