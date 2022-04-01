# the goal of this script is to test the methods of download_pubchem.py to determine to see if they work with one sdf file
import subprocess as sp
from rdkit import Chem
from rdkit.Chem import rdmolfiles as rd
from rdkit.Chem import MolToSmiles
from rdkit.Chem import QED
import numpy as np
import pandas as pd
import random


# The index interval of the pubcem sdf files
delta=500000

ind_i=1
ind_f=500000
# Don't need this since I already have a sdf file in this exerceise
# sp.call("gunzip Compound_"+f'{ind_i:09}'+"""_"""+f'{ind_f:09}'+".sdf.gz",shell=True)
# f'{ind_f:09}' - whatver number ind_f is, zeroes will added to the front of it such that the whole string is 9 characters long


################## Reading in the PubChem data #############################
### the goal here is to end up with a clean python list of rdkit Mol objects
suppl = rd.SDMolSupplier('data/Compound_000000001_000500000.sdf')
suppl_indices=np.arange(len(suppl))

# these should be the same
len(suppl_indices) # 443086
len(suppl) # 443086
# NOTE: CHANGE THE SIZE OF n HERE TO SOMETHING SMALLER IF YOU'RE GOING TO MESS WITH THIS SCRIPT OTHERWISE IT WILL TAKE FOREVER TO RUN
n=1000
rand_indices = np.random.choice(suppl_indices, size=n, replace=False)

# Now the task is to create a more manageable subset of the data by writing a for loop that will loop through ran_indices
# (so that the size is n), and add suppl[random_indices[i]] to the subset each iteration.

# you need the int() function to convert i from an object type 'numpy.int32' into a plain basic python int
# type(rand_indices[0]) # will return 'numpy.int32'
# type(int(rand_indices[0])) # will return 'int'
ranSamp=[]
for i in rand_indices:
   ranSamp.append(suppl[int(i)]) 

len(ranSamp) # should be equal to n
ranSamp[0] # should be the first item
type(ranSamp[0]) # the objects in this array are rdkit Mol objects


# This was throwing an error because  for whatever reason, rdkit couldn't make some of the odjects from our .sdf file into Mol objects
# for mol in ranSamp:
#    mol.GetNumAtoms()

# to figure out just how many or our objects were correctly converted to Mol objects, I created this simple counter:
# x = 0
# for i in ranSamp:
#   if i is not None:
#       x += 1
#       print(x)
# which returns 9982 iterations, which means that 18/10000 objects were not succesfully converted to Mol objects. Not bad.
# We'll just have to make a loop to have to remove the 18 objects of type None
theMols = []
for i in np.arange(len(ranSamp)):
   if ranSamp[i] is not None:
      theMols.append(ranSamp[i])
len(theMols) # should be a bit less than n

# This saves the mol object list as an .npy object that can be loaded into other scripts.
np.save("CS_numpy_objs/theMols.npy",theMols, allow_pickle=True)
# note that when you load in the .npy object, it is now an np.array rather than a list -> so np.save must convert a list to array.
x=np.load("CS_numpy_objs/theMols.npy", allow_pickle=True)
print(x[0])

# Now we have a list (NOT numpy array) of our complete Mol objects called theMols

################### Calculating QED Properties #############################
### make a dataframe of QED properties for each mol
##### will probaly want to add more properties later for Ghose method

# df=pd.DataFrame(columns=["MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS","QED"])

# for i in np.arange(len(theMols)):
#    row=list(QED.properties(theMols[i]))
#    row.append(QED.default(theMols[i]))
#    df.loc[len(df.index)]=row

# df.to_csv("PubChemData.csv",index=False)
# print("done")

