# imports
import os
import csv
import sys
from typing import List

from rdkit import Chem
from crem.crem import mutate_mol

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# checkpoints directory
database_dir = os.path.abspath(os.path.join(root, "..", "..", "databases"))

# read checkpoints (here, simply an integer number: 42)
database_path = os.path.join(database_dir, "replacements02_sc2.db")

# model to be run (here, calculate the Molecular Weight and add ckpt (42) to it)
def my_model(smiles_list: List[str], database_path: str) -> List[List[str]]:
    generated_smiles = list()
    input_mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
    for mol in input_mols:
        generated_smiles.append(list(mutate_mol(mol=mol, db_name=database_path)))
    
    return generated_smiles
    
# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader) # skip header
    smiles_list = [r[0] for r in reader]
    
# run model
outputs = my_model(smiles_list, database_path)

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["value"]) # header
    for o in outputs:
        writer.writerow(o)