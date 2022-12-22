# imports
import os
import csv
import sys
from typing import List

from rdkit import Chem
from crem.crem import mutate_mol, grow_mol

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# Path to ChEMBL dataset with fragments of SC score 2
database_dir = os.path.abspath(os.path.join(root, "..", "..", "databases"))
database_path = os.path.join(database_dir, "replacements02_sc2.db")

# Generate unique molecules with replacement fragments from the database
def my_model(smiles_list: List[str], database_path: str) -> List[List[str]]:
    output_smiles = list()
    input_mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
    for mol in input_mols:
        try:
            mutation_result = list(mutate_mol(mol=mol, db_name=database_path))
            growth_result = list(grow_mol(mol=mol, db_name=database_path))
            generated_smiles = list(set(mutation_result+growth_result)) # Keep unique smiles
        except KeyError:
            generated_smiles = list(set(mutation_result)) # Ignore grow_mol result if it runs into an error
        output_smiles.append(generated_smiles)
    
    return output_smiles
    
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
    writer.writerow(["generated_molecules"]) # header
    for o in outputs:
        writer.writerow(o)