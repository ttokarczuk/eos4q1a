# imports
import csv
import os
import sys
from typing import List

import numpy as np
import rdkit
from crem.crem import grow_mol, mutate_mol
from rdkit import DataStructs
from rdkit.Chem import MolFromSmiles, rdMolDescriptors
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import pairwise_distances_argmin_min

MAX_RETURNED_MOLS = 100  # Cap molecules to be returned

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# Path to ChEMBL dataset with fragments of SC score 2
database_dir = os.path.abspath(os.path.join(root, "..", "..", "databases"))
database_path = os.path.join(database_dir, "replacements02_sc2.db")


def generate_fingerprint(mol: rdkit.Chem.rdchem.Mol, nbits, radius=3) -> np.ndarray:
    """
    Return Morgan Fingerprint for a molecule as a fixed shape Numpy array.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit Chem.rdchem.Mol object

    Returns:
        np.ndarray: Morgan Fingerprint as numpy.ndarray
    """
    arr = np.zeros((1, nbits), dtype=np.int32)
    fp = rdMolDescriptors.GetHashedMorganFingerprint(mol, radius=radius, nBits=nbits)
    assert fp.GetLength() == nbits
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr.reshape(1, nbits)


def sample_molecules(smiles: List[str]) -> List[str]:
    """
    Samples generated molecules into 100 clusters by fitting a K-Means clustering model
    and returns the molecules closest to each cluster centre.
    Args:
        mols (List[str]): List of generated molecule smiles

    Returns:
        List[str]: List of 100 sampled molecule smiles
    """
    features = 2048
    max_sampled = 2000
    sampled_smiles: List[str] = []
    num_smiles = len(smiles)

    # Randomly choose 2000 molecules to sample with equal probability
    # This is done to ensure this sampling does not become a performance
    # bottleneck.
    if num_smiles > max_sampled:
        smiles = np.random.choice(smiles, max_sampled).tolist()

    mol_fps = np.zeros((len(smiles), features))
    for iter, smi in enumerate(smiles):
        mol = MolFromSmiles(smi)
        mol_fps[iter] = generate_fingerprint(mol, features)

    mbkm_estimator = MiniBatchKMeans(n_clusters=100)  # Default batch size: 1024
    mbkm_estimator.fit(mol_fps)
    closest, _ = pairwise_distances_argmin_min(
        mbkm_estimator.cluster_centers_, mol_fps, metric="euclidean"
    )
    for idx in closest:
        sampled_smiles.append(smiles[idx])
    return sampled_smiles


# Generate unique molecules with replacement fragments from the database
def my_model(smiles_list: List[str], database_path: str) -> List[List[str]]:
    output_smiles = list()
    input_mols = [MolFromSmiles(smi) for smi in smiles_list]
    for mol in input_mols:
        try:
            mutation_result = list(mutate_mol(mol=mol, db_name=database_path, ncores=2))
        except KeyError:
            mutation_result = []  # Nothing generated
        try:
            growth_result = list(grow_mol(mol=mol, db_name=database_path, ncores=2))
        except KeyError:
            growth_result = []  # Nothing generated

        generated_smiles = list(
            set(mutation_result + growth_result)
        )  # Keep unique smiles

        num_generated_smiles = len(generated_smiles)

        # Keep consistent structure in generated data
        if num_generated_smiles > 100:
            generated_smiles = sample_molecules(smiles=generated_smiles)
        elif num_generated_smiles < 100:
            generated_smiles = generated_smiles + [" "]*(100-num_generated_smiles)

        output_smiles.append(generated_smiles)

    return output_smiles

# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# run model
outputs = my_model(smiles_list, database_path)

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow([f"gen_mol_{i}" for i in range(100)])  # header
    for o in outputs:
        writer.writerow(o)
