#!/usr/bin/env python3
# # coding=utf-8
import pybel
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

def get_para(T=298, eps=None, sigma=None, lamb=None):
    f = 1 + lamb * 0.01 * (T - 298)
    eps = eps * f ** 2 * 4.184
    sigma = (2 * f) ** (-1/6) * sigma
    return eps, sigma

def get_para_298(T=298, eps=None, sigma=None, lamb=None):
    f = 1 + lamb * 0.01 * (T - 298)
    eps /= (f ** 2 * 4.184)
    sigma /= ((2 * f) ** (-1 / 6))
    return sigma / 2**(1/6), eps * 16 / 27

def get_canonical_smiles(smiles):
    py_mol = pybel.readstring("smi", smiles)
    return py_mol.write('can', opt={'n': None}).strip()

def get_rdkit_smiles(smiles):
    rdk_mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(rdk_mol)

def get_atom_numbers(smiles):
    mol = pybel.readstring("smi", smiles)
    mol.addh()
    return len(mol.atoms)

def get_heavy_atom_numbers(smiles):
    mol = pybel.readstring("smi", smiles)
    return len(mol.atoms)

def get_charge(smiles):
    mol = pybel.readstring("smi", smiles)
    return mol.charge

def get_stereo_isomer(smiles):
    rdk_mol = Chem.MolFromSmiles(smiles)
    isomers = tuple(EnumerateStereoisomers(rdk_mol))
    smiles_list = []
    for smi in sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
        smiles_list.append(smi)
    return smiles_list

def is_mol_stable(smiles):
    rdk_mol = Chem.MolFromSmiles(smiles)
    # double bond check
    atom_idx_list = []
    for bond in rdk_mol.GetBonds():
        if bond.GetBondTypeAsDouble()==2:
            atom_idx_list.append(bond.GetBeginAtomIdx())
            atom_idx_list.append(bond.GetEndAtomIdx())
    if len(set(atom_idx_list)) != len(atom_idx_list):
        return False
    else:
        return True

def is_aromatic(smiles):
    rdk_mol = Chem.MolFromSmiles(smiles)
    if rdk_mol.GetAromaticAtoms() == 0:
        return False
    else:
        return True

def remove_chirality(smiles):
    s = ''.join(smiles.split('@'))
    return get_canonical_smiles(s)

def similarity_comparison(smiles1, smiles2):
    from rdkit.Chem import AllChem as Chem
    from rdkit import DataStructs
    rdk_mol1 = Chem.MolFromSmiles(smiles1)
    fp1 = Chem.GetMorganFingerprintAsBitVect(rdk_mol1, 2)
    rdk_mol2 = Chem.MolFromSmiles(smiles2)
    fp2 = Chem.GetMorganFingerprintAsBitVect(rdk_mol2, 2)
    # print(smiles1, smiles2)
    return DataStructs.DiceSimilarity(fp1, fp2)

# a switch function of similarity_comparison
def similarity_score(smiles1, smiles2):
    a = similarity_comparison(smiles1, smiles2)
    if a < 0.6:
        return 0.
    else:
        return (a - 0.6) / 0.4

def get_similarity_score(smiles, smiles_list, cutoff=None):
    score = 0.
    for s in smiles_list:
        score += similarity_score(smiles, s)
        if cutoff != None and score > cutoff:
            return score
    return score

def is_similar(smiles, smiles_list, cutoff):
    score = 0.
    for s in smiles_list:
        score += similarity_score(smiles, s)
        if score > cutoff:
            return True
    return False