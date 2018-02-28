from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Crippen, MolFromSmiles, MolToSmiles
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops
from rdkit import rdBase

import numpy as np
import networkx as nx


import sascorer

rdBase.DisableLog('rdApp.error')

# from https://github.com/gablg1/ORGAN/blob/master/organ/mol_metrics.py#L83
def verify_sequence(smile):
    mol = Chem.MolFromSmiles(smile)
    return smile != '' and mol is not None and mol.GetNumAtoms() > 1

# from grammar VAE
logP_values = np.loadtxt('logP_values.txt')
SA_scores = np.loadtxt('SA_scores.txt')
cycle_scores = np.loadtxt('cycle_scores.txt')

def calc_score(smiles):
    if verify_sequence(smiles):
        molecule = MolFromSmiles(smiles)
        current_log_P_value = Descriptors.MolLogP(molecule)
        current_SA_score = -sascorer.calculateScore(molecule)
        cycle_list = nx.cycle_basis(nx.Graph(rdmolops.GetAdjacencyMatrix(molecule)))
        if len(cycle_list) == 0:
            cycle_length = 0
        else:
            cycle_length = max([ len(j) for j in cycle_list ])
        if cycle_length <= 6:
            cycle_length = 0
        else:
            cycle_length = cycle_length - 6
        current_cycle_score = -cycle_length


        current_SA_score_normalized = (current_SA_score - np.mean(SA_scores)) / np.std(SA_scores)
        current_log_P_value_normalized = (current_log_P_value - np.mean(logP_values)) / np.std(logP_values)
        current_cycle_score_normalized = (current_cycle_score - np.mean(cycle_scores)) / np.std(cycle_scores)

        score = (current_SA_score_normalized + current_log_P_value_normalized + current_cycle_score_normalized)
        return score
    else:
        raise ValueError("Error in calc_score: smiles is invalid.")