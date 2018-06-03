'''
@summary: Implementation of the NSGA-II algorithm in Python.
@version: 1.0
@since: 2011-01-10
@author: Marcelo Pita, http://marcelopita.wordpress.com
@contact: marcelo.souza.pita <at> gmail.com
@copyright: Copyright 2011 Marcelo Pita
@license:

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
from __future__ import print_function
import argparse
import copy
import random
import nltk
import random
import threading

import numpy as np
from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import Descriptors
from rdkit.Chem import QED

from nsga2 import Solution
from nsga2 import NSGAII

import cfg_util
import score_util
import zinc_grammar

rdBase.DisableLog('rdApp.error')
GCFG = zinc_grammar.GCFG

def CFGtoGene(prod_rules, max_len=-1):
    gene = []
    for r in prod_rules:
        lhs = GCFG.productions()[r].lhs()
        possible_rules = [idx for idx, rule in enumerate(GCFG.productions())
                          if rule.lhs() == lhs]
        gene.append(possible_rules.index(r))
    if max_len > 0:
        if len(gene) > max_len:
            gene = gene[:max_len]
        else:
            gene = gene + [np.random.randint(0, 256)
                           for _ in range(max_len-len(gene))]
    return gene


def GenetoCFG(gene):
    prod_rules = []
    stack = [GCFG.productions()[0].lhs()]
    for g in gene:
        try:
            lhs = stack.pop()
        except Exception:
            break
        possible_rules = [idx for idx, rule in enumerate(GCFG.productions())
                          if rule.lhs() == lhs]
        rule = possible_rules[g % len(possible_rules)]
        prod_rules.append(rule)
        rhs = filter(lambda a: (type(a) == nltk.grammar.Nonterminal)
                     and (str(a) != 'None'),
                     zinc_grammar.GCFG.productions()[rule].rhs())
        stack.extend(list(rhs)[::-1])
    return prod_rules

seed_smiles = []
with open('250k_rndm_zinc_drugs_clean.smi') as f:
    for line in f:
        smiles = line.rstrip()
        seed_smiles.append(smiles)

class ChemGE(Solution):
    '''
    Solution for the Molecule optimization.
    '''
    def __init__(self):
        '''
        Constructor.
        '''
        Solution.__init__(self, 2)

        s = np.random.choice(seed_smiles)
        self.attributes = CFGtoGene(cfg_util.encode(s), max_len=300)
        
        self.evaluate_solution()
        
    def evaluate_solution(self):
        '''
        Implementation of method evaluate_solution() for T1 function.
        '''
        self.smiles = cfg_util.decode(GenetoCFG(self.attributes))
        mol = Chem.MolFromSmiles(self.smiles)

        try:
            self.objectives[0] = Descriptors.MolLogP(mol)
            self.objectives[1] = QED.qed(mol)
        except:
            self.objectives[0] = -100
            self.objectives[1] = -100

        
    def crossover(self, other):
        '''
        Crossover of T1 solutions.
        '''
        child_solution = ChemGE()
        
        for i in range(300):
            child_solution.attributes[i] = self.attributes[i]

        return child_solution
    
    def mutate(self):
        '''
        Mutation of T1 solution.
        '''
        self.attributes[random.randint(0, len(self.attributes)-1)] = random.randint(0, 256)
    
if __name__ == '__main__':
    nsga2 = NSGAII(2)
    
    P = []
    for i in range(500):
        P.append(ChemGE())
    
    nsga2.run(P, 50, 100)
    
    csv_file = open('nsga2_out.csv', 'w')
    
    for i in range(len(P)):
        csv_file.write(P[i].smiles + str(P[i].objectives[0]) + ", " + str(P[i].objectives[1]) + "\n")
        
    csv_file.close()
