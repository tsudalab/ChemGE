from __future__ import print_function
import argparse
import threading
import nltk
import copy
import sys
import numpy as np
import cfg_util
import rdock_util
import score_util
import zinc_grammar
import time
from rdkit import Chem

import multiprocessing
import functools
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

GCFG = zinc_grammar.GCFG

def CFGtoGene(prod_rules, max_len=-1):
    gene = []
    for r in prod_rules:
        lhs = GCFG.productions()[r].lhs()
        possible_rules = [idx for idx, rule in enumerate(GCFG.productions()) if rule.lhs() == lhs]
        gene.append(possible_rules.index(r))
    if max_len > 0:
        if len(gene) > max_len:
            gene = gene[:max_len]
        else:
            gene = gene + [np.random.randint(0, 256) for _ in range(max_len-len(gene))]
    return gene

def GenetoCFG(gene):
    prod_rules = []
    stack = [GCFG.productions()[0].lhs()]
    for g in gene:
        try:
            lhs = stack.pop()
        except:
            break
        possible_rules = [idx for idx, rule in enumerate(GCFG.productions()) if rule.lhs() == lhs]
        rule = possible_rules[g%len(possible_rules)]
        prod_rules.append(rule)
        rhs = filter(lambda a: (type(a) == nltk.grammar.Nonterminal) 
                                and (str(a) != 'None'),
                                zinc_grammar.GCFG.productions()[rule].rhs())
        stack.extend(list(rhs)[::-1])
    return prod_rules


from rdkit.Chem import AllChem as Chem
def is_valid_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return smiles != '' and mol is not None and Descriptors.MolWt(mol) < 500


def selectParent(population, tournament_size=3):
    idx = np.random.randint(len(population), size=tournament_size)
    best = population[idx[0]]
    for i in idx[1:]:
        if population[i][0] > best[0]:
            best = population[i]
    return best

def crossover(p1_gene, p2_gene):
    crossover_point = np.random.choice(len(p1_gene))
    c1_gene = p1_gene[:crossover_point] + p2_gene[crossover_point:]
    c2_gene = p2_gene[:crossover_point] + p1_gene[crossover_point:]
    return (c1_gene, c2_gene)

def mutation(gene):
    idx = np.random.choice(len(gene))
    gene_mutant = copy.deepcopy(gene)
    gene_mutant[idx] = np.random.randint(0, 256)
    return gene_mutant

def canonicalize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if smiles != '' and mol is not None and mol.GetNumAtoms() > 1:
        return Chem.MolToSmiles(mol)
    else:
        return smiles

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smifile', default='250k_rndm_zinc_drugs_clean.smi')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--mu', type=int, default=32)
    parser.add_argument('--lam', type=int, default=64)
    parser.add_argument('--generation', type=int, default=1000)
    args = parser.parse_args()

    np.random.seed(args.seed)

    gene_length = 300
    p_crossover = 0

    N_mu = args.mu
    N_lambda = args.lam

    # initialize population
    seed_smiles = []
    with open(args.smifile) as f:
        for line in f:
            smiles = line.rstrip()
            seed_smiles.append(smiles)
    #for i in range(1000):
    #    smiles_list = np.random(seed_smiles, size=100, replace=False)
    #    scores = pool.map(func, smiles_list)
    #    for score, smiles in zip(scores, smiles_list)
    #        print("{},{}".format(score, smiles))
    #exit()

    start_time = time.time()

    initial_smiles = np.random.choice(seed_smiles, N_mu+N_lambda)
    initial_smiles = [s for s in initial_smiles]
    initial_genes = [CFGtoGene(cfg_util.encode(smiles), max_len=gene_length) for smiles in initial_smiles]
    initial_scores = rdock_util.score_qsub(initial_smiles)

    population = []
    for score, gene, smiles in zip(initial_scores, initial_genes, initial_smiles):
        population.append((score, smiles, gene))
    
    population = sorted(population, key=lambda x:x[0])[:N_mu]

    all_smiles = [canonicalize(p[1]) for p in population]
    all_result = [(p[0], s) for p, s in zip(population, all_smiles)]

    scores = [p[0] for p in population]
    max_score = np.max(scores)
    elapsed_time = time.time() - start_time
    print("%{},{},{}".format(0, max_score, elapsed_time))
    for p in population:
        print("{},{}".format(p[0], p[1]))

    for generation in range(args.generation):
        new_population_smiles = []
        new_population_genes = []
        for _ in range(N_lambda//2):
            p1 = population[np.random.randint(len(population))]
            p2 = population[np.random.randint(len(population))]

            p1_gene = p1[2]
            p2_gene = p2[2]

            if np.random.uniform() < p_crossover:
                c1_gene, c2_gene = crossover(p1_gene, p2_gene)
            else:
                c1_gene = mutation(p1_gene)
                c2_gene = mutation(p2_gene)

            c1_smiles = canonicalize(cfg_util.decode(GenetoCFG(c1_gene)))
            if c1_smiles != '' and c1_smiles not in all_smiles:
                new_population_smiles.append(c1_smiles)
                new_population_genes.append(c1_gene)
                all_smiles.append(c1_smiles)

            c2_smiles = canonicalize(cfg_util.decode(GenetoCFG(c2_gene)))
            if c2_smiles != '' and c2_smiles not in all_smiles:
                new_population_smiles.append(c2_smiles)
                new_population_genes.append(c2_gene)
                all_smiles.append(c2_smiles)
        new_population_scores = rdock_util.score_qsub(new_population_smiles)
        for score, gene, smiles in zip(new_population_scores, new_population_genes, new_population_smiles):
            population.append((score, smiles, gene))
            all_result.append((score, smiles))
        population = sorted(population, key=lambda x:x[0])[:N_mu]
        scores = [p[0] for p in population]
        max_score = np.max(scores)
        elapsed_time = time.time() - start_time
        print("%{},{},{}".format(generation+1, max_score,elapsed_time))
        for p in population:
            print("{},{}".format(p[0], p[1]))

    print("list of generated smiles:")
    for r in all_result:
        print("{},{}".format(r[0], r[1]))

if __name__ == "__main__":
    main()
