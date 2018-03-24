from __future__ import print_function
import argparse
import copy
import multiprocessing
import nltk
import time

import numpy as np
from rdkit import Chem
from rdkit import rdBase

import cfg_util
import rdock_util
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


def selectParent(population, tournament_size=3):
    idx = np.random.randint(len(population), size=tournament_size)
    best = population[idx[0]]
    for i in idx[1:]:
        if population[i][0] > best[0]:
            best = population[i]
    return best


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

    N_mu = args.mu
    N_lambda = args.lam

    # initialize population
    seed_smiles = []
    with open(args.smifile) as f:
        for line in f:
            smiles = line.rstrip()
            seed_smiles.append(smiles)

    start_time = time.time()

    pool = multiprocessing.Pool()

    initial_smiles = np.random.choice(seed_smiles, N_mu+N_lambda)
    initial_smiles = [s for s in initial_smiles]
    print(initial_smiles)
    initial_genes = [CFGtoGene(cfg_util.encode(smiles), max_len=gene_length)
                     for s in initial_smiles]
    print("initial score caluculating")
    initial_scores = pool.map(rdock_util.score, initial_smiles)
    print("initial score caluculated")

    population = []
    for score, gene, smiles in zip(initial_scores, initial_genes,
                                   initial_smiles):
        population.append((score, smiles, gene))

    population = sorted(population, key=lambda x: x[0])[:N_mu]

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
        for _ in range(N_lambda):
            p = population[np.random.randint(len(population))]
            p_gene = p[2]
            c_gene = mutation(p_gene)

            c_smiles = canonicalize(cfg_util.decode(GenetoCFG(c_gene)))
            if c_smiles != '' and c_smiles not in all_smiles:
                new_population_smiles.append(c_smiles)
                new_population_genes.append(c_gene)
                all_smiles.append(c_smiles)

        new_population_scores = pool.map(rdock_util.score,
                                         new_population_smiles)
        for score, gene, smiles in zip(new_population_scores,
                                       new_population_genes,
                                       new_population_smiles):
            population.append((score, smiles, gene))
            all_result.append((score, smiles))
        population = sorted(population, key=lambda x: x[0])[:N_mu]
        scores = [i[0] for i in population]
        min_score = np.min(scores)
        elapsed_time = time.time() - start_time
        print("%{},{},{}".format(generation+1, min_score, elapsed_time))
        for p in population:
            print("{},{}".format(p[0], p[1]))

    print("list of generated smiles:")
    for r in all_result:
        print("{},{}".format(r[0], r[1]))

if __name__ == "__main__":
    main()
