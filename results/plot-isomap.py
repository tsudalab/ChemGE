#! -*- coding: utf-8 -*-
from __future__ import print_function
import csv
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit import rdBase
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE, MDS, Isomap
import seaborn as sns

sns.set()

rdBase.DisableLog('rdApp.error')

np.random.seed(0)

# convert rdkit fingerprint to numpy array
def fp2arr(fp):
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr
mol_zinc = []
with open('zinc10000.txt', 'r') as f:
    for line in f:
        smiles = line.rstrip()
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            mol_zinc.append(mol)
        else:
            print(smiles)

mol_active = []
with open('actives_final.ism', 'r') as f:
    for line in f:
        smiles = line.split()[0]
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            mol_active.append(mol)
        else:
            print(smiles)


mols0 = []
mols100 = []
mols1000 = []
generation = -1
with open('log-rdock-population', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[0][0] == '%':
            generation += 1
            continue
        mol = Chem.MolFromSmiles(row[1])
        if mol is not None:
            if generation == 0:
                mols0.append(mol)
            if generation == 100:
                mols100.append(mol)
            if generation == 1000:
                mols1000.append(mol)

bestmol = Chem.MolFromSmiles("CP(N)C(=O)c1ccc([C-]2NCCN2C(=N)CSC2=NS(O)=CS2)s1")

fp_zinc = np.array([fp2arr(MACCSkeys.GenMACCSKeys(mol)) for mol in mol_zinc])
fp_active = np.array([fp2arr(MACCSkeys.GenMACCSKeys(mol)) for mol in mol_active])
fp_all = np.array([fp2arr(MACCSkeys.GenMACCSKeys(mol)) for mol in mol_zinc+mol_active])
fp0 = np.array([fp2arr(MACCSkeys.GenMACCSKeys(mol)) for mol in mols0])
fp100 = np.array([fp2arr(MACCSkeys.GenMACCSKeys(mol)) for mol in mols100])
fp1000 = np.array([fp2arr(MACCSkeys.GenMACCSKeys(mol)) for mol in mols1000])

fpbest = np.array([fp2arr(MACCSkeys.GenMACCSKeys(bestmol))])

isomap = Isomap(n_components=2).fit(fp_all)

X_zinc = isomap.transform(fp_zinc)
X_active = isomap.transform(fp_active)
X0 = isomap.transform(fp0)
X100 = isomap.transform(fp100)
X1000 = isomap.transform(fp1000)
Xbest = isomap.transform(fpbest)

plt.figure()
plt.scatter(X_zinc[:, 0], X_zinc[:, 1], marker='.', alpha=0.1)
plt.scatter(X_active[:, 0], X_active[:, 1], color='r', marker='o')
plt.xlim(-15, 30)
plt.ylim(-25, 25)
plt.savefig("isomap_plot_ZINC_active.pdf")

plt.figure()
plt.scatter(X_zinc[:, 0], X_zinc[:, 1], marker='.', alpha=0.1)
plt.scatter(X0[:, 0], X0[:, 1], marker='x',color='r')
plt.xlim(-15, 30)
plt.ylim(-25, 25)
plt.savefig("isomap_generation0.pdf")

plt.figure()
plt.scatter(X_zinc[:, 0], X_zinc[:, 1], marker='.', alpha=0.1)
plt.scatter(X100[:, 0], X100[:, 1], marker='x',color='r')
plt.xlim(-15, 30)
plt.ylim(-25, 25)
plt.savefig("isomap_generation100.pdf")

plt.figure()
plt.scatter(X_zinc[:, 0], X_zinc[:, 1], marker='.', alpha=0.1)
plt.scatter(X1000[:, 0], X1000[:, 1], marker='x',color='r')
plt.scatter(Xbest[:, 0], Xbest[:, 1], marker='*',color='r')
plt.xlim(-15, 30)
plt.ylim(-25, 25)
plt.savefig("isomap_generation1000.pdf")
