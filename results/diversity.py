from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

ms1 = []
with open("chemge-active.csv") as f:
    for line in f:
        ms1.append(Chem.MolFromSmiles(line.rstrip()))

ms2 = []
with open("zinc1000.txt") as f:
    for line in f:
        ms2.append(Chem.MolFromSmiles(line.rstrip()))

fps1 = [FingerprintMols.FingerprintMol(x) for x in ms1]
fps2 = [FingerprintMols.FingerprintMol(x) for x in ms2]
t = 0.0
for x in fps1:
    for y in fps2:
        t += 1-DataStructs.FingerprintSimilarity(x, y, metric=DataStructs.TanimotoSimilarity)
print(t/(len(fps1)*len(fps2)))
