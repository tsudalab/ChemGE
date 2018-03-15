from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.ML.Cluster import Butina
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

import csv

def ClusterFps(fps,cutoff=0.2):
    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs

ms = []
smiles = []
with open("mol-rdock-chemge-active.csv") as f:
    reader = csv.reader(f)
    for row in reader:
        s = row[1]
        smiles.append(s)
        try:
            m = Chem.MolFromSmiles(s)
        except:
            pass
        if m is not None:
            ms.append(m)

for m in ms:
    rdDepictor.Compute2DCoords(m)
fps = [AllChem.GetMorganFingerprintAsBitVect(x,2,1024) for x in ms]
clusters=ClusterFps(fps,cutoff=0.72)
print("There are {} clusters".format(len(clusters)))
for n, c in enumerate(clusters):
    idx = c[0]
    mol = ms[idx]
    Chem.Kekulize(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(400,200)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    with open("known-active-centroid-{}.svg".format(n), 'w') as f:
        f.write(svg)
