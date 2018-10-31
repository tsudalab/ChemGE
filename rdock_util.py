from __future__ import print_function
import hashlib
import os
import subprocess

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors


def score_qsub(smiles_list, num_docking=3):
    md5_list = []
    procs = []
    for smiles in smiles_list:
        smiles_md5 = str(hashlib.md5(smiles.encode('utf-8')).hexdigest())
        md5_list.append(smiles_md5)
        sdf_name = '{}.sdf'.format(smiles_md5)
        docking_result_file = '{}_out'.format(smiles_md5)

        # Translation from SMILES to sdf
        if smiles == '':
            mol = None
        else:
            mol = Chem.MolFromSmiles(smiles)
        try:
            if mol is not None and Descriptors.MolWt(mol) < 500:
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol)
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                fw = Chem.SDWriter(sdf_name)
                fw.write(mol)
                fw.close()

                # rdock calculation
                cmd = '''#!/bin/bash
export RBT_ROOT="$HOME/rDock_2013.1_src"
export RBT_HOME=$HOME/rDock_2013.1_src
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$RBT_ROOT/lib"
export PATH="$PATH:$RBT_ROOT/bin"
'''
                cmd += '$RBT_ROOT/bin/rbdock -r cavity.prm '\
                       '-p $RBT_ROOT/data/scripts/dock.prm '\
                       '-i {} -o {} -T 1 -n {} > /dev/null'\
                       .format(sdf_name, docking_result_file, num_docking)
                qsub_filename = 'qsub-{}.sh'.format(smiles_md5)
                with open(qsub_filename, 'w') as f:
                    f.write(cmd)
                procs.append(subprocess.Popen(
                             ['qsub -cwd -sync y {} > /dev/null'
                              .format(qsub_filename)], shell=True))
        except Exception:
            pass

    [p.wait() for p in procs]
    scores = []
    score_name = '<SCORE.INTER>'  # <SCORE> or <SCORE.INTER>
    for md5 in md5_list:
        min_score = 1e10
        path = '{}_out.sd'.format(md5)
        # find the minimum score of rdock from multiple docking results
        if os.path.exists(path):
            with open(path, 'r') as f:
                lines = f.readlines()
            isScore = False
            for line in lines:
                if isScore:
                    min_score = min(float(line), min_score)
                    isScore = False
                if score_name in line:  # next line has score
                    isScore = True
        scores.append(min_score)
    assert(len(smiles_list) == len(scores))
    return scores


def score(smiles, num_docking=3):
    smiles_md5 = str(hashlib.md5(smiles.encode('utf-8')).hexdigest())
    docking_result_file = '{}_out'.format(smiles_md5)
    sdf_name = '{}.sdf'.format(smiles_md5)
    score_name = '<SCORE.INTER>'  # <SCORE> or <SCORE.INTER>

    min_score = 1e10

    # Translation from SMILES to sdf
    if smiles == '':
        mol = None
    else:
        mol = Chem.MolFromSmiles(smiles)
    try:
        if mol is not None and Descriptors.MolWt(mol) < 500:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
            fw = Chem.SDWriter(sdf_name)
            fw.write(mol)
            fw.close()

            # rdock calculation
            cmd = '$RBT_ROOT/bin/rbdock -r cavity.prm '\
                  '-p $RBT_ROOT/data/scripts/dock.prm '\
                  '-i {} -o {} -T 1 -n {} > /dev/null'\
                  .format(sdf_name, docking_result_file, num_docking)
            path = docking_result_file+'.sd'
            if not os.path.exists(path):
                subprocess.call(cmd, shell=True)

            # find the minimum score of rdock from multiple docking results
            if os.path.exists(path):
                with open(path, 'r') as f:
                    lines = f.readlines()
                isScore = False
                for line in lines:
                    if isScore:
                        min_score = min(float(line), min_score)
                        isScore = False
                    if score_name in line:  # next line has score
                        isScore = True
    except Exception:
        pass
    return min_score
