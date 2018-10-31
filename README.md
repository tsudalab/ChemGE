# ChemGE
Molecule design using grammatical evolution.
The paper is available at https://arxiv.org/abs/1804.02134

The advantages of ChemGE are:

- Faster SMILES generation
- Inherent paralleism
- Novelty and diversity in designed molecules

In this repository, we provide the code used in our experiment.

1. Benchmark on druglikeness score (`optimizeJ.py`)
1. Design of high-scoring molecules for thymidine kinase (`oprimize-rdock.py`)
1. Scalability of ChemGE on parallel environment (`optimize-rdock-qsub.py`)

# Requirements
1. Python
1. RDKit
1. rDock

## How to set up on Ubuntu 16.04
Compile of rDock may fail in new compilers.
I recommend to use Ubuntu 16.04.

### Install Python libraries

```
git clone https://github.com/pyenv/pyenv.git $HOME/.pyenv
echo 'export PYENV_ROOT="$HOME/.pyenv"' >> $HOME/.bashrc
echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> $HOME/.bashrc
echo -e 'if command -v pyenv 1>/dev/null 2>&1; then\n  eval "$(pyenv init -)"\nfi' >> $HOME/.bashrc
source $HOME/.bashrc
pyenv install anaconda3-5.0.1
pyenv global anaconda3-5.0.1
git clone https://github.com/pyenv/pyenv-virtualenv.git $(pyenv root)/plugins/pyenv-virtualenv
exec "$SHELL"
conda create -c rdkit -n my-rdkit-env rdkit
source activate my-rdkit-env
pip install nltk networkx
```

### Install rDock

```
sudo apt-get update
sudo apt-get install build-essential libcppunit-dev libpopt-dev
cd /home/ubuntu/
wget https://sourceforge.net/projects/rdock/files/rDock_2013.1_src.tar.gz 
tar xf rDock_2013.1_src.tar.gz
cd rDock_2013.1_src/build/
make linux-g++-64
echo 'export RBT_ROOT="$HOME/rDock_2013.1_src"' >> $HOME/.bashrc
echo 'export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$RBT_ROOT/lib"' >> $HOME/.bashrc
echo 'export PATH="$PATH:$RBT_ROOT/bin"' >> $HOME/.bashrc
source $HOME/.bashrc
```

### Download ChemGE

```
git clone https://github.com/tsudalab/ChemGE.git
cd ChemGE
```

# Usage
Please execute in `my-rdkit-env` environment (please execute `source activate my-rdkit-env`)
## Optimize J

```
python -u optimize-J.py > log-file &
```

## Optimize docking score (multi-thread)

```
python -u optimize-rdock.py > log-file &
```

If all score is `10000000000.0`, installation of rDock may have failed.
Please check installation directory.
## Optimize docking score (qsub)
You need to set up parallel environment with qsub to execute this program.
Using [CfnCluster](https://github.com/awslabs/cfncluster) is recommended.

```
python -u optimize-rdock-qsub.py > log-file &
```

# Explanation of files

- `results/`: log files of our experiment
- `cavity.as`, `cavity.prm`, `receptor.mol2`: Information of thymidine kinase, which is required to run rDock. Detailed explanation is below.
- `250k_rndm_zinc_drugs_clean.smi`: ZINC dataset, which is from [mkusner/grammarVAE](https://github.com/mkusner/grammarVAE)
- `fpscores.pkl.gz`: Used to calculate J score, which is from [mkusner/grammarVAE](https://github.com/mkusner/grammarVAE)
- Python files: Used in experiment. A part of `zinc_grammar.py` and `cfg_util.py` is from [mkusner/grammarVAE](https://github.com/mkusner/grammarVAE).

# How to generate files required to execute rDock
Assume that rDock is installed following above step.

1. Download `receptor.pdb` and `crystal_ligand.mol2` from [KITH (DUD-E)](http://dude.docking.org/targets/kith)

2. Execute following commands
```
$ $RBT_ROOT/util/pdb2mol ./receptor.pdb  # pdb -> mol2
$ $RBT_ROOT/util/mol2sd crystal_ligand.mol2 # mol2 -> sd(f)
$ $RBT_ROOT/util/gen_prm crystal_ligand.sd receptor.mol2 > cavity.prm
$ $RBT_ROOT/build/exe/rbcavity -r cavity.prm -W
```

# License
This project is licensed under the terms of the MIT license.
