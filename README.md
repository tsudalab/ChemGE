# ChemGE
Molecule design with grammatical evolution.
The paper is available at ?????.

# Requirements
1. Python
1. RDKit
1. rDock

## How to set up on Ubuntu 16.04
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
```

### Download ChemGE

```
git clone https://github.com/tsudalab/ChemGE.git
cd ChemGE
```

# Usage
## Optimize J

```
python -u optimizeJ.py > log-file &
```

## Optimize docking score (multi-thread)

```
python -u optimize-rdock.py > log-file &
```

## Optimize docking score (qsub)

```
python -u optimize-rdock-qsub.py > log-file &
```

Using [CfnCluster](https://github.com/awslabs/cfncluster) is recommended.
