#!/bin/sh
. /local/env/envconda.sh


WD=$(pwd)
export CONDA_ALWAYS_YES="true"

# init env
conda create -p $WD"/.env" python=3.10

# activate env to install packages
conda activate $WD"/.env"

# installing required python packages
python setup.py install

unset CONDA_ALWAYS_YES
conda deactivate
