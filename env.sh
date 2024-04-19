#!/bin/sh
. /local/env/envconda.sh

desired_env="pancat"

if conda env list | grep -q "\b$desired_env\b"; then
    if conda info --envs | grep -q "^\s*$desired_env\b"; then
        echo "Error: The '$desired_env' environment should be active."
        echo "Please activate the '$desired_env' environment using:"
        echo "conda activate $desired_env"
        read -p "Do you want to activate environment $desired_env? ([Y]/n): " choice
        if [[ "$choice" == "y" || -z "$choice" ]]; then
            conda activate $desired_env
        else
            exit 1
        fi
    else
        echo "The '$desired_env' environment is currently active."
    fi
else
    read -p "The '$desired_env' environment does not exist. Do you want to create it? ([Y]/n): " choice
    if [[ ${choice:-y} == y ]]; then
        export CONDA_ALWAYS_YES="true"
        conda create -n $desired_env python=3.11
        conda activate $desired_env
        python setup.py install
        unset CONDA_ALWAYS_YES
        conda deactivate
        echo "The '$desired_env' environment has been created. Please activate it using:"
        echo "conda activate $desired_env"
    else
        exit 1
    fi
fi
