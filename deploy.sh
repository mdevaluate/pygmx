#!/bin/bash -l

module purge
module load gromacs/5.1.2
module load mdevaluate/dev

python setup.py install --prefix=/data/niels/modules/mdevaluate-dev
