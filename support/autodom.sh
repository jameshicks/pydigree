#!/bin/bash
# Simulate data under an autosomal dominant model
# Usage: autodom.sh template.ped chromosome.template constraintfile

../scripts/simulate_pedigree_data.py --template $1 \
    --chromosome $2 --freq 0 19999 0.0 \
    --method constrained \
    --constraintfile $3 \
    --effect 0 19999 2 2 1 \
    --effect 0 19999 2 1 1 \
    --liability-threshold 0.5 \
    --prefix autodom
