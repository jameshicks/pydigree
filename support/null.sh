#!/bin/bash
# Simulate data under a null model by gene dropping

../scripts/simulate_pedigree_data.py --template $1 \
    --chromosome $2 --freq 0 19999 0.0 \
    --method genedrop \
    --prefix null \
    --replications 1000 
