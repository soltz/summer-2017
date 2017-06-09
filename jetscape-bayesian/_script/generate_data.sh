#!/usr/local/bin/bash

# Run from main jetscape-bayesian directory via:
#   ./_script/generate_data.sh

# Assumes:
#  - trento data is in a file located in data/AuAu_200GeV_100k.txt

SIZE=100

# Create full data
echo "# Input parameters: to be listed here" > data/full_data.txt
echo "#       trento-stat          |      pythia-truth                 | slowjet" >> data/full_data.txt
echo "# iev mult xgj quench b Np e2   part1 pT eta phi part2 pT eta phi   Njets pT eta phi pT eta phi" >> data/full_data.txt

QUENCH_FACTOR=.05
while [[ $(echo $QUENCH_FACTOR'<1.0' | bc -l) = 1 ]]; do
    python3 src/output_pythia_trento.py -f data/AuAu_200GeV_100k.txt -o data/"0$QUENCH_FACTOR".txt -u "0$QUENCH_FACTOR" -n "$SIZE" -c -w -x data/filtered_"0$QUENCH_FACTOR".txt
    tail -n "$SIZE" data/"0$QUENCH_FACTOR".txt >> data/full_data.txt
    echo "0$QUENCH_FACTOR complete"
    QUENCH_FACTOR=$(echo 0$QUENCH_FACTOR'+0.05'| bc)
done
