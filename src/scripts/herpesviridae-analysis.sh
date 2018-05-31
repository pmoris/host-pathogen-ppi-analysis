#!/usr/bin/env bash

cd ../..

python src/scripts/pathogen_selection.py -i data/raw/ppi_data/ -n data/raw/taxdump/ -t 10292 -o data/interim/

python src/scripts/filter_and_remap.py -i data/interim/10292/ppi_data/10292-ppi-merged.tsv -t data/interim/10292/taxonid -o data/interim/10292

