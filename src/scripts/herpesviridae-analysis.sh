#!/usr/bin/env bash

cd $(dirname "$0")
cd ../..

python src/scripts/pathogen_selection.py -i data/raw/ppi_data/ -n data/raw/taxdump/ -t 10292 -o data/interim/ | tee $(dirname "$0")/10292_output_pathogen_selection.log

python src/scripts/filter_and_remap.py -i data/interim/10292/ppi_data/10292-ppi-merged.tsv -t data/interim/10292/taxonid -o data/interim/10292 2>&1 | tee $(dirname "$0")/10292_output_filter_and_remap.log
