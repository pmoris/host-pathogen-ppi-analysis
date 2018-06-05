#!/usr/bin/env bash

cd $(dirname "$0")
cd ../..

python src/scripts/pathogen_selection.py -i data/raw/ppi_data/ -n data/raw/taxdump/ -t 10292 -o data/interim/ | tee $(dirname "$0")/10292_output_pathogen_selection.log

python src/scripts/filter_and_remap.py -i data/interim/10292/ppi_data/10292-ppi-merged.tsv -t data/interim/10292/taxonid -o data/interim/10292 2>&1 | tee $(dirname "$0")/10292_output_filter_and_remap.log

# local version
# python src/scripts/filter_and_remap.py -i data/interim/10292/ppi_data/10292-ppi-merged.tsv -t data/interim/10292/taxonid -o data/interim/10292-local 2>&1 | tee $(dirname "$0")/10292_output_filter_and_remap.log

# filter interpro
python src/scripts/filter_interpro.py -i /media/pieter/Seagate\ Red\ Pieter\ Moris/workdir/iprextract/protein2ipr.v68.0.dat -p data/interim/10292/ppi_data/uniprot_identifiers.txt -o data/interim/10292/interpro/protein2ipr.dat

# filter gaf
python src/scripts/filter_gaf.py -i /media/pieter/Seagate\ Red\ Pieter\ Moris/workdir/uniprot-GOA/goa_uniprot_all.gaf -p data/interim/all/all_identifiers.txt -o data/interim/all/go_data/goa_uniprot_all.gaf

# for local version, first run
# python src/scripts/extract_all_uniprot.py -i raw/ppi_data/ -m data/raw/idmapping.dat -o data/interim/all
# then repeat filter steps.
