#!/usr/bin/env bash

cd $(dirname "$0")
cd ../../..

# extract all identifiers
python src/scripts/extract_all_uniprot.py -i data/raw/ppi_data -o data/interim/all/ 2>&1 | tee logs/all_extract_all_uniprot.output

# filter interpro
python src/scripts/filter_interpro.py -i /media/pieter/Seagate\ Red\ Pieter\ Moris/workdir/iprextract/protein2ipr.v68.0.dat -p data/interim/all/all_identifiers.txt -o data/interim/all/interpro/protein2ipr.dat 2>&1 | tee logs/all_filter_interpro.output

# filter gaf
python src/scripts/filter_gaf.py -i /media/pieter/Seagate\ Red\ Pieter\ Moris/workdir/uniprot-GOA/goa_uniprot_all.gaf -p data/interim/all/all_identifiers.txt -o data/interim/all/go_data/goa_uniprot.gaf 2>&1 | tee logs/all_filter_gaf.output
