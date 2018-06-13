#!/usr/bin/env bash

cd $(dirname "$0")
cd ../..

python src/scripts/pathogen_selection.py -i data/raw/ppi_data/ -n data/raw/taxdump/ -t 10292 -o data/interim/ 2>&1 | tee $(dirname "$0")/10292_output_pathogen_selection.log

python src/scripts/filter_and_remap.py -i data/interim/10292/ppi_data/10292-ppi-merged.tsv -t data/interim/10292/taxonid -o data/interim/10292 2>&1 | tee $(dirname "$0")/10292_output_filter_and_remap.log

# local version of remapping
# python src/scripts/filter_and_remap.py -i data/interim/10292/ppi_data/10292-ppi-merged.tsv -t data/interim/10292/taxonid -o data/interim/10292-local 2>&1 | tee $(dirname "$0")/10292_output_filter_and_remap.log

# filter interpro
python src/scripts/filter_interpro.py -i /media/pieter/Seagate\ Red\ Pieter\ Moris/workdir/iprextract/protein2ipr.v68.0.dat -p data/interim/10292/ppi_data/uniprot_identifiers.txt -o data/interim/10292/interpro/protein2ipr.dat 2>&1 | tee $(dirname "$0")/10292_interpro_filter.log

# filter gaf
python src/scripts/filter_gaf.py -i /media/pieter/Seagate\ Red\ Pieter\ Moris/workdir/uniprot-GOA/goa_uniprot_all.gaf -p data/interim/10292/ppi_data/uniprot_identifiers.txt -o data/interim/10292/go_data/goa_uniprot.gaf 2>&1 | tee $(dirname "$0")/10292_gaf_filter.log

# to create mapping file for all ppi identifiers (instead of a pathogen group), first create a list of identifiers
python src/scripts/extract_all_uniprot.py -i data/raw/ppi_data/ -o data/interim/all 2>&1 | tee $(dirname "$0")/all_online_extract_uniprot.log
# then repeat filter steps.
python src/scripts/filter_interpro.py -i /media/pieter/Seagate\ Red\ Pieter\ Moris/workdir/iprextract/protein2ipr.v68.0.dat -p data/interim/all/all_identifiers.txt -o data/interim/all/interpro/protein2ipr.dat 2>&1 | tee $(dirname "$0")/all_online_interpro_filter.log
python src/scripts/filter_gaf.py -i /media/pieter/Seagate\ Red\ Pieter\ Moris/workdir/uniprot-GOA/goa_uniprot_all.gaf -p data/interim/all/all_identifiers.txt -o data/interim/all/go_data/goa_uniprot.gaf 2>&1 | tee $(dirname "$0")/all_online_gaf_filter.log

# annotate ppi
python src/scripts/annotate.py -i data/interim/10292/ppi_data/ppi-filter-remap.tsv -b data/raw/go_data/go.obo -g data/interim/10292/go_data/goa_uniprot.gaf -p data/interim/10292/interpro/protein2ipr.dat -o data/interim/10292/ppi_data/ppi-annotations.tsv 2>&1 | tee $(dirname "$0")/10292_output_annotate.log