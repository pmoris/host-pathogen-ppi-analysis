# Scripts

Taken together, these scripts can perform the full analysis pipeline:

1. `pathogen_selection.py` can extract the protein-protein interactions for a desired pathogen.
2. `filter_and_remap.py` can filter PPIs and remap identifiers to UniProt AC.
3. `filter_interpro.py` and `filter_gaf.py` can extract the relevant entries from files that map UniProt AC to InterPro/Gene Ontology terms.
4. `annotate.py` can add these annotation labels to the PPI dataset.
5. `pairwise_mining.py` can mine for pairs of annotation labels, across host-pathogen PPIs, that co-occur more frequently than expected (i.e. are associated).

# bash pipeline

Inside the `analysis-*` directories, bash scripts can be found that call these scripts in the correct order and with the required arguments and input files. These bash scripts also create log files for all generated output and errors.
