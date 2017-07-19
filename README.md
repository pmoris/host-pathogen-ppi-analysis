# host-pathogen-ppi-fim

Python scripts to create transaction datasets of protein-protein interaction (PPI) pairs using their gene ontology and other (orthogonal) data sources as the input for frequent item set mining.

The data directories contain information on where the PPI datasets and annotations were obtained.

---
# Dependencies
- pandas
- numpy
- GO tools https://github.com/pmoris/go-tools (The scripts expect to find a directory named "go_tools" containing these modules.)
