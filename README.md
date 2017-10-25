# host-pathogen-ppi-fim

In this project we utilise frequent item set mining and association rule mining to analyse and visualise the protein-protein interaction (PPI) between pathogens and their hosts (more specifically, Herpesviridae). The aim was to obtain an overview of the high-level types (strategies, life cycle stages, etc.) that these PPIs describe based on their labels.

## Contents

* [host-pathogen-ppi-fim](#host-pathogen-ppi-fim)
	* [Contents](#contents)
		* [data](#data)
		* [src](#src)
		* [notebooks](#notebooks)
		* [cytoscape](#cytoscape)
		* [poster](#poster)
* [Dependencies](#dependencies)

<!-- /code_chunk_output -->

### data
The [data directories](./data/) contain information on where the raw data and annotations were obtained. Only interim and processed data are stored in the repository.

### src
Contains python scripts to:
- Pre-process raw annotation (InterPro) and taxon id data.
- create transaction datasets of protein-protein interaction (PPI) pairs using their gene ontology and other (orthogonal) data sources (e.g. InterPro domains) as the input for frequent item set mining.
- Mine the PPI pairs for frequently occurring combinations of labels. The software to mine these item sets and association rules was developed in our research group by Danh Bui (Adrem Data Lab).

### notebooks
Jupyter notebooks that explore the data and summarise some of the findings.

### cytoscape
Cytoscape sessions that visualise the PPI network and the obtained association rules.

### poster
A preliminary version of this work was presented as a poster during the [ISMB/ECCB 2017 conference](https://www.iscb.org/ismbeccb2017).

---

# Prerequisites
- python3.6
- pandas
- numpy
- urllib3
- GO-tools https://github.com/pmoris/go-tools (The scripts expect to find a directory named "go_tools" containing these modules.)

# License
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

# Acknowledgements
My supervisors (Kris Laukens & Pieter Meysman) and other colleagues at the ADREM Data Lab (University of Antwerp) and biomina.
