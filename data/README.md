**The following sections describe how the different data sources were obtained.**

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

* [Host-pathogen protein-protein interaction data sets](#host-pathogen-protein-protein-interaction-data-sets)
* [NCBI taxonomy data](#ncbi-taxonomy-data)
* [Gene Ontology data](#gene-ontology-data)
	* [.GAF Gene Association File](#gaf-gene-association-file)
	* [.OBO Gene Ontology file](#obo-gene-ontology-file)
* [InterPro data](#interpro-data)
	* [Output format](#output-format)

<!-- /code_chunk_output -->

----

# Host-pathogen protein-protein interaction data sets

- The HPIDB 2.0 data set was retrieved from: [http://www.agbase.msstate.edu/hpi/downloads/hpidb2.mitab.zip].[^fn1]
- The VirHostNet 2.0 data set is available from: [http://virhostnet.prabi.fr/] (release January 2017).[^fn2]
- The PHISTO dataset was obtained from [http://www.phisto.org/index.xhtml] by using the browse utility (Data update: January 19, 2017).[^fn3]
- The PSI-MI ontology .obo file was downloaded from [http://ontologies.berkeleybop.org/mi.obo] (data-version: releases/2017-06-14).[^fn4]

[^fn1]: Ammari, Mais G., Cathy R. Gresham, Fiona M. McCarthy, and Bindu Nanduri. 2016. "HPIDB 2.0: A Curated Database For Host–Pathogen Interactions". Database 2016: baw103. doi:10.1093/database/baw103.
[^fn2]: Guirimand, T., S. Delmotte, and V. Navratil. 2014. "Virhostnet 2.0: Surfing On The Web Of Virus/Host Molecular Interactions Data". Nucleic Acids Research 43 (D1): D583-D587. doi:10.1093/nar/gku1121.
[^fn3]: Durmuş Tekir, Saliha, Tunahan Çakır, Emre Ardıç, Ali Semih Sayılırbaş, Gökhan Konuk, Mithat Konuk, and Hasret Sarıyer et al. 2013. "PHISTO: Pathogen–Host Interaction Search Tool". Bioinformatics 29 (10): 1357-1358. doi:10.1093/bioinformatics/btt137.
[^fn4]: Smith, Barry, Michael Ashburner, Cornelius Rosse, Jonathan Bard, William Bug, Werner Ceusters, and Louis J Goldberg et al. 2007. "The OBO Foundry: Coordinated Evolution Of Ontologies To Support Biomedical Data Integration". Nature Biotechnology 25 (11): 1251-1255. doi:10.1038/nbt1346.

----

# NCBI taxonomy data
The entire NCBI taxonomy tree was retrieved via:
```shell
$ wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
$ mkdir taxdump
$ tar -xzvf taxdump.tar.gz -C taxdump/
$ rm taxdump.tar.gz
```

The [retrieve_taxids.py script](../src/data_prep/retrieve_taxids.py) was called using the following arguments in order to extract all the child taxids of the Herpesviridae (10292):

```shell
$ python3 retrieve_taxids.py src/data/raw/taxdump taxdump 10292
```

----

# Gene Ontology data

## .GAF Gene Association File

The gene association file containing all annotations for human (9606) and Herpesviridae (10292) proteins was retrieved using the following QuickGO query:

[https://www.ebi.ac.uk/QuickGO/annotations?format=gaf&tax=9606,10292&aspect=molecular_function,cellular_component,biological_process](https://www.ebi.ac.uk/QuickGO/annotations?format=gaf&tax=9606,10292&aspect=molecular_function,cellular_component,biological_process)

For all hosts present in the databases:
```
code 	taxid 	primary taxid 	taxname
1 	9606 	9606 	Homo sapiens
1 	10090 	10090 	Mus musculus
1 	10116 	10116 	Rattus norvegicus
1 	34839 	34839 	Chinchilla lanigera
1 	3702 	3702 	Arabidopsis thaliana
1 	60711 	60711 	Chlorocebus sabaeus
1 	6239 	6239 	Caenorhabditis elegans
1 	7955 	7955 	Danio rerio
1 	7998 	7998 	Ictalurus punctatus
1 	8364 	8364 	Xenopus tropicalis
1 	9031 	9031 	Gallus gallus
1 	9534 	9534 	Chlorocebus aethiops
1 	9541 	9541 	Macaca fascicularis
1 	9598 	9598 	Pan troglodytes
1 	9601 	9601 	Pongo abelii
1 	9606 	9606 	Homo sapiens
1 	9796 	9796 	Equus caballus
1 	9823 	9823 	Sus scrofa
1 	9913 	9913 	Bos taurus
1 	9940 	9940 	Ovis aries
1 	10292 	10292 	Herpesviridae
```

[https://www.ebi.ac.uk/QuickGO/GAnnotation?format=gaf&aspect=molecular_function%2Ccellular_component%2Cbiological_process&tax=9606%2C10090%2C10116%2C34839%2C3702%2C60711%2C6239%2C7955%2C7998%2C8364%2C9031%2C9534%2C9541%2C9598%2C9601%2C9606%2C9796%2C9823%2C9913%2C9940%2C10292&limit=-1](https://www.ebi.ac.uk/QuickGO/GAnnotation?format=gaf&aspect=molecular_function%2Ccellular_component%2Cbiological_process&tax=9606%2C10090%2C10116%2C34839%2C3702%2C60711%2C6239%2C7955%2C7998%2C8364%2C9031%2C9534%2C9541%2C9598%2C9601%2C9606%2C9796%2C9823%2C9913%2C9940%2C10292&limit=-1)

## .OBO Gene Ontology file

The gene ontology the *go-basic.obo* (daily release data-version: releases/2017-05-30) was downloaded from [http://purl.obolibrary.org/obo/go/go-basic.obo](http://purl.obolibrary.org/obo/go/go-basic.obo) .

----

# InterPro data

InterPro entries for all UniProtKBs were obtained from [https://www.ebi.ac.uk/interpro/download.html](https://www.ebi.ac.uk/interpro/download.html).[^fn5]

UniProtKB proteins | All UniProtKB proteins and the InterPro entries and individual signatures they match, in a tab-delimited format. | protein2ipr.dat.gz| gzipped 	

Due to the large size this file should be filtered using the script provided in `src/data_prep/filter_interpro.py`

[^fn5]: Philip Jones, David Binns, Hsin-Yu Chang, Matthew Fraser, Weizhong Li, Craig McAnulla, Hamish McWilliam, John Maslen, Alex Mitchell, Gift Nuka, Sebastien Pesseat, Antony F. Quinn, Amaia Sangrador-Vegas, Maxim Scheremetjew, Siew-Yit Yong, Rodrigo Lopez, and Sarah Hunter (2014). InterProScan 5: genome-scale protein function classification. Bioinformatics, Jan 2014; doi:10.1093/bioinformatics/btu031

## Output format

Information about the InterPro .tsv output format can be found here:
https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats
