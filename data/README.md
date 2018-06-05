**The following sections describe how the different data sources were obtained.**

<!-- @import "[TOC]" {cmd="toc" depthFrom=1 depthTo=6 orderedList=false} -->
<!-- code_chunk_output -->

* [Host-pathogen protein-protein interaction datasets](#host-pathogen-protein-protein-interaction-datasets)
	* [Important notes](#important-notes)
* [NCBI taxonomy data](#ncbi-taxonomy-data)
	* [Example output for the _Herpesviridae_ (taxid:10292)](#example-output-for-the-_herpesviridae_-taxid10292)
* [Protein identifier remapping](#protein-identifier-remapping)
* [Gene Ontology data](#gene-ontology-data)
	* [.GAF Gene Association File](#gaf-gene-association-file)
		* [Creating .gaf file without specifying host](#creating-gaf-file-without-specifying-host)
		* [Creating a .gaf file from a larger database instead of via QuickGO](#creating-a-gaf-file-from-a-larger-database-instead-of-via-quickgo)
	* [.OBO Gene Ontology file](#obo-gene-ontology-file)
* [InterPro data](#interpro-data)
	* [Output format](#output-format)

<!-- /code_chunk_output -->

----

# Host-pathogen protein-protein interaction datasets

- The HPIDB 3.0 dataset was retrieved in PSI-MITAB 2.5 format from: [hpidb.igbb.msstate.edu/downloads/hpidb2.mitab.zip].[^fn1]
- The VirHostNet 2.0 dataset is available in PSI-MITAB 2.5 format from: [http://virhostnet.prabi.fr/] (release January 2018).[^fn2]
- The IntAct PSI-MI TAB dataset was retrieved from the computationally maintained dataset: Virus [https://www.ebi.ac.uk/intact/downloads] (release date 22-03-2018 on FTP site) and was downloaded in the MI-TAB 2.5 format [https://www.ebi.ac.uk/intact/query/annot:%22dataset:virus%22?conversationContext=4](https://www.ebi.ac.uk/intact/query/annot:%22dataset:virus%22?conversationContext=4).
- The BIOGRID dataset was retrieved in PSI-MITAB 2.5 format from: [https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.4.160/](https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.4.160/) (release 3.4.160 April 23th, 2018).[^fn3]
- The PHISTO dataset was obtained from [http://www.phisto.org/index.xhtml] by using the browse utility and selecting all viruses (Data update: June 12, 2017).[^fn4]
- The PSI-MI ontology .obo file was downloaded from [http://ontologies.berkeleybop.org/mi.obo] (data-version: releases/2018-05-08).[^fn5] This file was required to parse the PHISTO dataset.

[^fn1]: Ammari, Mais G., Cathy R. Gresham, Fiona M. McCarthy, and Bindu Nanduri. 2016. "HPIDB 2.0: A Curated Database For Host–Pathogen Interactions". Database 2016: baw103. doi:10.1093/database/baw103.
[^fn2]: Guirimand, Thibaut, Stéphane Delmotte, and Vincent Navratil. 2014. "Virhostnet 2.0: Surfing On The Web Of Virus/Host Molecular Interactions Data". Nucleic Acids Research 43 (D1): D583-D587. doi:10.1093/nar/gku1121.
[^fn3]: Stark, C. 2006. "Biogrid: A General Repository For Interaction Datasets". Nucleic Acids Research 34 (90001): D535-D539. doi:10.1093/nar/gkj109.
[^fn4]: Durmuş Tekir, Saliha, Tunahan Çakır, Emre Ardıç, Ali Semih Sayılırbaş, Gökhan Konuk, Mithat Konuk, and Hasret Sarıyer et al. 2013. "PHISTO: Pathogen–Host Interaction Search Tool". Bioinformatics 29 (10): 1357-1358. doi:10.1093/bioinformatics/btt137.
[^fn5]: Smith, Barry, Michael Ashburner, Cornelius Rosse, Jonathan Bard, William Bug, Werner Ceusters, and Louis J Goldberg et al. 2007. "The OBO Foundry: Coordinated Evolution Of Ontologies To Support Biomedical Data Integration". Nature Biotechnology 25 (11): 1251-1255. doi:10.1038/nbt1346.

More information about the PSI-MITAB format can be found at the [PSIQUIC documentation](https://psicquic.github.io/MITAB25Format.html).

## Important notes

The HPIDB dataset contains multiple identifiers in the unique reference columns for protein A and B (e.g. `dip:DIP-928N|refseq:NP_001188|uniprotkb:Q13323`).

----

# NCBI taxonomy data

The entire NCBI taxonomy tree was retrieved via:

```shell
$ wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
$ mkdir taxdump
$ tar -xzvf taxdump.tar.gz -C taxdump/
$ rm taxdump.tar.gz
```

The [pathogen_selection.py script](../src/scripts/pathogen_selection.py) was called using the following arguments in order to extract all the child taxonIDs of the _Herpesviridae_ (10292) ([NCBI Taxonomy link](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Tree&id=10292&lvl=3&p=mapview&p=has_linkout&p=blast_url&p=genome_blast&keep=1&srchmode=1&unlock)):

## Example output for the _Herpesviridae_ (taxid:10292)

```shell
python src/scripts/pathogen_selection.py -i data/raw/ppi_data/ -n data/raw/taxdump/ -t 10292 -o data/interim/
```

> Size of dataframe (30332, 19)
> origin
> BIOGRID-ALL-3.4.160.mitab        1427
> hpidb2-19.02.2018.mitab          9129
> intact-virus-22.03.2018.mitab    4529
> phi_data.csv                     8118
> virhostnet-01.2018.mitab         7129
> dtype: int64

This script also extracts a list of all associated proteins' taxa.

> The following taxon IDs were found associated with 10292:

> 9598 Pan troglodytes 2
> 10036 Mesocricetus auratus 1
> 9796 Equus caballus 15
> 9601 Pongo abelii 8
> 11926 Human T-cell lymphotrophic virus type 1 (strain ATK) 3
> 648242 Adeno-associated virus 2 Srivastava/1982 12
> 34839 Chinchilla lanigera 2
> 3702 Arabidopsis thaliana 3
> 10090 Mus musculus 157
> 375286 Janthinobacterium sp. Marseille 9
> 8364 Xenopus tropicalis 2
> 9534 Chlorocebus aethiops 5
> 559292 Saccharomyces cerevisiae S288C 1
> 9913 Bos taurus 62
> 748671 Lactobacillus crispatus ST1 6
> 7955 Danio rerio 2
> 10684 Bacillus phage PBS2 1
> 9541 Macaca fascicularis 4
> 7998 Ictalurus punctatus 2
> 60711 Chlorocebus sabaeus 2
> 10116 Rattus norvegicus 127
> 83333 Escherichia coli K-12 4
> 9031 Gallus gallus 36
> 258594 Rhodopseudomonas palustris CGA009 1
> 6239 Caenorhabditis elegans 4
> 83334 Escherichia coli O157:H7 1
> 9606 Homo sapiens 25784
> 9823 Sus scrofa 15
> 9940 Ovis aries 2

|Dataset|intra-virus interactions|intra-host interactions|intra-pathogen interactions (different species)|pathogen column|
| --- | --- | --- | --- | --- |
|BIOGRID-ALL-3.4.160.mitab|17|0|0|xref_A + xref_B|
|hpidb2-19.02.2018.mitab|0|0|0|xref_B|
|intact-virus-22.03.2018.mitab|1241|0|476|xref_A + xref_B|
|virhostnet-01.2018.mitab|2801|0|624|xref_A + xref_B|
|phi_data.csv|0|0|0|xref_B|

----

# Protein identifier remapping

The `interim/pathogen-taxonid/mapping` directory contains files that remap protein identifiers to UniProt accession numbers. These can be generated using the `filter_and_remap.py` script. _WARNING_: it seems that the EBI mapping services returns an empty list every now and then. This could be due to too many requests being submitted. No warning is raised however...

For example, running the same script a few times in a row:

```
# first and second run
refseq: 0 out of identifiers 10 were succesfully remapped to UniProt accession numbers.

# third run
refseq: 2 out of identifiers 10 were succesfully remapped to UniProt accession numbers.
```

----

# Gene Ontology data

## .GAF Gene Association File

The gene association file containing all annotations for human (9606) and _Herpesviridae_ (10292) proteins can retrieved using the following QuickGO query, but due to the download limit opposed by QuickGo it cannot be downloaded in its entirety of it's too large.

~~[https://www.ebi.ac.uk/QuickGO/annotations?format=gaf&tax=9606,10292&aspect=molecular_function,cellular_component,biological_process](https://www.ebi.ac.uk/QuickGO/annotations?format=gaf&tax=9606,10292&aspect=molecular_function,cellular_component,biological_process)~~

**NOTE** url query types were changed! The new format for taxid and descendants follows this style:

[https://www.ebi.ac.uk/QuickGO/annotations?format=gaf&aspect=molecular_function,cellular_component,biological_process&taxonId=10292,9606&taxonUsage=descendants](https://www.ebi.ac.uk/QuickGO/annotations?format=gaf&aspect=molecular_function,cellular_component,biological_process&taxonId=10292,9606&taxonUsage=descendants)

----

### Creating .gaf file without specifying host

To retrieve all hosts present in a PPI database, the following method can be used:

* Use the `retrieve_taxid` python script to generate a list of all pathogen taxonids of interest, as described [above](#ncbi-taxonomy-data).
* Retrieve all taxids from all databases and omit the already retrieved pathogen ids:

		grep -v -f <(sed 's/|.*//g' ../../interim/child_taxids_of_10292.txt | sort -u) <( { cut -f10 *.mitab; cut -f11 *.mitab;} | grep -Po "^taxid:\d+" | sed 's/.*://g' | sort -u) > host_taxids.txt

To find the correct column numbers, use this one-liner `awk -F$'\t' ' { for (i = 1; i <= NF; ++i) print i, $i; exit } ' intact-virus-22.03.2018.mitab`

Problem: this list will now also include all other pathogens in your datasets. Solution: do this process on an already trimmed interaction file, for example, after selecting entries containing your pathogen of interest. This can also easily be done in pandas at this stage.

*_Note_: it is much easier to simply define the two taxon groups of interest beforehand.*

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

----

### Creating a .gaf file from a larger database instead of via QuickGO

To avoid the download limits imposed by QuickGO, a gene association file can be created by extracting the relevant entries from the entire UniProt gene association file, available at the Gene Ontology Annotation consortium: https://www.ebi.ac.uk/GOA/downloads.

	UniProtKB 	Gene association 	UniProtKB ID mapping | Readme

Be warned: unpacking this file takes ~100 GB of space.

Additional information about this file format is available in the readme: ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/README

Version information

> !Generated: 2018-04-23 14:22
> !GO-version: http://purl.obolibrary.org/obo/go/releases/2018-04-20/extensions/go-plus.owl

Due to the large size of this file, the relevant entries should be filtered prior to downstream analysis. A filtering script is provided in `src/data_prep/filter_gaf.py`. The scripts expects a unique list of protein identifiers (UniProt Accession Numbers). This file is automatically created by running the `pathogen_selection.py` and `filter_and_remap.py` scripts. This will ensure that identifiers have been remapped to UniProt ACs.

`python src/scripts/filter_gaf.py -i /data/raw/goa_uniprot_all.gaf -p data/interim/10292/ppi_data/uniprot_identifiers.txt -o data/interim/10292/go_data/goa_uniprot_10292.gaf`

Alternatively, the gene associations for all proteins in a given set of (unfiltered) mitab ppi files can be created manually by first invoking the `extract_all_uniprot.py` script. This script allows extracting and remapping all identifiers from a directory with `.mitab` files, prior to any selection of pathogens or filtering. Naturally, this will still result in rather large file. Moreover, because remapping using the online service would take a long time or even fail because of the larger number of requests that would be sent, this script expects a local `idmapping.dat` file (also available at EBI GOA).

`python src/scripts/extract_all_uniprot.py -i raw/ppi_data/ -m data/raw/idmapping.dat -o data/interim/all`

## .OBO Gene Ontology file

The gene ontology the *go-basic.obo* (data-version: releases/2018-05-07) was downloaded from [http://purl.obolibrary.org/obo/go/go-basic.obo](http://purl.obolibrary.org/obo/go/go-basic.obo) .

----

# InterPro data

InterPro entries for all UniProtKBs were obtained from [https://www.ebi.ac.uk/interpro/download.html](https://www.ebi.ac.uk/interpro/download.html).[^fn5] Version 68.0 was used (see ftp site). For more information see: [https://www.ebi.ac.uk/interpro/user_manual.html](https://www.ebi.ac.uk/interpro/user_manual.html).

    UniProtKB proteins | All UniProtKB proteins and the InterPro entries and individual signatures they match, in a tab-delimited format. | protein2ipr.dat.gz| gzipped

Due to the large size of this file, the relevant entries should be filtered prior to downstream analysis. A filtering script is provided in `src/data_prep/filter_interpro.py`. The scripts expects a unique list of protein identifiers (UniProt Accession Numbers). This file is automatically created by running the `pathogen_selection.py` and `filter_and_remap.py` scripts. This will ensure that identifiers have been remapped to UniProt ACs.

`python src/scripts/filter_interpro.py -i /data/raw/protein2ipr.v68.0.dat -p data/interim/10292/ppi_data/uniprot_identifiers.txt -o data/interim/10292/interpro/protein2ipr.dat`

Alternatively, a list can be obtained using a bash script as follows. Note however, that this list might contain non UniProt identifiers that could still potentially be remapped, and thus, no InterPro annotations will be obtained for them.

`cat <(tail -n +2 data/raw/ppi_data/*.mitab | cut -f1) <(tail -n +2 data/raw/ppi_data/*.mitab | cut -f2) | sed -r 's/^[^:]*?:(.*)/\1/g' | sort -u | grep -v "|"`

To alleviate this issue, another script was created `extract_all_uniprot.py` that allows extracting and remapping all identifiers from a directory with `.mitab` files, prior to any selection of pathogens or filtering. Naturally, this will still result in rather large file. Moreover, because remapping using the online service would take a long time or even fail because of the larger number of requests that would be sent, this script expects a local `idmapping.dat` file (also available at EBI GOA).

`python src/scripts/extract_all_uniprot.py -i raw/ppi_data/ -m data/raw/idmapping.dat -o data/interim/all`

[^fn5]: Philip Jones, David Binns, Hsin-Yu Chang, Matthew Fraser, Weizhong Li, Craig McAnulla, Hamish McWilliam, John Maslen, Alex Mitchell, Gift Nuka, Sebastien Pesseat, Antony F. Quinn, Amaia Sangrador-Vegas, Maxim Scheremetjew, Siew-Yit Yong, Rodrigo Lopez, and Sarah Hunter (2014). InterProScan 5: genome-scale protein function classification. Bioinformatics, Jan 2014; doi:10.1093/bioinformatics/btu031

## Output format

Information about the InterPro .tsv output format can be found here:
https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats
