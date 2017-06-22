The entire NCBI taxonomy tree was retrieved via: 
```shell
$ wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
$ mkdir taxdump
$ tar -xzvf taxdump.tar.gz -C taxdump/
$ rm taxdump.tar.gz
```
The [retrieve_taxids.py script](../ppi_scripts/retrieve_taxids.py) was called using the following arguments in order to extract all the child taxids of the Herpesviridae (10292): 
```shell 
$ python3 retrieve_taxids.py taxdump 10292
``` 
