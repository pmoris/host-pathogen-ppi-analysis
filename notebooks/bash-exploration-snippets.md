# Collection of `bash` snippets to quickly explore the data

## List all taxids in the PPI datasets
```bash {cmd=true}
cat <(tail -n +2  ../data/raw/ppi_data/hpidb2_March14_2017_mitab.txt | cut -f10) <(tail -n +2 ../data/raw/ppi_data/hpidb2_March14_2017_mitab.txt | cut -f11) | sed -E 's/taxid:([[:digit:]]+).*/\1/g' | sort -u | head
```

## Number of columns in data set

```bash {cmd=true}
head -1 ../data/raw/ppi_data/hpidb2_March14_2017_mitab_plus.txt | awk -F'\t' '{print NF }'
```

### What if there are some lines with more/fewer entries?
```bash {cmd=true}
awk -F'\t' '{print NF}' ../data/raw/ppi_data/hpidb2_March14_2017_mitab_plus.txt | sort -nu | tail -n 1
```

### List names as well
```bash {cmd=true}
awk -F'\t' ' { for (i = 1; i <= NF; ++i) print i, $i; exit } ' ../data/raw/ppi_data/hpidb2_March14_2017_mitab_plus.txt
```

## Retrieve data sources
```bash {cmd=true}
cat <(tail -n +2 ../data/raw/ppi_data/hpidb2_March14_2017_mitab.txt | cut -f1) <(tail -n +2 ../data/raw/ppi_data/hpidb2_March14_2017_mitab.txt | cut -f2) | sed -r 's/(^.*?):.*/\1/g' | sort -u
```

## Number of Entrez genes
```bash {cmd=true}
cat <(cut -f1 ../data/raw/ppi_data/hpidb2_March14_2017_mitab_plus.txt) <(cut -f2 ../data/raw/ppi_data/hpidb2_March14_2017_mitab_plus.txt) | grep entrez | sort -u | sed -rn 's/^.*:(.*)/\1/p' | wc -l
```

## Check overlap between data sets
```bash {cmd=false}
comm <(cut -f3 -d, ../data/raw/ppi_data/phisto_Jan19_2017.csv | sed 's/"//g' | sort -u ) <(cut -f2 ../data/raw/ppi_data/hpidb2_March14_2017_mitab.txt | sed s/uniprotkb://g | sort -u)
```
