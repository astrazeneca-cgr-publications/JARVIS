## Current OMIM-HI gene list used:
`omim_haploinsufficient_genes.txt`


#### Get alternative list of Haploinsufficient genes: 
1. [not useful - deprecated]  
```
wget http://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/omim/attribute_set_library_crisp.gmt.gz
```
Then grep for _"autosomal dominant"_ (and maybe for _"autosomal recessive"_ too).

```
cat attribute_set_library_crisp.gmt | grep dominant | cut -f1 | sort > altern_beta_omim_haploinsufficient_genes.txt
```
----

2.
```
cat doi_ejhg.2008.111_haploinsufficiency_table.csv | cut -d ',' -f1,3,4,5,6 > doi_ejhg.2008.111_haploinsufficient_genes.csv.tmp
```

