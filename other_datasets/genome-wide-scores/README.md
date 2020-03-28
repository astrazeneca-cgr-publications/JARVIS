# For variant annotation (pathogenic/benign)

Need to provide for each `score` BED files with the intersection of the score's values with pathogenic and benign variants from a certain `resource_set`.

The naming standard for the files is:
```
genome-wide-scores/${score}/${score}.${resource_set}_[pathogenic|benign].bed
```

e.g.
cadd.hgmd_pathogenic.bed
cadd.clinvar_pathogenic.bed
cadd.clinvar_benign.bed

Available genome-wide scores
----------------------------

- cadd
- orion
- phyloP46way
- phyloP100way
- phastCons46way
- phastCons100way
- dann
- gerp

Download scores
---------------

## PhyloP
- Download phyloP 46- and 100-way
```
./download_phyloP.sh 46
./download_phyloP.sh 100
```


## PhastCons 46- and 100- way
- Download phastCons 46- and 100-way
```
./download_phastCons.sh 46
./download_phastCons.sh 100
```



## EIGEN
```
https://xioniti01.u.hpc.mssm.edu/v1.1/
```


# Use gwRVIS scores with genomic coordinates, stratified by genomic class
- Location:
```
out/[out_dir]/full_genome_out/BED

# full_genome.intergenic.bed
# mid_window_point.intergenic.bed
# etc.
```
