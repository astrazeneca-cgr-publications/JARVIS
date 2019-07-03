# Sort BED files of Ensembl regulatory elements
```
for f in Monocytes_CD14plus.*; do cat "$f" | sort -k1,1 -k2,2n > "${f/bed/sorted.bed}"; done
```
