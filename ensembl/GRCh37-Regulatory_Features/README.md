### Sort BED files of Ensembl regulatory elements
```
for f in Monocytes_CD14plus.*; do echo "$f"; cat "$f" | sort -k1,1 -k2,2n | sed 's/^/chr/g' | sed 's/ /_/g' | sed 's/\tHistone//g' > "${f/bed/sorted.bed}"; done
```
