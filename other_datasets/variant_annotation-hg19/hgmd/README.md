# 1. Parse HGMD VCF into BED file
```
python process_HGMD.py
```

# 2. Get file with unique HGMD variants
```
cat hgmd.pathogenic.bed | awk '!seen[$1"_"$2]  {print $0} {++seen[$1"_"$2]}' > hgmd.pathogenic.uniq.bed

cat hgmd.pathogenic.no_chr_prefix.bed | awk '!seen[$1"_"$2]  {print $0} {++seen[$1"_"$2]}' > hgmd.pathogenic.uniq.no_chr_prefix.bed
```
