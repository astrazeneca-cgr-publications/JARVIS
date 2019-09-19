### QC for intervals distribution

#### CCDS
v=19
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
zcat gencode.v$v.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="CDS" {print $1,$4-1,$5}' | sortBed | mergeBed -i - > gencode_v${v}_cds_merged.bed


#### Get total region length from BED files:
```
cat gencode_v19_intron.bed | awk '{sum+=$3-$2} END {print sum}';
cat gencode_v19_intergenic.bed | awk '{sum+=$3-$2} END {print sum}';
```


#### Get intervals from BED files:
```
 cat gencode_v19_intron.bed | awk '{ print ($3 - $2) }' > gencode_v19_intron.bed.intervals
 cat gencode_v19_intergenic.bed | awk '{ print ($3 - $2) }' > gencode_v19_intergenic.bed.intervals
```

#### intergenic
```
>>> s =  pd.read_table('gencode_v19_intergenic.bed.intervals', header=None)
>>> s.columns = ['dist'] 
>>> np.median(s.dist) 
13388.0 
>>> np.std(s.dist)
335242.718956013 
>>> np.mean(s.dist)
46027.50019054319 
```

#### introns
```
>>> s =  pd.read_table('gencode_v19_intron.bed.intervals', header=None) 
s.columns = ['dist'] 
>>> np.median(s.dist) 
1416.0 
>>> np.std(s.dist) 
13295.886782352303 
>>> np.mean(s.dist) 
5025.9008046141835
```

### Merge VISTA and PHANTOM5 enhancers
```
cat vista_enhancers_genes_list.bed | cut -f1,2,3 > vista_enhancers.bed;  

cat vista_enhancers.bed phantom5_human_permissive_enhancers_phase_1_and_2.bed > merged_vista_and_phantom.bed;   

sort -k1,1 -k2,2n merged_vista_and_phantom.bed > merged_vista_and_phantom.sorted.bed;
```
