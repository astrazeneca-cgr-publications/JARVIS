### QC for intervals distribution

# Get genomic annotation files:
```
./get_annotations.sh
```
# Based on: https://davetang.org/muse/2013/01/18/defining-genomic-regions/


#### Get total region length from BED files:
```
cat gencode_v${v}_intron.bed | awk '{sum+=$3-$2} END {print sum}';
cat gencode_v${v}_intergenic.bed | awk '{sum+=$3-$2} END {print sum}';
```


#### Get intervals from BED files:
```
 cat gencode_v${v}_intron.bed | awk '{ print ($3 - $2) }' > gencode_v${v}_intron.bed.intervals
 cat gencode_v${v}_intergenic.bed | awk '{ print ($3 - $2) }' > gencode_v${v}_intergenic.bed.intervals
```

#### intergenic
```
>>> s =  pd.read_table('gencode_v${v}_intergenic.bed.intervals', header=None)
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
>>> s =  pd.read_table('gencode_v${v}_intron.bed.intervals', header=None) 
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
