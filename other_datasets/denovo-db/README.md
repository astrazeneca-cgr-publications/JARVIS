# Process original de-novo tsv files into BED files with selected columns

```
zcat < denovo-db.non-ssc-samples.variants.tsv.gz | grep -v '#' | awk '{ print "chr"$9"\t"$10-1"\t"$10"\t"$7"\t"$22 "\t"$3"\t"$4"\t"$5"\t"$6"\t"$19"\t"$20"\t"$21"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31 }' > bed/denovo-db.non-ssc-samples.variants.bed;

zcat < denovo-db.ssc-samples.variants.tsv.gz | grep -v '#' | awk '{ print "chr"$9"\t"$10-1"\t"$10"\t"$7"\t"$22 "\t"$3"\t"$4"\t"$5"\t"$6"\t"$19"\t"$20"\t"$21"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31 }' > bed/denovo-db.ssc-samples.variants.bed
```

> Get header of output BED files:
```
zcat < denovo-db.non-ssc-samples.variants.tsv.gz | grep -v '##' | grep '#' | awk '{ print "chr"$9"\t"$10-1"\t"$10"\t"$7"\t"$22 "\t"$3"\t"$4"\t"$5"\t"$6"\t"$19"\t"$20"\t"$21"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31 }'
# Result:
#chrChr  -1      Position        PrimaryPhenotype        FunctionClass   PubmedID	NumProbands	NumControls	SequenceType	Transcript	codingDnaSize	Gene	PolyPhen(HDiv)	PolyPhen(HVar)	SiftScore	CaddScore	LofScore	LrtScore

# (First three columns are: chrChr  -1      Position  --> chr	star	end)
```

