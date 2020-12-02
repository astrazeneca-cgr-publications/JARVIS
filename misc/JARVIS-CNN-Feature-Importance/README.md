1. Explore_JARVIS_CNN_Feature_Importance.ipynb

2. cdhit:
```
cd-hit-est -i JARVIS-CNN_pattern_seqs.fasta -o cdhit_seq_clusters.id80.W7.fa -c 0.8 -n 7
```

3. process_cdhit_clusters.py


4. Tomtom run on web-server:
- Target motifs: 1,808
```
tomtom -no-ssc -oc . -verbosity 1 -min-overlap 5 -dist pearson -evalue -thresh 10.0 query_motifs db/EUKARYOTE/jolma2013.meme db/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme db/MOUSE/uniprobe_mouse.meme
```

5. Download results:
- tomtom.tsv
- tomtom.html

6. Process TomTom output (tomtom.tsv):
Process_TomTom_output.ipynb
