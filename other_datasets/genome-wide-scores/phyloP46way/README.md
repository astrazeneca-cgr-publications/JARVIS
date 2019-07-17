tabix All_chromosomes.phyloP46way.bed.bgz -B ../../clinvar/clinvar.pathogenic.bed > phylop46way.clinvar_pathogenic.bed &
tabix All_chromosomes.phyloP46way.bed.bgz -B ../../clinvar/clinvar.benign.bed > phylop46way.clinvar_benign.bed &

cat phyloP46way.clinvar_benign.bed phyloP46way.clinvar_pathogenic.bed > phyloP46way.clinvar.bed
