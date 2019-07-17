tabix All_chromosomes.phastCons46way.bed.bgz -B ../../clinvar/clinvar.pathogenic.bed > phastCons46way.clinvar_pathogenic.bed &
tabix All_chromosomes.phastCons46way.bed.bgz -B ../../clinvar/clinvar.benign.bed > phastCons46way.clinvar_benign.bed &

cat phastCons46way.clinvar_pathogenic.bed phastCons46way.clinvar_benign.bed > phastCons46way.clinvar.bed
