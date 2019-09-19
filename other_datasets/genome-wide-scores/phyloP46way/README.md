# > ClinVar
tabix All_chromosomes.phyloP46way.bed.bgz -B ../../clinvar/clinvar.pathogenic.bed > phyloP46way.clinvar_pathogenic.bed &
tabix All_chromosomes.phyloP46way.bed.bgz -B ../../clinvar/clinvar.benign.bed > phyloP46way.clinvar_benign.bed &

cat phyloP46way.clinvar_benign.bed phyloP46way.clinvar_pathogenic.bed > phyloP46way.clinvar.bed


# > HGMD
tabix All_chromosomes.phyloP46way.bed.bgz -B ../../hgmd/hgmd.pathogenic.bed > phyloP46way.hgmd_pathogenic.bed &

# Combine HGMD pathogenic with ClinVar benign (or other versions of benign later on):
cat phyloP46way.hgmd_pathogenic.bed phyloP46way.clinvar_benign.bed > phyloP46way.hgmd.bed
