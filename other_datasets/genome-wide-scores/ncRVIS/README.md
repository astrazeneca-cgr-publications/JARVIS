# 1. Download genomic boundaries file
https://doi.org/10.1371/journal.pgen.1005492.s015  -->  journal.pgen.1005492.s015

mv journal.pgen.1005492.s015 ncrvis_genomic_boundaries.txt

# 2. Table with ncRVIS and RVIS scores:
ncRVIS_scores_table.csv


# 3. Convert boundaries files to BED file with ncRVIS annotation (merged)
python parse_genomic_boundaries.py	# --> All_chromosomes.ncRVIS.bed


# Run
sbatch ./prepare_variant_specific_files.sh
