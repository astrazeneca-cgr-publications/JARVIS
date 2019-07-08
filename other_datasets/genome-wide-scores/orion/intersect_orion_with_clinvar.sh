#!/bin/sh

non_cod_include_classes="intergenic intron lincrna mature_mirna ucne utr vista"
coding_include_classes="omim-HI tolerant intolerant ccds"

for cl in $non_cod_include_classes $coding_include_classes; do
	echo $cl

	pathogenic_clinvar_genomic_class_bed="../../out/gnomad-regression_beta-winlen_3000.MAF_0.001.varType_snv.Pop_all-PASS_ONLY-NO_SEGDUP-NO_LCR-high_conf_regions/full_genome_out/BED/${cl}.clinvar_pathogenic.bed"
	benign_clinvar_genomic_class_bed="../../out/gnomad-regression_beta-winlen_3000.MAF_0.001.varType_snv.Pop_all-PASS_ONLY-NO_SEGDUP-NO_LCR-high_conf_regions/full_genome_out/BED/${cl}.clinvar_benign.bed"

	#intersectBed -wo -a orion.ClinVar_benign.bed -b clinvar_genomic_class.bed | cut -f1,2,3,4,9 > ${cl}.orion.Clinvar.bed

	intersectBed -wo -a orion.ClinVar_pathogenic.bed -b $pathogenic_clinvar_genomic_class_bed | cut -f1,2,3,4,9 > ${cl}.orion.Clinvar.pathogenic.bed
	intersectBed -wo -a orion.ClinVar_benign.bed -b $benign_clinvar_genomic_class_bed | cut -f1,2,3,4,9 > ${cl}.orion.Clinvar.benign.bed
	
	echo "python plot_clinvar_densities.py ${cl}.orion.Clinvar.pathogenic.bed ${cl}.orion.Clinvar.benign.bed $cl"
	#python plot_clinvar_densities.py ${cl}.orion.Clinvar.pathogenic.bed ${cl}.orion.Clinvar.benign.bed $cl

	echo "python run_logistic_regression.py ${cl}.orion.Clinvar.pathogenic.bed ${cl}.orion.Clinvar.benign.bed $cl"
	python run_logistic_regression.py ${cl}.orion.Clinvar.pathogenic.bed ${cl}.orion.Clinvar.benign.bed $cl

done
