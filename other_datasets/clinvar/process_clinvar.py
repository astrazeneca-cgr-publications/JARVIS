import sys, os
import gzip
import re

clinvar_vcf = 'clinvar.vcf.gz' #sys.argv[1] 

#1       949422  475283  G       A       .       .       AF_ESP=0.00546;AF_EXAC=0.00165;AF_TGP=0.00619;ALLELEID=446939;CLNDISDB=MedGen:C4015293,OMIM:616126,Orphanet:ORPHA319563;CLNDN=Immunodeficiency_38_with_basal_ganglia_calcification;CLNHGVS=NC_000001.10:g.949422G>A;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Benign;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=ISG15:9636;MC=SO:0001583|missense_variant;ORIGIN=1;RS=143888043

pathogenic_out = open('clinvar.pathogenic.bed', 'w')
benign_out = open('clinvar.benign.bed', 'w')
pathogenic_no_chr_prefix_out = open('clinvar.pathogenic.no_chr_prefix.bed', 'w')
benign_no_chr_prefix_out = open('clinvar.benign.no_chr_prefix.bed', 'w')

clinsigs = dict()

with gzip.open(clinvar_vcf, 'r') as fh:
	for line in fh:
		line = line.decode('utf-8')
		if line.startswith('#'):
			continue

		line = line.rstrip()

		vals = line.split('\t')
		chrom, position = vals[0], int(vals[1])

		info_vals = vals[7].split(';')
		info_dict = dict(item.split("=") for item in info_vals)

		if 'CLNSIG' not in info_dict:
			continue

		cln_sig = info_dict['CLNSIG']

		# convert to 0-based BED format
		coords = 'chr' + chrom + '\t' + str(position - 1) + '\t' + str(position)
		no_chr_prefix_coords = chrom + '\t' + str(position - 1) + '\t' + str(position)

		#if re.match('pathogenic', cln_sig, re.IGNORECASE):
		if 'Pathogenic' in cln_sig:
			pathogenic_out.write(coords + '\t' + cln_sig + '\n')
			pathogenic_no_chr_prefix_out.write(no_chr_prefix_coords + '\t' + cln_sig + '\n')
		#elif re.match('benign', cln_sig, re.IGNORECASE):
		elif 'Benign' in cln_sig:
			benign_out.write(coords + '\t' + cln_sig + '\n')
			benign_no_chr_prefix_out.write(no_chr_prefix_coords + '\t' + cln_sig + '\n')
		
pathogenic_out.close()
benign_out.close()
pathogenic_no_chr_prefix_out.close()
benign_no_chr_prefix_out.close()

		
os.system("cat clinvar.pathogenic.bed | mergeBed | sed 's/$/\tPathogenic/g' > clinvar.pathogenic.merged.bed ")
os.system("cat clinvar.benign.bed | mergeBed | sed 's/$/\tBenign/g' > clinvar.benign.merged.bed ")

os.system("cat clinvar.pathogenic.no_chr_prefix.bed | mergeBed | sed 's/$/\tPathogenic/g' > clinvar.pathogenic.no_chr_prefix.merged.bed ")
os.system("cat clinvar.benign.no_chr_prefix.bed | mergeBed | sed 's/$/\tBenign/g' > clinvar.benign.no_chr_prefix.merged.bed ")
