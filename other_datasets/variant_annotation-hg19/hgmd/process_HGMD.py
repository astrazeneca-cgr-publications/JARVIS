import sys
import re

hgmd_vcf = 'hgmd_pro_2019.2_hg19.vcf' #sys.argv[1] 

#1       865595  CM1613956       A       G       .       .       CLASS=DM?;MUT=ALT;GENE=SAMD11;STRAND=+;DNA=NM_152486.2:c.133A>G;PROT=NP_689699.2:p.K45E;DB=rs903331232;PHEN="Retinitis_pigmentosa";RANKSCORE=0.21


pathogenic_out = open('hgmd.pathogenic.bed', 'w')
pathogenic_no_chr_prefix_out = open('hgmd.pathogenic.no_chr_prefix.bed', 'w')

cnt = 0
with open(hgmd_vcf, 'r') as fh:
	for line in fh:
		if line.startswith('#'):
			continue

		line = line.rstrip()

		vals = line.split('\t')
		chrom, position = vals[0], int(vals[1])

		info_vals = vals[7].split(';')
		class_field = info_vals[0]
		_, variant_class = class_field.split('=')
	

		# convert to 0-based BED format
		coords = 'chr' + chrom + '\t' + str(position - 1) + '\t' + str(position)
		no_chr_prefix_coords = chrom + '\t' + str(position - 1) + '\t' + str(position)

		#if re.match('pathogenic', cln_sig, re.IGNORECASE):
		if variant_class == 'DM':
			pathogenic_out.write(coords + '\tPathogenic\n')
			pathogenic_no_chr_prefix_out.write(no_chr_prefix_coords + '\tPathogenic\n')
		

		
