from sys import argv, exit
from subprocess import call
import gzip
import re


chrom = argv[1]
population = argv[2]
in_vcf_dir = argv[3]
out_vcf_dir = argv[4]


input_file = in_vcf_dir + '/gnomad.genomes.r2.1.1.sites.' + chrom + '.vcf.bgz'
out_vcf_file = out_vcf_dir + '/gnomad.genomes.r2.1.1.sites.' + chrom + '.vcf.bgz'
target_substr = 'AF_' + population.lower() + '=' # gnomad r2.1.1 has change population identifiers with lowercase letters

out_fh = gzip.open(out_vcf_file, 'w')

cnt=0
with gzip.open(input_file) as gz:
	for line in gz:
		line = line.decode("utf-8")
		if line.startswith('#'):
			out_fh.write(line.encode())
		else:
			vals = line.split('\t')
	
			info_fields = vals[7].split(';')
			target_field = [s for s in info_fields if re.match(target_substr, s)]
			try:
				target_field = target_field[0]
			except:
				continue

			alt_alleles = vals[4].split(',')

			alt_alleles_idxs_to_include = []
			alt_afs_to_include = []

			af_subpop_field = target_field.split('=')
			af_vals = af_subpop_field[1].split(',')
			for i, v in enumerate(af_vals):
				# filter out entries with zero Allele Frequency
				if (v != '.') and (float(v) > 0.0):
					alt_alleles_idxs_to_include.append(i)
					alt_afs_to_include.append(v)
			

			alt_alleles_passed = [alt_alleles[i] for i in alt_alleles_idxs_to_include]
			passed_alleles_field = ",".join(alt_alleles_passed)
			passed_afs_field = ",".join(alt_afs_to_include)

			if len(alt_alleles_passed) >= 1: # qualifying viarants in specified subpopulation
				new_line = "\t".join(vals[:4])
				new_line += "\t" + passed_alleles_field
				new_line += "\t" + vals[5]
				new_line += "\t" + vals[6]

				info_vals = vals[7].split(';')
				info_vals[2] = "AF=" + passed_afs_field

				dp_index = [i for i in range(len(info_vals)) if 'DP=' in info_vals[i]]

				new_info_field = ';'.join(info_vals[:3]) + ';' + info_vals[dp_index[0]] # replaced total AF values with AF_{SUBPOPULATION} ones
				new_line += "\t" + new_info_field
				out_fh.write((new_line + "\n").encode())

				#print("--------------------------------\n")
				#print(alt_alleles_idxs_to_include)
				#print(alt_alleles_passed)
				#print(passed_alleles_field)
				#print(passed_afs_field)

				#print(new_line)
				#print("\n\n")
				#cnt += 1

			#if cnt > 5:
			#	out_fh.close()
			#	exit()

out_fh.close()
