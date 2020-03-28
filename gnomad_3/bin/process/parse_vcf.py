import sys
import gzip
import re
import numpy as np



def parse_info_fields(info_pairs_dict):
	""" 
	Fields to retain: 
			AC 
			AN 
			DP
	"""

	other_fields_str = ''

	fields_to_retain = ['AC', 'AF', 'AN', 'DP']

	# Also reconstruct 'AF' from 'AC' and 'AN' (absent in gnomad v3.0 ?):
	# AF = AC / AN
	#try:
	#	info_pairs_dict['AF'] = float(info_pairs_dict['AC']) / float(info_pairs_dict['AN'])
	#except:
	#	info_pairs_dict['AF'] = 0

	for field in fields_to_retain:

		value = info_pairs_dict[field]
		other_fields_str += '\t' + str(value)
		print(field, value)

	#if verbose:
	#	print('other_fields_str:\n' + other_fields_str)


	return other_fields_str




def filter_based_on_flags(FILTER, INFO, KEEP_PASS_ONLY, FILTER_LCR):
	"""
	    Fitler out entries based on FILTER annotation (e.g. PASS, RF, etc.)
	    and any 'lcr' (low complexity region) flag in the INFO field
	"""

	valid = True
	if KEEP_PASS_ONLY:
		if FILTER != 'PASS':
			valid = False

	if FILTER_LCR:
		if 'lcr' in INFO:
			valid = False
	
	return valid
			



def process_variant_entry(line, out_fh):

	CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = line.split('\t')

	valid_entry = filter_based_on_flags(FILTER, INFO, KEEP_PASS_ONLY, FILTER_LCR)
	if not valid_entry:
		#print("Not valid entry")
		return valid_entry

	#print(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)

	base_columns = [POS, REF, ALT, QUAL]
	base_output_str = '\t'.join(base_columns)

	if verbose:
		print('\n\n=============================================')
		print('base_output_str:', base_output_str)

	

	# - Parse paired "field=value" fields
	info_pairs_dict = dict([ f.split('=') for f in INFO.split(';') if '=' in f])
	print(info_pairs_dict)


	info_str = ''
	# > Get all relevant fields of interest
	info_str += parse_info_fields(info_pairs_dict)
	if verbose:
		print('info_str:', info_str)



	full_out_str = base_output_str + info_str + '\n'
	if verbose:
		print(full_out_str)


	out_fh.write(full_out_str)
	
	return valid_entry





if __name__ == '__main__':

	# ======================
	dbg=True
	verbose=True 
	# ======================

	dataset = sys.argv[1]
	chrom = sys.argv[2]
	input_file = sys.argv[3]
	KEEP_PASS_ONLY = bool(int(sys.argv[4]))
	FILTER_LCR = bool(int(sys.argv[5]))
	out_dir = sys.argv[6]

	print('Input file:', input_file)

	out_file = out_dir + '/chr' + chrom + '_' + dataset + '_table.all.txt.filtered' 
	print('Output file:', out_file)

	out_fh = open(out_file, 'w')
	out_fh.write("POS\tREF\tALT\tQUAL\tAC\tAF\tAN\tDP\n")


	cnt = 0
	with gzip.open(input_file) as fh:
		for line in fh:
			line = line.decode('utf-8')
			line = line.rstrip()

			if line.startswith('#'):
				continue	
			else:
				valid_entry = process_variant_entry(line, out_fh)

				if dbg:
					if valid_entry: 
						cnt += 1
						if verbose: print(line)

			if dbg:
				if cnt > 5: sys.exit()

	out_fh.close()
