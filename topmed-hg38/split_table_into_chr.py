import sys
import os

table_file = sys.argv[1]
filter_arg = bool(int(sys.argv[2]))
snv_only_arg = bool(int(sys.argv[3]))

print('Table file:', table_file)
print('PASS only variants:', filter_arg)
print('SNVs only:', snv_only_arg)

filter_annot = ''
if filter_arg:
	filter_annot = '-FILTERED'
if snv_only_arg:
	filter_annot = '-SNV_only' + filter_annot

out_dir = 'filtered_variant_tables' + filter_annot
if not os.path.exists(out_dir):
	os.makedirs(out_dir)


# Create file handles for each chromosome and keep them in a list
header = "POS\tREF\tALT\tQUAL\tAC\tAF\tAN\tDP\n"
chr_fh_dict = dict()
for chr in range(1,23):
	fh = open(out_dir + '/chr' + str(chr) + '_topmed-hg38_table.all.txt.filtered', 'w')
	chr_fh_dict['chr' + str(chr)] = fh
	fh.write(header)

for chr in ['X', 'Y']:
	fh = open(out_dir + '/chr' + chr + '_topmed-hg38_table.txt.filtered', 'w')
	chr_fh_dict['chr' + chr] = fh
	fh.write(header)



cnt = 0
with open(table_file) as lines:
	for l in lines:
		if cnt == 0:
			cnt += 1
			continue

		l = l.rstrip()

		vals = l.split('\t')
		chr = vals[0]
		cur_fh = chr_fh_dict[chr]

		cur_fh.write('\t'.join(vals[1:]) + '\n')

		#['chr1', '10484', 'C', 'T', 'PASS', '1', '7.96381e-06', '125568', '500000']

		#for v in vals[1:-1]:
		#	cur_fh.write(v + "\t")
		#af_vals = vals[-1].split('=')
		#cur_fh.write(af_vals[1] + "\n")


		cnt += 1
		if cnt % 10000000 == 0:
			print(chr)

# close file handles
for chr, fh in chr_fh_dict.items():
	fh.close()
