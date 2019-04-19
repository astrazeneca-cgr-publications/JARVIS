from sys import argv, exit
from subprocess import call
import gzip


chrom = argv[1]
population = argv[2]
in_vcf_dir = argv[3]
out_vcf_dir = argv[4]

#21      9411211 .       A       G       43.66   RF      AC=1;AF=4.81047e-05;AN=20788;BaseQRankSum=3.69000e-01;ClippingRankSum=-1.73900e+00;DP=263129;FS=0.00000e+00;InbreedingCoeff=-6.50000e-03;MQ=5.04000e+01;MQRankSum=5.30000e-02;QD=2.57000e+00;ReadPosRankSum=7.91000e-01;SOR=6.93000e-01;VQSLOD=-8.17500e+00;VQSR_culprit=MQ;GQ_HIST_ALT=0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|1;DP_HIST_ALT=0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;AB_HIST_ALT=0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0;GQ_HIST_ALL=198|517|492|1261|1763|948|2000|1933|908|1568|1224|507|1089|109|322|121|242|21|98|102;DP_HIST_ALL=1023|3863|4861|3417|1619|415|187|32|5|1|0|0|0|0|0|0|0|0|0|0;AB_HIST_ALL=0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0;AC_Male=0;AC_Female=1;AN_Male=11692;AN_Female=9096;AF_Male=0.00000e+00;AF_Female=1.09938e-04;GC_Male=5846,0,0;GC_Female=4547,1,0;GC_raw=15422,1,0;AC_raw=1;AN_raw=30846;GC=10393,1,0;AF_raw=3.24191e-05;Hom_AFR=0;Hom_AMR=0;Hom_ASJ=0;Hom_EAS=0;Hom_FIN=0;Hom_NFE=0;Hom_OTH=0;Hom=0;Hom_raw=0;AC_AFR=0;AC_AMR=0;AC_ASJ=0;AC_EAS=0;AC_FIN=0;AC_NFE=1;AC_OTH=0;AN_AFR=6222;AN_AMR=572;AN_ASJ=206;AN_EAS=1064;AN_FIN=2468;AN_NFE=9668;AN_OTH=588;AF_AFR=0.00000e+00;AF_AMR=0.00000e+00;AF_ASJ=0.00000e+00;AF_EAS=0.00000e+00;AF_FIN=0.00000e+00;AF_NFE=1.03434e-04;AF_OTH=0.00000e+00;POPMAX=NFE;AC_POPMAX=1;AN_POPMAX=9668;AF_POPMAX=1.03434e-04;DP_MEDIAN=17;DREF_MEDIAN=2.51189e-11;GQ_MEDIAN=99;AB_MEDIAN=2.94118e-01;AS_RF=1.02455e-01;AS_FilterStatus=RF;CSQ=G|intergenic_variant|MODIFIER||||||||||||||||1||||SNV|1||||||||||||||||||||||||||||||||||||||||||||;GC_AFR=3111,0,0;GC_AMR=286,0,0;GC_ASJ=103,0,0;GC_EAS=532,0,0;GC_FIN=1234,0,0;GC_NFE=4833,1,0;GC_OTH=294,0,0;Hom_Male=0;Hom_Female=0;segdup

input_file = in_vcf_dir + '/gnomad.genomes.r2.0.2.sites.chr' + chrom + '.vcf.bgz'
out_vcf_file = out_vcf_dir + '/gnomad.genomes.r2.0.2.sites.chr' + chrom + '.vcf'
target_substr = 'AF_' + population

out_fh = open(out_vcf_file, 'w')

#cnt=0
with gzip.open(input_file) as gz:
	for line in gz:
		line = line.decode("utf-8")
		if line.startswith('#'):
			out_fh.write(line)
		else:
			vals = line.split('\t')
	
			info_fields = vals[7].split(';')
			target_field = [s for s in info_fields if target_substr in s]
			target_field = target_field[0]

			alt_alleles = vals[4].split(',')

			alt_alleles_idxs_to_include = []
			alt_afs_to_include = []

			af_subpop_field = target_field.split('=')
			af_vals = af_subpop_field[1].split(',')
			for i, v in enumerate(af_vals):
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
				new_info_field = "POPULATION=" + population + ";AF=" + passed_afs_field + ';' + vals[7]  # replaced total AF values with AF_{SUBPOPULATION} ones
				new_line += "\t" + new_info_field
				out_fh.write(new_line + "\n")


				"""
				print("--------------------------------\n")
				print(line)

				print(alt_alleles_idxs_to_include)
				print(alt_alleles_passed)
				print(passed_alleles_field)
				print(passed_afs_field)

				print(new_line)
				print("\n\n")
				
				cnt += 1
				"""
		

			#if cnt > 5:
			#	out_fh.close()
			#	exit()

out_fh.close()

