import sys
import os
from subprocess import call
from custom_utils import create_out_dir, get_config_params

config_file = sys.argv[1]	#'config.yaml'

# read run parameters from config file and store into a dictionary
config_params = get_config_params(config_file)
print(config_params)

out_dir = create_out_dir(config_file)
print(out_dir)

win_len = config_params['win_len']
all_variants_upper_thres = config_params['all_variants_upper_thres']       
#filter_outliers_before_regression = config_params['filter_outliers_before_regression']  

# Select R PATH installation on a different system, e.g. a cluster
#local_rscript_path = '/usr/bin/Rscript' 
#if not os.path.exists(local_rscript_path):         
#	rscript_x11_path = 'Rscript'
#print("Rscript path: " + rscript_x11_path)


# fit linear regression for all autocomal chromosomes
chr_type = 'autosomal'
print("Fitting linear regression for autosomal chromosomes only")
print('Rscript', 'full_genome_r_studres_glm.R', out_dir, all_variants_upper_thres, win_len, chr_type)
call(['Rscript', 'full_genome_r_studres_glm.R', out_dir, str(all_variants_upper_thres), str(win_len), chr_type])

# fit linear regression for sex chromosomes (only X in that case) separately
#chr_type = 'sex'
#print("Fitting linear regression for sex chromosomes only (X in this case)")
#print(rscript_x11_path, 'full_genome_r_studres_glm.R', out_dir, filter_outliers_before_regression, all_variants_upper_thres, win_len, chr_type)
#call([rscript_x11_path, 'full_genome_r_studres_glm.R', out_dir, str(filter_outliers_before_regression), str(all_variants_upper_thres), str(win_len), chr_type])
