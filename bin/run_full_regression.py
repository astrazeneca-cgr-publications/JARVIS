import sys
import os
from subprocess import call
from custom_utils import create_out_dir, get_run_parameters

args = sys.argv
config_file = args[1]	#'config.log'

# read run parameters from config file and store into a dictionary
run_params = get_run_parameters(config_file)
print(run_params)

out_dir = create_out_dir(config_file)
print(out_dir)

win_len = run_params['win_len']
all_variants_upper_thres = run_params['all_variants_upper_thres']       
filter_outliers_before_regression = run_params['filter_outliers_before_regression']  
rscript_x11_path = run_params['rscript_x11_path']

# Select R PATH installation on a different system, e.g. a cluster
local_rscript_path = '/usr/bin/Rscript' 
if not os.path.exists(local_rscript_path):         
	rscript_x11_path = 'Rscript'
print("Rscript path: " + rscript_x11_path)

#rvis_dir = out_dir + '/rvis_scores'

# fit linear regression for all autocomal chromosomes
chr_type = 'autosomal'
print("Fitting linear regression for autosomal chromosomes only")
print(rscript_x11_path, 'full_genome_r_studres_glm.R', out_dir, filter_outliers_before_regression, all_variants_upper_thres, win_len, chr_type)
sys.exit()
call([rscript_x11_path, 'full_genome_r_studres_glm.R', out_dir, str(filter_outliers_before_regression), str(all_variants_upper_thres), str(win_len), chr_type])

# fit linear regression for sex chromosomes (only X in that case) separately
chr_type = 'sex'
print("Fitting linear regression for sex chromosomes only (X in this case)")
print(rscript_x11_path, 'full_genome_r_studres_glm.R', out_dir, filter_outliers_before_regression, all_variants_upper_thres, win_len, chr_type)
call([rscript_x11_path, 'full_genome_r_studres_glm.R', out_dir, str(filter_outliers_before_regression), str(all_variants_upper_thres), str(win_len), chr_type])
