#!/bin/env python3

########################
############ Table filtering
########################

# Librairies
import sys
import os.path
import pandas as pd
import matplotlib.pyplot as plt
import subprocess

# Architecture
dir = '/shared/projects/gametic_segregation_distortion'
data_root = dir + '/' + 'data'
data_file = data_root + '/' + 'genealogies_table_TRUE.txt'
root = dir + '/' + 'results'

pathways_dic = {
	'table' : root + '/' + 'table'
}

external_prg = {
	'script_cM' : dir + '/' + 'ext_script/linear_interpolation_onto_the_map.pl'
}

cross = str(sys.argv[1])
print(cross)

# fonction
def genealogy_table_to_dict(input_file):
	'''
	Read the data table and make it into a dictionnary
	'''
	major_dic = {}	# Output
	
	with open(input_file, mode = 'r') as genealogy_file:	# File created before the run
		lines = genealogy_file.readlines() 
		
		for line in lines:
			line = line.strip('\n')
			my_array = line.split('\t')
			
			if len(my_array) == 4:
				
				cross = my_array[0]
				indiv_type = my_array[1]
				indiv_name = my_array[2]
				raw_fastq_path = my_array[3]
				
				if cross not in major_dic.keys():
					major_dic[cross] = { (indiv_type, indiv_name) : [raw_fastq_path]}
				
				elif (indiv_type, indiv_name) not in major_dic[cross].keys():
					major_dic[cross][(indiv_type, indiv_name)] = [raw_fastq_path]
				
				else:
					major_dic[cross][(indiv_type, indiv_name)].append(raw_fastq_path)

			else:
				continue

	return major_dic



# Setup
data_dict = genealogy_table_to_dict(data_file)

input_table = pathways_dic['table'] + '/' + cross  + '.table'
output_table_no100bp = pathways_dic['table'] + '/' + cross + '.table.no100bp.out'
output_table = pathways_dic['table'] + '/' + cross + '.table.out'

individuals = [ x[1] for x in data_dict[cross].keys() ]

ind_parent1 = [ x[1] for x in data_dict[cross].keys() if x[0] == 'P01' ][0]
ind_parent2 = [ x[1] for x in data_dict[cross].keys() if x[0] == 'P02' ][0]
ind_pollen = [ x[1] for x in data_dict[cross].keys() if x[0] == 'p' ][0]
ind_leaves = [ x[1] for x in data_dict[cross].keys() if x[0] == 'L' ][0]

ind_parent1_AD = ind_parent1 + '.AD'
ind_parent2_AD = ind_parent2 + '.AD'
ind_pollen_AD = ind_pollen + '.AD'
ind_leaves_AD = ind_leaves + '.AD'

ind_parent1_GT = ind_parent1 + '.GT'
ind_parent2_GT = ind_parent2 + '.GT'
ind_pollen_GT = ind_pollen + '.GT'
ind_leaves_GT = ind_leaves + '.GT'

ind_parent1_DP = ind_parent1 + '.DP'
ind_parent2_DP = ind_parent2 + '.DP'
ind_pollen_DP = ind_pollen + '.DP'
ind_leaves_DP = ind_leaves + '.DP'

ind_parent1_GQ = ind_parent1 + '.GQ'
ind_parent2_GQ = ind_parent2 + '.GQ'
ind_pollen_GQ = ind_pollen + '.GQ'
ind_leaves_GQ = ind_leaves + '.GQ'

data = pd.read_table(input_table)

list_ind_DP = [ind_parent1_DP, ind_parent2_DP, ind_pollen_DP, ind_leaves_DP]
list_ind_GQ = [ind_parent1_GQ, ind_parent2_GQ, ind_pollen_GQ, ind_leaves_GQ]
list_ind_GT = [ind_parent1_GT, ind_parent2_GT, ind_pollen_GT, ind_leaves_GT]
list_ind_AD = [ind_parent1_AD, ind_parent2_AD, ind_pollen_AD, ind_leaves_AD]

# Counter
cnt = []
cnt.append(data.shape[0])

print('Starting...')

percentile = []

# Formatting data
for indiv in list_ind_GT:
	data[indiv] = data[indiv].str.replace('|', '/', regex = False).str.split('/')   # Format the GT field of each individual

for indiv in list_ind_AD:
	data[indiv] = data[indiv].str.split(',')

print('Formatted...')

# Keep only scaffold 1 to 8
data = data[ data['CHROM'].str.match('^scaffold_[1-8]$') ]

cnt.append(data.shape[0])

print('Scaffold removed...')

# Filter on SNP
data = data[ data['TYPE'] == 'SNP' ]

cnt.append(data.shape[0])

print('SNP only kept...')

# Filter on GQ
for indiv in list_ind_GQ:
	data = data[data[indiv] >= 30]

cnt.append(data.shape[0])

print('GQ filtered...')

# Filter on DP
sum_DP = data.loc[:, list_ind_DP].sum(axis = 1)

inf = sum_DP.quantile(0.1)
sup = sum_DP.quantile(0.9)

plt.hist(sum_DP, bins = 'auto')
plt.savefig(''.join([cross, '.rawDP.pdf']))
plt.close()

inf_bool = sum_DP > inf
sup_bool = sum_DP < sup

data2 = data[sup_bool]
sum_DP2 = data2.loc[:, list_ind_DP].sum(axis = 1)
plt.hist(sum_DP2, bins = 'auto')
plt.savefig(''.join([cross, '.filteredDPsup.pdf']))
plt.close()

sub_data = data[inf_bool & sup_bool]
sum_DP3 = sub_data.loc[:, list_ind_DP].sum(axis = 1)
plt.hist(sum_DP3, bins = 'auto')
plt.savefig(''.join([cross, '.filteredDPboth.pdf']))
plt.close()

percentile.append(inf)
percentile.append(sup)
cnt.append(sub_data.shape[0])

print('DP filtered...')

# Check that each variant position is consistent with mendelian expectations (somatic only)
P01_hom = sub_data[ind_parent1_GT].str[0] == sub_data[ind_parent1_GT].str[1]	# Parent 1 homozygote
P02_hom = sub_data[ind_parent2_GT].str[0] == sub_data[ind_parent2_GT].str[1]	# Parent 2 homozygote
Parents_het = sub_data[ind_parent1_GT] != sub_data[ind_parent2_GT]					# Parents heterozygote
L_het = sub_data[ind_leaves_GT].str[0] != sub_data[ind_leaves_GT].str[1]		# Leaf sample is heterozygous

similar_leaves_parents = sub_data.apply( lambda x: \
		(x[ind_parent1_GT][0] in set(x[ind_leaves_GT]) \
		and x[ind_parent2_GT][0] in set(x[ind_leaves_GT])), axis = 1 )  # Variants present in leaf sample are also present in both parents     

all_conditions = P01_hom & P02_hom & Parents_het & L_het & similar_leaves_parents 

osub_data = sub_data[ all_conditions ] # Applying the conditions to filter

df = pd.DataFrame()
df['CHROM'] = osub_data['CHROM']
df['POS'] = osub_data['POS']
df['somatic_parental_library_1'] = osub_data.apply(lambda x: x[ind_leaves_AD][0] if x[ind_leaves_GT][0] in x[ind_parent1_GT] else x[ind_leaves_AD][1], axis = 1)
df['somatic_parental_library_2'] = osub_data.apply(lambda x: x[ind_leaves_AD][0] if x[ind_leaves_GT][0] in x[ind_parent2_GT] else x[ind_leaves_AD][1], axis = 1)
df['germline_parental_library_1'] = osub_data.apply(lambda x: x[ind_pollen_AD][0] if x[ind_pollen_GT][0] in x[ind_parent1_GT] else x[ind_pollen_AD][1], axis = 1)
df['germline_parental_library_2'] = osub_data.apply(lambda x: x[ind_pollen_AD][0] if x[ind_pollen_GT][0] in x[ind_parent2_GT] else x[ind_pollen_AD][1], axis = 1)
df.to_csv(output_table_no100bp, header=None, index=None, sep='\t', mode='w')

print('Non-mendelian filtered...')

cnt.append(osub_data.shape[0])

# Less than 100bp from each other 
duplicated = list( ~osub_data.duplicated( subset=['CHROM'] ) )
position = list( osub_data['POS'] )

lag = [int] * ( len(osub_data) + 1 )
lag[0:150] = [ position[x] - position[0] for x in range(150) if x < len(position) ]

my_vec = []

for pos in range(len(osub_data)):
	if duplicated[pos]:
		my_vec.append(True)
		lag[pos:(pos + 150)] = [ position[x] - position[pos] for x in range(pos, (pos + 150)) if x < len(position) ]
	else:
		if lag[pos] > 150:
			lag[pos:(pos + 150)] = [ position[x] - position[pos] for x in range(pos, (pos + 150)) if x < len(position) ]
			my_vec.append(True)
		else:
			my_vec.append(False)

final_data = osub_data[my_vec]

cnt.append(final_data.shape[0])

print(cross, ":", cnt)

df = pd.DataFrame()
df['CHROM'] = final_data['CHROM']
df['POS'] = final_data['POS']
df['somatic_parental_library_1'] = final_data.apply(lambda x: x[ind_leaves_AD][0] if x[ind_leaves_GT][0] in x[ind_parent1_GT] else x[ind_leaves_AD][1], axis = 1)
df['somatic_parental_library_2'] = final_data.apply(lambda x: x[ind_leaves_AD][0] if x[ind_leaves_GT][0] in x[ind_parent2_GT] else x[ind_leaves_AD][1], axis = 1)
df['germline_parental_library_1'] = final_data.apply(lambda x: x[ind_pollen_AD][0] if x[ind_pollen_GT][0] in x[ind_parent1_GT] else x[ind_pollen_AD][1], axis = 1)
df['germline_parental_library_2'] = final_data.apply(lambda x: x[ind_pollen_AD][0] if x[ind_pollen_GT][0] in x[ind_parent2_GT] else x[ind_pollen_AD][1], axis = 1)
df.to_csv(output_table, header=None, index=None, sep='\t', mode='w')

print('Less than 100bp from each other filtered...')

# Preparing data for plotting

my_output = output_table + '.toPlot.txt'

data_table = pd.read_table(output_table, header = None)

dic_output = {}

dic_output['scaffold'] = []
dic_output['somatic_mean'] = []
dic_output['germline_mean'] = []
dic_output['mean_window'] = []
dic_output['baseline'] = []

for scaffold in data_table.loc[:,0].unique():
	sub_data = data_table[data_table.iloc[:, 0] == scaffold]
	df_somatic = sub_data.iloc[:,3] / ( sub_data.iloc[:,2] + sub_data.iloc[:,3] )
	df_position = sub_data.iloc[:,1]
	somatic_res = []
	scaffold_res = []
	baseline_res = []
	mean_window = []
	for x in range(0, len(df_somatic), 500):
		mn = x
		mx = (x + 499)
		mean_somatic = df_somatic.iloc[mn:mx].mean()
		mean_position = df_position.iloc[mn:mx].mean()
		scaffold_res.append(scaffold)
		somatic_res.append(mean_somatic)
		baseline_res.append(0)
		mean_window.append(mean_position)
	dic_output['scaffold'].extend(scaffold_res)
	dic_output['baseline'].extend(baseline_res)
	dic_output['somatic_mean'].extend(somatic_res)
	dic_output['mean_window'].extend(mean_window)

for scaffold in data_table.loc[:,0].unique():
	sub_data = data_table[data_table.iloc[:, 0] == scaffold]
	df_germline = sub_data.iloc[:,5] / ( sub_data.iloc[:,4] + sub_data.iloc[:,5] )
	germline_res = []
	for x in range(0, len(df_germline), 500):
		mn = x
		mx = (x + 499)
		mean_germline = df_germline.iloc[mn:mx].mean()
		germline_res.append(mean_germline)
	dic_output['germline_mean'].extend(germline_res)

zip_object = zip(dic_output['somatic_mean'], dic_output['germline_mean'])
difference = [ x[0] - x[1] for x in zip_object ]

dic_output['difference_mean'] = difference

df = pd.DataFrame.from_dict(dic_output)

df.to_csv(my_output, index = None, sep = '\t', mode = 'w')

print('plot-ready data...')

# Preparing data for MAPSD

my_output = output_table + '.4MAPSD'

cmd = 'perl' + ' ' + external_prg['script_cM'] + ' ' + output_table + ' > ' + my_output + '\n'

subprocess.call(cmd, shell=True)

print('End !')
