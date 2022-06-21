#!/bin/env python3

########################
############ Table analysis
########################

# Librairies
import sys
import subprocess
import pandas as pd

# Architecture
dir = '/shared/projects/gametic_segregation_distortion'
root = dir + '/' + 'results'

pathways_dic = {
	'table' : root + '/' + 'table'
}

external_prg = {
	'script_MAPSD' : dir + '/' + 'ext_script/MAP_SD-master/src/map_sd',
	'script_R_before_MAPSD' : dir + '/' + 'ext_script/plot_beforeMAPSD.R',
	'script_R_after_MAPSD' : dir + '/' + 'ext_script/plot_afterMAPSD.R'
}


cross = str(sys.argv[1])
print(cross)

# Setup

df_MAPSD =  pathways_dic['table'] + '/' + cross + '.table.out.4MAPSD'

df_out_main = pathways_dic['table'] + "/" + cross + "_MAPSD.main.txt"
df_out_bootstrap = pathways_dic['table'] + "/" + cross + "_MAPSD.bootstrap.txt"
df_out_windows = pathways_dic['table'] + "/" + cross + "_MAPSD.windows.txt"

# Creation of plots of the difference between germinal and somatic

print('Starting...')

cmd = ""
cmd += "Rscript --vanilla" + " " + external_prg['script_R_before_MAPSD'] + " " + cross + "\n"

subprocess.call(cmd, shell=True)
print('Plot of the difference between somatic and germinal created...')

# Start MAPSD

cmd = ""
cmd += external_prg['script_MAPSD'] + " " + "-d" + " " + df_MAPSD + " > " + df_out_main + " " + "&" + "\n"

subprocess.call(cmd, shell=True)
print('MAPSD main completed...')

df_table = pd.read_table(df_MAPSD, header = None)

cmd = ""
for scaffold in df_table.loc[:,0].unique():
    df_subtable = df_table[df_table.loc[:,0] == scaffold]
    nameFile_df_subtable = pathways_dic['table'] + "/" + cross + "_" + str(scaffold)
    nameFile_df_subtable_BS = pathways_dic['table'] + "/" + cross + "_" + str(scaffold) + ".BS"
    
    df_subtable.to_csv(nameFile_df_subtable, header=None, index=None, sep='\t', mode='w')
    cmd += external_prg['script_MAPSD'] + " " + "-d" + " " + nameFile_df_subtable + " " + "-b 1000" + " > " + nameFile_df_subtable_BS + " " + "&" + "\n"

cmd += "wait" + "\n"

subprocess.call(cmd, shell=True)
print('MAPSD bootstrap by scaffold completed...')

cmd = ""
for scaffold in df_table.loc[:,0].unique():
	nameFile_df_subtable = pathways_dic['table'] + "/" + cross + "_" + str(scaffold)
	nameFile_df_subtable_W = pathways_dic['table'] + "/" + cross + "_" + str(scaffold) + ".W"
	cmd += external_prg['script_MAPSD'] + " " + "-d" + " " + nameFile_df_subtable + " " + "-w 100000 -e 1e-2 -t 1e-3" + " > " + nameFile_df_subtable_W + " " + "&" + "\n"

cmd += "wait" + "\n"

subprocess.call(cmd, shell=True)
print('MAPSD windows by scaffold completed...')

cmd = ""
cmd += "cat" + " " + pathways_dic['table'] + "/" + cross + "_" + "scaffold" + "*" + ".BS" + " > " + df_out_bootstrap + "\n"
cmd += "rm" + " " + pathways_dic['table'] + "/" + cross + "_" + "scaffold" + "*" + ".BS" + "\n"

cmd += "cat" + " " + pathways_dic['table'] + "/" + cross + "_" + "scaffold" + "*" + ".W" + " > " + df_out_windows + "\n"
cmd += "rm" + " " + pathways_dic['table'] + "/" + cross + "_" + "scaffold" + "*" + ".W" + "\n"

cmd += "rm" + " " + pathways_dic['table'] + "/" + cross + "_" + "scaffold" + "*" + "\n"
cmd += '\n' + 'wait' + '\n'

subprocess.call(cmd, shell=True)
print('MAPSD completed...')

# Creation of plots for MAPSD

cmd = ""
cmd += "Rscript --vanilla" + " " + external_prg['script_R_after_MAPSD'] + " " + cross + "\n"

subprocess.call(cmd, shell=True)
print('End !')
