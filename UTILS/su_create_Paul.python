#!/usr/bin/env python
#
#  Python code to generate Seismic Unix files for Ux and Uz and plot them
#
from os import *
import os.path, string
#
def createSU():
	SEM=getcwd()
	filename=SEM+'/DATA/Par_file'
	#
	# Variables to be put in SU header
	# 
	if path.exists(filename):
		variables=['nt','dt','sismostype']
	else:
		print 'No Par_file found !'
		return 
	#     
	# Open the file and get the lines
	# 
	f = file(filename,'r')
	lignes= f.readlines()
	f.close()
	# Get the title
	for ligne in lignes:
		lsplit=string.split(ligne)
		if lsplit!= []:
			if lsplit[0]=='title':
		   		title=' '.join(lsplit[2:])
		   		break 
	print '#'*50
	print '#  SU file creation for '+title
	print  '#'*50
	#  Get the variables
	for var in variables:
		for ligne in lignes:
			lsplit=string.split(ligne)
			if lsplit!= []:
				if lsplit[0]==var:
					exec var+'='+string.replace(string.split(''.join(ligne))[2],'d','e') 
					break
	# 
	print sismostype
	chdir(SEM+'/OUTPUT_FILES')
	labels='label1="Time" label2="Receivers"'
	# Create headers and Su file
	if sismostype==4:
		ordres=['suaddhead < pressure_file_single.bin ns='+str(nt)+' | sushw key=dt a='+str(int(dt*1e6))+' > pressure_file.su']
		ordres.append('suxwigb < pressure_file.su perc=96 '+labels+' title=" Pressure : '+title+'"&')
	else:
		ordres=['suaddhead < Ux_file_single.bin ns='+str(nt)+' | sushw key=dt a='+str(int(dt*1e6))+' > Ux_file.su']
		ordres.append('suaddhead < Uz_file_single.bin ns='+str(nt)+' | sushw key=dt a='+str(int(dt*1e6))+' > Uz_file.su')
		ordres.append('suxwigb < Ux_file.su perc=96 '+labels+' title=" Ux : '+title+'"&')
		ordres.append('suxwigb < Uz_file.su xbox=600 perc=96 '+labels+' title=" Uz : '+title+'"&')
	# 
	for i in range(len(ordres)):
		 system(ordres[i]) 
if __name__=='__main__':
	createSU()
