#!/usr/bin/env python

import sys,os
import string
import subprocess
#######################job name #################################
jobname = sys.argv[1]
scriptdir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(scriptdir)

jobdir = rootdir + "/example/" + jobname + "/"
#######################job name #################################

#######################run folder################################
docking =scriptdir + "/docking/"
#######################run folder################################

if __name__ == "__main__":
	print jobdir
	grid = open(jobdir + "/grid_output","r")
	gridreader=grid.read()
	grid.close()
	sphere=open(jobdir + "/sphere_result","r")
	spherereader=sphere.read()
	sphere.close()
	dockingnum=list()
	gridfile=list()
	spherefile=list()
	for line in gridreader.splitlines():
		gridfile.append(line.strip().split(".bmp")[0])
		dockingnum.append(line.strip().split(".bmp")[0][-1])
	
	for line in spherereader.splitlines():
		spherefile.append(line.strip())
	
	outputname = "site1"
	command0 ="cd " + jobdir 	
	#edock_rigid_webserver outputname mintemperature maxtemperature energy_grid_file shperefile vdw_weight bindsite_weight
	command1 = docking+"/edock_rigid_webserver" + " " + outputname + " " + "1 60" + " " + gridfile[0] + " " + spherefile[0] + " 1 1 " + "\n"
	stdout,stderr = subprocess.Popen(';'.join([command0,command1]),shell=True).communicate()



