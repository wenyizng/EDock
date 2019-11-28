#!/usr/bin/env python

import sys,os
import string
import subprocess
#######################job name #################################
jobname = sys.argv[1]
scriptdir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(scriptdir)
jobdir = "../example/" + jobname + "/"
#######################job name #################################

#######################run folder################################
#chimera folder should be changed as your chimera bin folder
beforedocking =scriptdir + "/beforedocking/"
chimera=beforedocking +"/chimera/bin/" 
addmol2_1 = beforedocking + "addmol2_1/"
bindingsitescenter_3 = beforedocking +"bindingsitescenter_3/"
negativeimage_4 = beforedocking + "negativeimage_4/"

#######################run folder################################
if __name__ == "__main__":
	print jobdir
	#check files
	checkfiles = ['receptor.pdb','lig_charge.mol2','coach_cluster.dat']
	checkflag = 0
	for checking in checkfiles:
		if not os.path.exists(jobdir + checking):
			print checking +" is missing!!"
		else:
			checkflag = checkflag +1
	if(checkflag < 3):
		exit()
	
	#step 0 : reindex user PDB file
	command1 = "cd "+jobdir
	command2 = scriptdir +"/reindex_pdb.py 1 receptor.pdb receptor.pdb"
	stdout,stderr = subprocess.Popen(';'.join([command1,command2]),shell=True).communicate()
	
	#step 1 : add the charge for receptor.pdb
	command0 = "cp " + addmol2_1 + "dockpre.py" + " " + jobdir
	stdout,stderr = subprocess.Popen(';'.join([command0]),shell=True).communicate()
	print "cp dockpre.py finished"
	
	command1 = "cd "+jobdir
	command2 = chimera+"chimera"+" --nogui pdb:receptor.pdb dockpre.py"
	command3="mv dp.mol2 receptor.mol2"
	stdout,stderr = subprocess.Popen(';'.join([command1,command2,command3]),shell=True).communicate()
	print "add charge finished"
	
	#step 2 : read the binding sites
	command1 = "cd " + jobdir 
	command2 = beforedocking+"/bindingsitescenter_3/read_binding_site"
	stdout,stderr = subprocess.Popen(';'.join([command1,command2]),shell=True).communicate()
	print "calcute the docking center finished"
	
	#step 3 : generate the negative images
	coach_result = jobdir + "coach_result"
	fp_coach=open(coach_result,"r")
	fp_coach_reader=fp_coach.read()
	fp_coach.close()
	coachnum =list()
	output = open(jobdir + "negativeimage_result","w")
	for line in fp_coach_reader.splitlines():
		if("site" in line):
			tmparray=list()
			tmparray = line.split("\t")
			coachnum.append(int(tmparray[0][-1]))
	command0= "cd " + jobdir 
	for i in range(len(coachnum)):
		command1 =  beforedocking+"/negativeimage_4/GenerateBindingSiteNegativeImage.py" + " " + "receptor" + " " + "center" + str(coachnum[i])
		stdout,stderr = subprocess.Popen(';'.join([command0,command1]),shell=True).communicate()
		checkfile = jobdir + "negativeimage_model"+str(coachnum[i])+".pdb"
		if os.path.exists(checkfile):
			output.write("negativeimage_model"+str(coachnum[i])+".pdb"+"\n")
	output.close()
	print "generate the negative images finished"
	
	#step 4 : generate the docking box
	negima_result = jobdir + "negativeimage_result"
	fp_neg=open(negima_result,"r")
	fp_neg_reader=fp_neg.read()
	fp_neg.close()
	neg =list()
	output = open(jobdir + "box_result","w")
	command0 = "cd " + jobdir
	for line in fp_neg_reader.splitlines():
		negfile = line.split(".pdb")[0]
		neg.append(int(negfile[-1]))
		print negfile
		command1 = beforedocking+"/negativeimage_4/generate_box.py" + " " + negfile + " " + "lig_charge.mol2"
		stdout,stderr = subprocess.Popen(';'.join([command0,command1]),shell=True).communicate()
		output.write("rec_box"+str(neg[len(neg)-1])+".pdb"+"\n")
	output.close()
	print "generate the docking box"

	#step 5 : write sph
	negima_result = jobdir + "negativeimage_result"
	fp_neg=open(negima_result,"r")
	fp_neg_reader=fp_neg.read()
	fp_neg.close()
	neg =list()
	output = open(jobdir + "sphere_result","w")
	command0 = "cd " + jobdir 
	for line in fp_neg_reader.splitlines():
		neg.append(int(line.split(".pdb")[0][-1]))
		#print line
		command1 = beforedocking+"/negativeimage_4/write_sph.py" + " " + line.strip()
		stdout,stderr = subprocess.Popen(';'.join([command0,command1]),shell=True).communicate()
		output.write("select_spheres"+str(neg[len(neg)-1])+".sph"+"\n")
	output.close()
	print "write sph file finished"
	
	#step 6 : write grid file
	box_result = jobdir + "box_result"
	fp_box=open(box_result,"r")
	fp_box_reader=fp_box.read()
	fp_box.close()
	box =list()
	grid_result = open(jobdir + "grid_result","w")
	for line in fp_box_reader.splitlines():
		box.append(line.strip())
	grid = open(negativeimage_4 + "grid.in","r")
	gridreader = grid.read()
	grid.close()
	for i in box:
		newgrid = "grid" + str(int(i.split(".pdb")[0][-1])) + ".in"
		output=open(jobdir + newgrid,"w")
	
		for line in gridreader.splitlines():
			if("box_file" in line):
				boxfile=i.strip()
				output.write("box_file                       "+boxfile+"\n")
				continue
			if("score_grid_prefix" in line):
				gridoutput="grid"+str(int(i.split(".pdb")[0][-1]))
				output.write("score_grid_prefix              "+gridoutput+"\n")
				continue
			else:
				output.write(line+"\n")
		output.close()
		grid_result.write("grid"+str(int(i.split(".pdb")[0][-1]))+".in"+"\n")
	grid_result.close()
	print "write grid.in file finised"
	
	#step 7 : run grid
	grid_result = jobdir + "grid_result"
	fp_grid=open(grid_result,"r")
	fp_grid_reader=fp_grid.read()
	fp_grid.close()
	grid =list()
	output = open(jobdir + "grid_output","w")
	command0 = "cd " + jobdir 
	for line in fp_grid_reader.splitlines():
		grid.append(int(line.split(".in")[0][-1]))
		print line
		command1 = beforedocking+"/negativeimage_4/grid" + " -i " + line.strip()
		stdout,stderr = subprocess.Popen(';'.join([command0,command1]),shell=True).communicate()
		checkfile = "grid"+str(grid[len(grid)-1])+".bmp"
		if os.path.exists(jobdir + checkfile):
			output.write(checkfile+"\n")
	output.close()	
	print "run grid file finished"
	#step 8 : docking by EDock
	#step 9 : runing XSCORE
	

