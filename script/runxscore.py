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
xscorefolder =scriptdir +"/afterdocking/"
#######################run folder################################


if __name__ == "__main__":
	print jobdir
	#check the simulation conformations file
	if not os.path.exists(jobdir + "site1.mol2"):
		print "no simulation result site1.mol2 file!"
		exit()
		
	command1 = "cp " + xscorefolder + "/score.input " + jobdir
	stdout,stderr = subprocess.Popen(';'.join([command1]),shell=True).communicate()
	
	# run xscore 
	command1 = "cd " + jobdir
	command2 = xscorefolder + "/xscore" + " " + "score.input"
	stdout,stderr = subprocess.Popen(';'.join([command1,command2]),shell=True).communicate()
	print "xscore finished"
