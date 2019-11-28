#!/usr/bin/env python

import sys,os
import string
import subprocess

chimera="/nfs/amino-home/wenyizng/chimera/bin/"
jobdir = "/nfs/amino-home/wenyizng/dock6/dock6/code/MC11_dataset/comet/COACH3/"
if __name__ =="__main__":
	command0 = "cp dockpre.py" + " " + jobdir + sys.argv[1].strip()
	
	stdout,stderr = subprocess.Popen(';'.join([command0]),shell=True).communicate()
	print "cp dockpre.py finished"
	
	command1 = "cd "+jobdir+sys.argv[1].strip()
	command2 = chimera+"chimera"+" --nogui pdb:model1.pdb dockpre.py"
	command3="mv dp.mol2 receptor.mol2"
	stdout,stderr = subprocess.Popen(';'.join([command1,command2,command3]),shell=True).communicate()
	print "add charge finished"


