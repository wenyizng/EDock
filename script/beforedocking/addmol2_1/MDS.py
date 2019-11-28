#!/usr/bin/env python

import sys,os
import string
import subprocess

chimera="/nfs/amino-home/wenyizng/chimera/bin/"
jobdir = "/nfs/amino-home/wenyizng/dock6/dock6/code/MC11_dataset/COACH2/"

if __name__ =="__main__":
        command0 = "cp writedms.py" + " " + jobdir + sys.argv[1].strip()
	
	stdout,stderr = subprocess.Popen(';'.join([command0]),shell=True).communicate()
	print "cp writedms.py finished"
	
	command1 = "cd " + jobdir + sys.argv[1]
	command2 = chimera + "chimera --nogui receptor.mol2 writedms.py"
	stdout,stderr = subprocess.Popen(';'.join([command1,command2]),shell=True).communicate()
