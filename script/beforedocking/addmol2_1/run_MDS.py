#!/usr/bin/env python
import sys,os
import string
import subprocess
run_dir = "/nfs/amino-home/wenyizng/dock6/dock6/code/MC_simulation11/beforedocking/addmol2_1/"
if __name__ =="__main__":
	#fpin=open("../../BSP_lig.txt","r")
	fpin=open("../../coach_work","r")
	#fpin = open("../../DUDE","r")
	fpin_reader=fpin.read()
	fpin.close()
	for line in fpin_reader.splitlines():
		jobname = line.split('\t')[0]
		print jobname
		#command1 = "cd " + tofile + jobname + "/"
		command2 = run_dir + "MDS.py" + " " + jobname
		stdout,stderr = subprocess.Popen(';'.join([command2]),shell=True).communicate()

