#!/usr/bin/env python
# 2017-3-2
'''
function: generate the .sph file for the negative image, it is for the maching in the eDock
exmaple: ./write_sph.py 
input:negativeimage.pdb
output: selected_shperes.sph
'''
import sys,os
import string
import math


if __name__ =="__main__":
	#fp = open("BSP_benchmark_ITASSER_new.txt","r")
	#benchmark1 BSP_benchmark_ITASSER_new
	#fp_reader = fp.read()
	#fp.close()
	#for line in fp_reader.splitlines():
	#dock_dir = dock_bench + sys.argv[1] + "/" + "select_spheres.sph"
	dock_dir = "select_spheres"+sys.argv[1].split(".pdb")[0][-1]+".sph"
	#dock_dir = "select_spheres.sph"
	#site_file = dock_bench + sys.argv[1] + "/" + "negativeimage_model.pdb"
	site_file = sys.argv[1]
	fp_sph = open(dock_dir , "w")
	fp_site = open(site_file,"r")
	fp_site_reader = fp_site.read()
	fp_site.close()
	num = 0
	for siteline in fp_site_reader.splitlines():
		if ("TER" in siteline):
			break
		num += 1
			
	fp_sph.write("DOCK spheres within 10.0 ang of ligands" + "\n")
	fp_sph.write("cluster     1   number of spheres in cluster " + str(num) + "\n")
		
	fp_site = open(site_file,"r")
	fp_site_reader = fp_site.read()
	fp_site.close()
	num = 0
	for siteline in fp_site_reader.splitlines():
		if ("TER" in siteline):
			break
		num += 1
		#print siteline[30:38].strip(' ')
		#print siteline[38:46].strip(' ')
		#print siteline[46:54].strip(' ')
		
		fp_sph.write(" " * (5-len(str(num))) + str(num))
		fp_sph.write(" " * (8-len(siteline[30:38].strip(' '))) + siteline[30:38].strip(' ') + "00")
		fp_sph.write(" " * (8-len(siteline[38:46].strip(' '))) + siteline[38:46].strip(' ') + "00")			
		fp_sph.write(" " * (8-len(siteline[46:54].strip(' '))) + siteline[46:54].strip(' ') + "00")
		fp_sph.write("   1.000 1111 0  0" + "\n")

			
		
	fp_sph.close()

	
		
