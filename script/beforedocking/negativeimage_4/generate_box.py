#!/usr/bin/env python
import sys,os
import string
import math
'''
example: generate.py negativeimage
'''
from GenerateBindingSiteNegativeImage_MOD import *

def getDistance2 (r1 ,r2 ):
	return ( ( float(r1[0]) - float(r2[0]))  **2  + ( float(r1[1]) - float(r2[1]))  ** 2 + ( float(r1[2]) - float(r2[2])) ** 2) 

def readpdb(filename):
	pdb = open(filename,"r")
	pdbreader = pdb.read()
	pdb.close()
	
	atoms = list()
	atom_num = 0
	for line in pdbreader.splitlines():
		if("ATOM" in line):
			atom_num += 1 
			tmpatom = [0.0,0.0,0.0]
			tmpatom[0] = float(line[30:38].strip())
			tmpatom[1] = float(line[38:46].strip())
			tmpatom[2] = float(line[46:54].strip())
			print tmpatom
			atoms.append(tmpatom)
	
	return atoms
def drowbox(atoms,box_file):
	fp_box = open(box_file,"w")
	xmin = 1000; xmax = -1000; ymin = 1000; ymax = -1000; zmin = 1000; zmax = -1000
	
	for i in range(len(atoms)):
		#print "*****",r_lig_atoms[i]
		if (atoms[i][0] < xmin):
			xmin = float(atoms[i][0])
		if (atoms[i][0] > xmax):
			xmax = float(atoms[i][0])
		if (atoms[i][1] < ymin):
			ymin = float(atoms[i][1])
		if (atoms[i][1] > ymax):
			ymax = float(atoms[i][1])
		if (atoms[i][2] < zmin):
			zmin = float(atoms[i][2])
		if (atoms[i][2] > zmax):
			zmax = float(atoms[i][2])
	#print "^^^^^^^^^^^^^^^^^^^^^^^6"
	#print xmin, xmax, ymin, ymax, zmin,zmax
	spacing = 5.0
        #spacing = 16.0;
	ligdis = getLongestDistanceFromCentroidForMol2File (sys.argv[2])
	print 'ligdis is',ligdis
	'''
	if(float(ligdis) > 8.0):
		spacing = float(ligdis)
	else:
		spacing = 8.0
	'''
	#print spacing
	minx = xmin-spacing; maxx = xmax + spacing; 
	miny = ymin-spacing; maxy = ymax +spacing;
	minz = zmin - spacing; maxz = zmax + spacing;
	centerx = (minx+ maxx)/2; centery = (miny + maxy)/2; centerz = (minz + maxz)/2
	denx = maxx-minx; deny = maxy-miny; denz = maxz-minz; 
	
	print minx,maxx,miny,maxy,minz,maxz
	### write box file
	fp_box.write("HEADER    CORNERS OF BOX" + "\n")
	fp_box.write("REMARK    CENTER (X Y Z) " + (8-len(str("%.3f" %centerx)))*" " + str("%.3f" %centerx) + (8-len(str("%.3f" %centery)))*" " + str("%.3f" %centery) + (8-len(str("%.3f" %centerz)))*" " + str("%.3f" %centerz) + "\n")
	fp_box.write("REMARK    DIMENSIONS (X Y Z) " + (8-len(str("%.3f" %denx)))*" " + str("%.3f" %denx) + (8-len(str("%.3f" %deny)))*" " + str("%.3f" %deny) + (8-len(str("%.3f" %denz)))*" " + str("%.3f" %denz) + "\n")
	'''
	#1: minx miny minz	
	fp_box.write("ATOM      1  C   BOX     1     " + (7-len(str(("%.3f" %minx))))*" " + str(("%.3f" %minx)) + " " + (7-len(str(("%.3f" %miny))))*" " + str(("%.3f" %miny)) + " " + (7-len(str(("%.3f" %minz))))*" " + str(("%.3f" %minz)) + " " + "\n" )
	#2: maxx miny minz
	fp_box.write("ATOM      2  C   BOX     1     " + (7-len(str(("%.3f" %maxx))))*" " + str(("%.3f" %maxx)) + " " + (7-len(str(("%.3f" %miny))))*" " + str(("%.3f" %miny)) + " " + (7-len(str(("%.3f" %minz))))*" " + str(("%.3f" %minz)) + " " + "\n" )
	#3: maxx miny maxz
	fp_box.write("ATOM      3  C   BOX     1     " + (7-len(str(("%.3f" %maxx))))*" " + str(("%.3f" %maxx)) + " " + (7-len(str(("%.3f" %miny))))*" " + str(("%.3f" %miny)) + " " + (7-len(str(("%.3f" %maxz))))*" " + str(("%.3f" %maxz)) + " " + "\n" )
	#4: minx miny maxz
	fp_box.write("ATOM      4  C   BOX     1     " + (7-len(str(("%.3f" %minx))))*" " + str(("%.3f" %minx)) + " " + (7-len(str(("%.3f" %miny))))*" " + str(("%.3f" %miny)) + " " + (7-len(str(("%.3f" %maxz))))*" " + str(("%.3f" %maxz)) + " " + "\n" )
	#5: minx maxy minz
	fp_box.write("ATOM      5  C   BOX     1     " + (7-len(str(("%.3f" %minx))))*" " + str(("%.3f" %minx)) + " " + (7-len(str(("%.3f" %maxy))))*" " + str(("%.3f" %maxy)) + " " + (7-len(str(("%.3f" %minz))))*" " + str(("%.3f" %minz)) + " " + "\n" )
	#6: maxx maxy minz
	fp_box.write("ATOM      6  C   BOX     1     " + (7-len(str(("%.3f" %maxx))))*" " + str(("%.3f" %maxx)) + " " + (7-len(str(("%.3f" %maxy))))*" " + str(("%.3f" %maxy)) + " " + (7-len(str(("%.3f" %minz))))*" " + str(("%.3f" %minz)) + " " + "\n" )
	#7: maxx maxy maxz
	fp_box.write("ATOM      7  C   BOX     1     " + (7-len(str(("%.3f" %maxx))))*" " + str(("%.3f" %maxx)) + " " + (7-len(str(("%.3f" %maxy))))*" " + str(("%.3f" %maxy)) + " " + (7-len(str(("%.3f" %maxz))))*" " + str(("%.3f" %maxz)) + " " + "\n" )
	#8: minx maxy maxz
	fp_box.write("ATOM      8  C   BOX     1     " + (7-len(str(("%.3f" %minx))))*" " + str(("%.3f" %minx)) + " " + (7-len(str(("%.3f" %maxy))))*" " + str(("%.3f" %maxy)) + " " + (7-len(str(("%.3f" %maxz))))*" " + str(("%.3f" %maxz)) + " " + "\n" )
	fp_box.write("CONECT    1    2    4    5" + "\n" + "CONECT    2    1    3    6" + "\n" + "CONECT    3    2    4    7" + "\n" + "CONECT    4    1    3    8"  + "\n" + "CONECT    5    1    6    8" + "\n" + "CONECT    6    2    5    7" + "\n" + "CONECT    7    3    6    8" + "\n" + "CONECT    8    4    5    7" + "\n" )	
	'''
	fp_box.write("ATOM      1  C   BOX     1    " + (8-len(str(("%.3f" %minx))))*" " + str(("%.3f" %minx)) + (8-len(str(("%.3f" %miny))))*" " + str(("%.3f" %miny)) + (8-len(str(("%.3f" %minz))))*" " + str(("%.3f" %minz)) + "\n" )
	#2: maxx miny minz
	fp_box.write("ATOM      2  C   BOX     1    " + (8-len(str(("%.3f" %maxx))))*" " + str(("%.3f" %maxx))  + (8-len(str(("%.3f" %miny))))*" " + str(("%.3f" %miny)) + (8-len(str(("%.3f" %minz))))*" " + str(("%.3f" %minz)) + "\n" )
	#3: maxx miny maxz
	fp_box.write("ATOM      3  C   BOX     1    " + (8-len(str(("%.3f" %maxx))))*" " + str(("%.3f" %maxx))  + (8-len(str(("%.3f" %miny))))*" " + str(("%.3f" %miny)) + (8-len(str(("%.3f" %maxz))))*" " + str(("%.3f" %maxz)) + "\n" )
	#4: minx miny maxz
	fp_box.write("ATOM      4  C   BOX     1    " + (8-len(str(("%.3f" %minx))))*" " + str(("%.3f" %minx)) + (8-len(str(("%.3f" %miny))))*" " + str(("%.3f" %miny)) + (8-len(str(("%.3f" %maxz))))*" " + str(("%.3f" %maxz)) + "\n" )
	#5: minx maxy minz
	fp_box.write("ATOM      5  C   BOX     1    " + (8-len(str(("%.3f" %minx))))*" " + str(("%.3f" %minx)) + (8-len(str(("%.3f" %maxy))))*" " + str(("%.3f" %maxy)) + (8-len(str(("%.3f" %minz))))*" " + str(("%.3f" %minz)) + "\n" )
	#6: maxx maxy minz
	fp_box.write("ATOM      6  C   BOX     1    " + (8-len(str(("%.3f" %maxx))))*" " + str(("%.3f" %maxx)) + (8-len(str(("%.3f" %maxy))))*" " + str(("%.3f" %maxy)) + (8-len(str(("%.3f" %minz))))*" " + str(("%.3f" %minz)) + "\n" )
	#7: maxx maxy maxz
	fp_box.write("ATOM      7  C   BOX     1    " + (8-len(str(("%.3f" %maxx))))*" " + str(("%.3f" %maxx)) + (8-len(str(("%.3f" %maxy))))*" " + str(("%.3f" %maxy)) + (8-len(str(("%.3f" %maxz))))*" " + str(("%.3f" %maxz)) + "\n" )
	#8: minx maxy maxz
	fp_box.write("ATOM      8  C   BOX     1    " + (8-len(str(("%.3f" %minx))))*" " + str(("%.3f" %minx)) + (8-len(str(("%.3f" %maxy))))*" " + str(("%.3f" %maxy)) + (8-len(str(("%.3f" %maxz))))*" " + str(("%.3f" %maxz)) + "\n" )
	fp_box.write("CONECT    1    2    4    5" + "\n" + "CONECT    2    1    3    6" + "\n" + "CONECT    3    2    4    7" + "\n" + "CONECT    4    1    3    8"  + "\n" + "CONECT    5    1    6    8" + "\n" + "CONECT    6    2    5    7" + "\n" + "CONECT    7    3    6    8" + "\n" + "CONECT    8    4    5    7" + "\n" )	
	fp_box.close()


if __name__ =="__main__":
	
	negfile  = sys.argv[1] + ".pdb"
	box_file = "rec_box"+sys.argv[1][-1]+".pdb"
	#box_file = "rec_box" + ".pdb"
	atoms=list
	atoms=readpdb(negfile)
	drowbox(atoms,box_file)


