#!/usr/bin/env python
import sys,os
import string
import math


fp=open("lig_charge.mol2","r")
fpreader=fp.read()
fp.close()
fout=open("lig_center.pdb","w")
nativelig = list()
flag = False
nativeatoms = list()
center_native_temp = [0.0,0.0,0.0]
center_native = [0.0,0.0,0.0]
nativeatomsnum = 0
for nativeline in fpreader.splitlines():
	if("@<TRIPOS>ATOM" in nativeline):
		flag = True
		continue
	if("@<TRIPOS>BOND" in nativeline):	
		flag = False	
	if(flag == True and nativeline[47:48] != 'H'):
		r_native =  [0.0,0.0,0.0]
		array = nativeline.split()
		r_native[0] = float(array[2])
		r_native[1] = float(array[3])
		r_native[2] = float(array[4])
		center_native_temp[0] += r_native[0]
		center_native_temp[1] += r_native[1]
		center_native_temp[2] += r_native[2]
		nativeatomsnum += 1
		nativeatoms.append(r_native)
center_native[0] = center_native_temp[0]/nativeatomsnum
center_native[1] = center_native_temp[1]/nativeatomsnum
center_native[2] = center_native_temp[2]/nativeatomsnum
print center_native

fout.write("ATOM      1  DUM PBS     1    ")
X = '%.3f' % float(center_native[0])
Y = '%.3f' % float(center_native[1])
Z = '%.3f' % float(center_native[2])
fout.write(" " * (8-len(X)) + str(X) + " " * (8-len(Y)) + str(Y) + " " * (8-len(Z)) + str(Z) + "\n")	

