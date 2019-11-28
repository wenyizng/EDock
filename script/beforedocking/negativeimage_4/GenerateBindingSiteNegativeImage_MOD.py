#!/usr/bin/env python
# 2016-2-3 wenyi Zhang

docstring='''
shuoming
shuoming
'''
import sys,os, math, string
class Atom():
	atom_name = ""        #three characters
	atom_type = ""
	residue_name = ""     #three characters
	sequence_number = ""  #three characters
	charge = 0.0
	Emin = 0.0
	Rmin = 0.0 

        #EEF1 parameters
	volume = 0.0 
	dG_ref = 0.0 
	dG_free = 0.0
	lambda_value = 0.0
	R = [0.0,0.0,0.0]
	energy = 0.0
	def __init__(self,atom = "", type_name = "", res="", seq = "", r = ['0.0','0.0','0.0']):
		self.atom_name = atom
		self.atom_type = type_name
		self.residue_name = res		
		self.sequence_number = seq
		self.R = ['0.0','0.0','0.0']
		self.R[0] = r[0]    
		self.R[1] = r[1]
		self.R[2] = r[2]

	def Atom(self):

		R[0] = 0.0
		R[1] = 0.0
		R[2] = 0.0
	def Atom(self,r):
		R[0] = r[0]
		R[1] = r[1]
		R[2] = r[2]
	'''
	def make_self(self):
		return Atom(
			atom=str(self.atom_name),
			type_name = str(self.atom_type),
			res=str(self.residue_name),
			seq =str(self.sequence_numer),
			r = list(self.R))
	'''
	
	

def getDistance (r1 ,r2 ):
# Calculate a distance between two atoms
	return math.sqrt( ( float(r1[0]) - float(r2[0]) ) * ( float(r1[0]) - float(r2[0]) ) + ( float(r1[1]) - float(r2[1]) ) * ( float(r1[1]) - float(r2[1]) ) + ( float(r1[2]) - float(r2[2]) ) * ( float(r1[2]) - float(r2[2]) ) )


def getAtomCoordinatesFromPDBFile (file_name):
# Get the coordinates of receptor from a PDB file
	vAtoms = list() #Atom()
	line = ""
	r = [0.0, 0.0, 0.0]
	fileinputstream = open(file_name, "r")
	lineall = fileinputstream.read()
	fileinputstream.close()
	
	for line in lineall.splitlines():
		if line.startswith("ATOM  "):
			r[0] = float(line[30:38].strip(' '))
			r[1] = float(line[38:46].strip(' '))
			r[2] = float(line[46:54].strip(' '))
			vAtoms.append( Atom( line[13:16], "" ,line[17:20], line[23:26], r ))
	
		if "TER" in line or "END" in line:
			break
	
	return vAtoms


def getCoordinatesFromPDBFile (file_name):
# get the coordinates from PDB file 
	vAtoms = list() #Atom()
	line = ""
	r = [0.0, 0.0, 0.0]
	fileinputstream = open(file_name, "r")
	lineall = fileinputstream.read()
	fileinputstream.close()
	
	for line in lineall.splitlines():
		if line.startswith("ATOM  ") or line.startswith("HETATM"):
			r[0] = float(line[30:38].strip(' '))
			r[1] = float(line[38:46].strip(' '))
			r[2] = float(line[46:54].strip(' '))
			vAtoms.append( Atom( line[13:16], "" ,line[17:20], line[23:26], r ))
		
		if "TER" in line or "END" in line:
			break
	return vAtoms 


def getHeavyAtomsFromFile (file_name):
#Read heavy atoms from PDB file
	vAtoms = list() #Atom()
	line = ""
	r = [0.0, 0.0, 0.0]
	fileinputstream = open(file_name, "r")
	lineall = fileinputstream.read()
	fileinputstream.close()
	num =0
	for line in lineall.splitlines():
		if line.startswith("ATOM  ") and not line[12:16].strip(" ").startswith('H') :
			
			r[0] = float(line[30:38].strip(' '))
			r[1] = float(line[38:46].strip(' '))
			r[2] = float(line[46:54].strip(' '))
			vAtoms.append( Atom( line[13:16], "" ,line[17:20], line[23:26], r ))
		#if "TER" in line:
		#	break
	return vAtoms


def setSelectedAtomsFromFile (atom_type, reader, vAtoms):
#Read heavy atoms fro PDB file
	line = ""
	r = [0.0,0.0,0.0]
	
	for line in reader.splitlines():
		if line.startswith("ATOM  ") and line[13:16].strip(' ').upper() == atom_type.upper(): #contain space 
			r[0] = float(line[30:38].strip(' '))
			r[1] = float(line[38:46].strip(' '))
			r[2] = float(line[46:54].strip(' '))
			vAtoms.append( Atom( line[13:16], "" ,line[17:20], line[23:26], r ))
		if "TER" in line:
			break

	

def isPassed(atom_type, res_type):
	if (atom_type.strip(' ').upper()== "CA" or 
	(res_type.strip(' ').upper() == "ALA" and atom_type.strip(' ').upper == "CB") or 
	(res_type.strip(' ').upper() == "ARG" and atom_type.strip(' ').upper == "NH1") or 
	(res_type.strip(' ').upper() == "ASN" and atom_type.strip(' ').upper == "ND2") or 
	(res_type.strip(' ').upper() == "ASP" and atom_type.strip(' ').upper == "OD1") or 
	(res_type.strip(' ').upper() == "CYS" and atom_type.strip(' ').upper == "SG") or 
	(res_type.strip(' ').upper() == "GLN" and atom_type.strip(' ').upper == "NE2") or 
	(res_type.strip(' ').upper() == "GLU" and atom_type.strip(' ').upper == "OE1") or 
	(res_type.strip(' ').upper() == "HIS" and atom_type.strip(' ').upper == "NE2") or 
	(res_type.strip(' ').upper() == "ILE" and atom_type.strip(' ').upper == "CD1") or 
	(res_type.strip(' ').upper() == "LEU" and atom_type.strip(' ').upper == "CD1") or 
	(res_type.strip(' ').upper() == "LYS" and atom_type.strip(' ').upper == "NZ") or 
	(res_type.strip(' ').upper() == "MET" and atom_type.strip(' ').upper == "CE") or 
	(res_type.strip(' ').upper() == "PHE" and atom_type.strip(' ').upper == "CZ") or 
	(res_type.strip(' ').upper() == "PRO" and atom_type.strip(' ').upper == "CG") or 
	(res_type.strip(' ').upper() == "SER" and atom_type.strip(' ').upper == "OG") or 
	(res_type.strip(' ').upper() == "THR" and atom_type.strip(' ').upper == "OG1") or 
	(res_type.strip(' ').upper() == "TRP" and atom_type.strip(' ').upper == "NE1") or 
	(res_type.strip(' ').upper() == "TYR" and atom_type.strip(' ').upper == "OH") or 
	(res_type.strip(' ').upper() == "VAL" and atom_type.strip(' ').upper == "CG1")):
		return True
	else:
		return False

def setFilteredAtomsFromFile (reader, vAtoms):
	line = ""
	atom_type = ""
	res_type = ""
	r = [0.0,0.0,0.0]
	for line in reader.splitlines():
		if "TER" in line:
			break
		
		atom_type = line[13:16].strip(' ')
		res_type = line[17:20]
		if (isPassed (atom_type, res_type)):
			r[0] = float(line[30:38].strip(' '))
			r[1] = float(line[38:46].strip(' '))
			r[2] = float(line[46:54].strip(' '))
			vAtoms.append( Atom( line[13:16], "" ,line[17:20], line[23:26], r ))



def getCACenterofMass (vAtoms):
#Calculate the Center of Geometry
	R = [0.0,0.0,0.0]
	atom_mass = 12.0
	mass =0
	atom_name = ""
	for i in range(len(vAtoms)):
		atom_name = vAtoms[i].atom_name
		if ("CA" in atom_name):
			R[0] += atom_mass * vAtom[i].R[0]
			R[1] += atom_mass * vAtom[i].R[1]
			R[2] += atom_mass * vAtom[i].R[2]
			mass += atom_mass


	R[0] = R[0] / mass
	R[1] = R[1] / mass
	R[2] = R[2] / mass
	
	return R
			

def horizotalTranslation (VAtoms, move):
# Traslate the trajectory in horizontal direction	
	for i in range(len(VAtoms)):
		vAtoms[i].R[0] = vAtoms[i].R[0] - move[0]
		vAtoms[i].R[1] = vAtoms[i].R[1] - move[1]
		vAtoms[i].R[2] = vAtoms[i].R[2] - move[2]
		
def horizotalTranslation (point, move):
	point[0] = point[0] - move[0]
	point[1] = point[1] - move[1]
	point[2] = point[2] - move[2]


def setRotationMatrix (Rx, Ry, Rz, x_angle, y_angle, z_angle):
#set rotation Matrix at given x, y and z angle
	alpha = float(x_angle * 3.141592/180.0)
	beta = float(y_angle * 3.141592/180.0)
	gamma = float(z_angle * 3.141592/180.0)
	
	#Rotation Matrix at the axis of X
	Rx[0][0] = 1; Rx[0][1] = 0; Rx[0][2] = 0
	Rx[1][0] = 0; Rx[1][1] = math.cos(alpha); Rx[1][2] = -math.sin(alpha)
	Rx[2][0] = 0; Rx[2][1] = math.sin(alpha); Rx[2][2] = math.cos(alpha)
	
	#Rotation Matrix at the axis of Y
	Ry[0][0] = math.cos(beta); Ry[0][1] = 0; Ry[0][2] = math.sin(beta)
	Ry[1][0] = 0; Ry[1][1] = 1; Ry[1][2] = 0
	Ry[2][0] = -math.sin(beta); Ry[2][1] = 0; Ry[2][2] = math.cos(beta)
	
	#Rotation Matrix at the axis of Z
	Rz[0][0] = math.cos(gamma); Rz[0][1] = -math.sin(gamma); Rz[0][2] = 0
	Rz[1][0] = math.sin(gamma); Rz[1][1] = math.cos(gamma); Rz[1][2] = 0
	Rz[2][0] = 0; Rz[2][1] = 0; Rz[2][2] = 1


def setRotatedCoordinates(Rx, Ry, Rz, vAtoms ):
#Set coordinates roatated by a given rotational matrix
	temp = [0.0,0.0,0.0]
	for i in range(len(vAtoms)):
		for j in range(3):
			temp[j] = Rx[j][0] * vAtoms[i].R[0] + Rx[j][1] * vAtoms[i].R[1] + Rx[j][2] * vAtoms[i].R[2] 
		vAtoms[i].R[0] = temp[0]; vAtoms[i].R[1] = temp[1]; vAtoms[i].R[2] = temp[2]

		for j in range(3):
			temp[j] = Ry[j][0] * vAtoms[i].R[0] + Ry[j][1] * vAtoms[i].R[1] + Ry[j][2] * vAtoms[i].R[2] 
		vAtoms[i].R[0] = temp[0]; vAtoms[i].R[1] = temp[1]; vAtoms[i].R[2] = temp[2]

		for j in range(3):
			temp[j] = Rz[j][0] * vAtoms[i].R[0] + Rz[j][1] * vAtoms[i].R[1] + Rz[j][2] * vAtoms[i].R[2] 
		vAtoms[i].R[0] = temp[0]; vAtoms[i].R[1] = temp[1]; vAtoms[i].R[2] = temp[2]


def setRotatedCoordinates(Rx, Ry, Rz, vAtoms, vAtoms_rotated):
	temp = [0.0,0.0,0.0]
	for i in range(len(vAtoms)):
		for j in range(3):
			temp[j] = Rx[j][0] * vAtoms[i].R[0] + Rx[j][1] * vAtoms[i].R[1] + Rx[j][2] * vAtoms[i].R[2] 
		vAtoms_rotated[i].R[0] = temp[0]; vAtoms_rotated[i].R[1] = temp[1]; vAtoms_rotated[i].R[2] = temp[2]		
		for j in range(3):
			temp[j] = Ry[j][0] * vAtoms_rotated[i].R[0] + Ry[j][1] * vAtoms_rotated[i].R[1] + Ry[j][2] * vAtoms_rotated[i].R[2] 
		vAtoms_rotated[i].R[0] = temp[0]; vAtoms_rotated[i].R[1] = temp[1]; vAtoms_rotated[i].R[2] = temp[2]

		for j in range(3):
			temp[j] = Rz[j][0] * vAtoms_rotated[i].R[0] + Rz[j][1] * vAtoms_rotated[i].R[1] + Rz[j][2] * vAtoms_rotated[i].R[2] 
		vAtoms_rotated[i].R[0] = temp[0]; vAtoms_rotated[i].R[1] = temp[1]; vAtoms_rotated[i].R[2] = temp[2]	
	
	
def setRotatedCoordinates(Rx, Ry, Rz, point):	
	temp = [0.0,0.0,0.0]
	# Rotation at the axis of X
	for j in range(3):
		temp[j] = Rx[j][0] * point[0] + Rx[j][1] * point[1] + Rx[j][2] * point[2] 
	point[0] = temp[0]; point[1] = temp[1]; point[2] = temp[2]
	#Rotation at the axis of Y
	for j in range(3):
		temp[j] = Ry[j][0] * point[0] + Ry[j][1] * point[1] + Ry[j][2] * point[2] 
	point[0] = temp[0]; point[1] = temp[1]; point[2] = temp[2]
	#Rotation at the axis of YRotation at the axis of Y
	for j in range(3):
		temp[j] = Rz[j][0] * point[0] + Rz[j][1] * point[1] + Rz[j][2] * point[2] 
	point[0] = temp[0]; point[1] = temp[1]; point[2] = temp[2]
		
	
def getRMSD (vAtoms1, vAtoms2):
# Calculating RMSD
	sum_value = 0.0
	for i in range(len(vAtoms1)):
		sum_value += (vAtoms1[i].R[0]-vAtoms2[i].R[0]) ** 2 + (vAtoms1[i].R[1]-vAtoms2[i].R[1]) ** 2 + (vAtoms1[i].R[2]-vAtoms2[i].R[2]) ** 2
	
	return math.sqrt(sum_value/float(len(vAtoms1)))
	
def writeAtomsInPDBFormat (vAtoms, chain_id, output):
# Write the PDB file format // format is wrong
	for i in range(len(vAtoms)):
		serial_num = str(i+1)
		output.write( "ATOM  " )
		output.write(' '*(5-len(serial_num))+serial_num)
		output.write( vAtoms[i].atom_name )
		output.write( " " )
		output.write( vAtoms[i].residue_name )
		output.write( " " )
		output.write( chain_id )
		output.write( "    " )
		output.write( vAtoms[i].sequence_number )
		output.write( "    " )
		X = '%.3f' % float(vAtoms[i].R[0])
		Y = '%.3f' % float(vAtoms[i].R[1])
		Z = '%.3f' % float(vAtoms[i].R[2])
		output.write(" " * (8-len(X)) + X + " " * (8-len(Y))+ Y + " " * (8-len(Z)) + Z + '\n')
		output.write("TER" + "\n")
def writeNegativeimageInPDBFormat (vGrids, output):
	for i in range(len(vGrids)):
		serial_num=str(i+1)

		output.write( "ATOM  "+' '*(5-len(serial_num))+serial_num+" DUM PBS  "+ \
			' '*(5-len(serial_num))+serial_num+'    ')
		X = '%.3f' % float(vGrids[i].R[0])
		Y = '%.3f' % float(vGrids[i].R[1])
		Z = '%.3f' % float(vGrids[i].R[2])
		output.write(" " * (8-len(X)) + X + " " * (8-len(Y))+ Y + " " * (8-len(Z)) + Z + '\n')
	output.write("TER" + "\n")

def writeXYZinPDBFormat (vAtoms, output):	
# Write the PDB file format
	for i in range(len(mRays)):
		serial_num=str(i+1)

		output.write( "ATOM  "+' '*(5-len(serial_num))+serial_num+" CA  GLY  "+ \
			' '*(5-len(serial_num))+serial_num+'    ')
		X = '%.3f' % float(mRays[i][0])
		Y = '%.3f' % float(mRays[i][1])
		Z = '%.3f' % float(mRays[i][2])
		output.write(" " * (8-len(X)) + X + " " * (8-len(Y))+ Y + " " * (8-len(Z)) + Z + '\n')
	output.write("TER" + "\n")
		
def writeHeteroAtomInPDBFormat (atom, output):
	output.write( "HETATM" )
	output.write( "       " )
	output.write( atom.atom_name )
	output.write( " " )
	output.write( atom.residue_name )
	output.write( " " )
	output.write( "A" )
	output.write( " " )
	output.write( atom.sequence_number )
	output.write( "    " )
	X = '%.3f' % float(atom.R[0])
	Y = '%.3f' % float(atom.R[1])
	Z = '%.3f' % float(atom.R[2])
	output.write(" " * (8-len(X)) + X + " " * (8-len(Y))+ Y + " " * (8-len(Z)) + Z)
	output.write( "  1.00  0.00           " )
	output.write( "\n" )
	
def writeAtomsInSDFormat(title, vAtoms, output):
	output.write( title + "\n" )
	output.write( "                    3D" + "\n" )
	output.write( "Structure written by HSLee" + "\n" )
	output.write(" " * (10-len(vAtoms)) + str(len(vAtoms))+ " 0  0     0  0  0  0  0  0999 V2000" + "\n")
	for i in range(len(vAtoms)):
		#print vAtoms[i].R
		#print i
		X = '%.4f' % float(vAtoms[i].R[0])
		Y = '%.4f' % float(vAtoms[i].R[1])
		Z = '%.4f' % float(vAtoms[i].R[2])
		output.write(" " * (10-len(X)) + X + " " * (10-len(Y))+ Y + " " * (10-len(Z)) + Z)
		output.write( " " + vAtoms[i].atom_name + " 0  0  0  0  0  0\n")
		#import pdb; pdb.set_trace()
	output.write("M  END" + "\n")
	output.write("\n$$$$\n")

def getCentroidFromPDBFile (coord_file):
#Caculate the Centroid from a PDB file
	line = ""
	R = [0.0,0.0,0.0]
	total_atom_num = 0 
	fileinput = open(coord_file, "r")
	inputall = fileinput.read()
	fileinput.close()
	for line in inputall.splitlines():
		if (line.startswith("ATOM  ") or line.startswith("HETATM")) and (not line[12:16].strip(' ').startswith("H")):
			R[0] += float(line[30:38].strip(' '))
			R[1] += float(line[38:46].strip(' '))
			R[2] += float(line[46:54].strip(' '))
			total_atom_num += 1
		if ("TER" in line or "END" in line):
			break
	
	R[0] = R[0]/total_atom_num
	R[1] = R[1]/total_atom_num
	R[2] = R[2]/total_atom_num
	return R

def getCentroidFromPDBFile (coord_file, residues):
#Calculate the Centroid from a PDB file using selected residues
	line = ""
	R = [0.0,0.0,0.0]
	res_num = 0
	isSelectResidue = False
	array_length = len(residues)
	total_atom_num = 0
	fileinput = open(coord_file, "r")
	inputall = fileinput.read()
	fileinput.close()
	for line in inputall.splitlines():
		if ("TER" in line or "END" in line):
			break
		res_num = int(line[22:26].strip(' '))
		for i in range(array_length):
			if(residues[i]==res_num):
				isSelectResidue = True
		if(line.startswith("ATOM  ") and isSelectResidue):
			R[0] += float(line[30:38].strip(' '))
			R[1] += float(line[38:46].strip(' '))
			R[2] += float(line[46:54].strip(' '))
			total_atom_num += 1
	
	
	R[0] = R[0]/total_atom_num
	R[1] = R[1]/total_atom_num
	R[2] = R[2]/total_atom_num
	return R

def getCentroidFromVector(vAtoms):
	R = [0.0,0.0,0.0]
	num_atoms = len(vAtoms)
	for i in range(num_atoms):
		R[0] += vAtoms[i].R[0]
		R[1] += vAtoms[i].R[1]
		R[2] += vAtoms[i].R[2]
		
	R[0] = R[0] / num_atoms
	R[1] = R[1] / num_atoms
	R[2] = R[2] / num_atoms
	
	return R

def getLongestXYZCoordinatesFromVector (vAtoms,cenroid):
	R = [[0.0,0.0],[0.0,0.0],[0.0,0.0]]
	dist_0 = 0.0 
	dist_1 = 0.0
	num_atoms = len(vAtoms)
	for i in range(num_atoms):
		for j in range(3):
			dist_0  = vAtoms[i].R[j] - centroid[j]
			if (dist_0 > R[j][0]):
				R[j][0] = dist_0
			
			dist_1  =  centroid[j] - vAtoms[i].R[j]
			if (dist_1 > R[j][1]):
				R[j][1] = dist_1

	return R

def getLongestDistanceFromCentroidForPDBFile(coord_file):
	distance = 0.0 
	longest_distance = 0.0
	line = ""
	centroid = [0.0,0.0,0.0]
	num_atoms = 0
	R = [0.0,0.0,0.0]
	
	fileinput = open(coord_file, "r")
	inputall = fileinput.read()
	fileinput.close()
	for line in inputall.splitlines():
		if line.startswith("ATOM  ") or line.startswith("HETATM") :
			num_atoms += 1
			centroid[0] += float(line[30:38].strip(' '))
			centroid[1] += float(line[38:46].strip(' '))
			centroid[2] += float(line[46:54].strip(' '))
		if "TER" in line or "END" in line:
			break
	
	centroid[0] = centroid[0]/num_atoms
	centroid[1] = centroid[1]/num_atoms
	centroid[2] = centroid[2]/num_atoms
	
	fileinput = open(coord_file, "r")
	inputall = fileinput.read()
	fileinput.close()
	for line in inputall.splitlines():
		if line.startswith("ATOM  ") or line.startswith("HETATM") :
			R[0] = float(line[30:38].strip(' '))
			R[1] = float(line[38:46].strip(' '))
			R[2] = float(line[46:54].strip(' '))
			distance = math.sqrt((R[0]-centroid[0])**2 + (R[1]-centroid[1])**2 + (R[2]-centroid[2])**2)
			if(distance > longest_distance):
				longest_distance = distance
			
		if "TER" in line or "END" in line:
			break

	return longest_distance

def getLongestDistanceFromCentroidForSDFFile (coord_file):
	distance = 0.0
	longest_distance = 0.0
	line = ""
	centroid = [0.0,0.0,0.0]
	num_atoms = 0
	R = [0.0,0.0,0.0]
	total_num_atoms = 0
	num_heavy_atoms = 0
	fileinput = open(coord_file,"r")
	inputline = fileinput.readline()
	inputline = fileinput.readline()
	inputline = fileinput.readline()
	inputline = fileinput.readline()
	array = inputline.split(' ')
	#print array
	total_num_atoms = int(array[1])
	inputline = fileinput.readline()
	while inputline:
		num_atoms += 1
		if (num_atoms > total_num_atoms):
			break
		if (len(inputline)<31):
			break
		if (inputline[31] != "H"):
			num_heavy_atoms += 1
			centroid[0] += float(inputline[0:10].strip(' '))
			centroid[1] += float(inputline[10:20].strip(' '))
			centroid[2] += float(inputline[20:30].strip(' '))
		inputline = fileinput.readline()

	fileinput.close()
	centroid[0] = centroid[0]/num_heavy_atoms
	centroid[1] = centroid[1]/num_heavy_atoms
	centroid[2] = centroid[2]/num_heavy_atoms
	num_atoms = 0;
	fileinput = open(coord_file,"r")
	inputline = fileinput.readline()
	inputline = fileinput.readline()
	inputline = fileinput.readline()
	inputline = fileinput.readline()
	inputline = fileinput.readline()
	while inputline:
		num_atoms += 1
		if (num_atoms > total_num_atoms):
			break
		if (len(inputline)<31):
			break
		if (inputline[31] != "H"):
			R[0] = float(inputline[0:10].strip(' '))
			R[1] = float(inputline[10:20].strip(' '))
			R[2] = float(inputline[20:30].strip(' '))
			distance = math.sqrt((R[0]-centroid[0])**2 + (R[1]-centroid[1])**2 + (R[2]-centroid[2])**2)		
			if(distance > longest_distance):
				longest_distance = float(distance)
		inputline = fileinput.readline()

	fileinput.close()
	return longest_distance
def getLongestDistanceFromCentroidForMol2File (coord_file):
	distance = 0.0
	longest_distance = 0.0
	line = ""
	centroid = [0.0,0.0,0.0]
	num_atoms = 0
	R = [0.0,0.0,0.0]
	total_num_atoms = 0
	num_heavy_atoms = 0
	fileinput = open(coord_file,"r")
	inputall = fileinput.read()
	fileinput.close()
	flag = False
	for line in inputall.splitlines():
		if("@<TRIPOS>ATOM" in line):
			flag = True
			continue
		if(flag == True):
			if("@<TRIPOS>BOND" in line):
				flag = False
				break
			else:
				if (line[47] != "H"):
					#print line[47]
					num_heavy_atoms += 1
					#centroid[0] += float(line[18:27].strip(' '))
					#centroid[1] += float(line[28:37].strip(' '))
					#centroid[2] += float(line[38:47].strip(' '))	
					centroid[0] += float(line.split()[2])
					centroid[1] += float(line.split()[3])
					centroid[2] += float(line.split()[4])
				
		else:
			continue
	
	centroid[0] = centroid[0]/num_heavy_atoms
	centroid[1] = centroid[1]/num_heavy_atoms
	centroid[2] = centroid[2]/num_heavy_atoms
	for line in inputall.splitlines():
		if("@<TRIPOS>ATOM" in line):
			flag = True
			continue
		if(flag == True):
			if("@<TRIPOS>BOND" in line):
				flag = False
				break
			else:
				#print line[47]
				if (line[47] != "H"):
					#R[0] = float(line[18:27].strip(' '))
					#R[1] = float(line[28:37].strip(' '))
					#R[2] = float(line[38:47].strip(' '))

					R[0] = float(line.split()[2])
					R[1] = float(line.split()[3])
					R[2] = float(line.split()[4])

					distance = math.sqrt((R[0]-centroid[0])**2 + (R[1]-centroid[1])**2 + (R[2]-centroid[2])**2)		
					if(distance > longest_distance):
						longest_distance = float(distance)
				
				
		else:
			continue
	#print longest_distance
	#print num_heavy_atoms
	return longest_distance
		

def getLongestDistanceFromCentroid (vAtoms, centroid):
	distance = 0.0
	longest_distance = -10
	R = [0.0,0.0,0.0]
	total_atom_num = 0
	for i in range(len(vAtoms)):
		R[0] = vAtoms[i].R[0]
		R[1] = vAtoms[i].R[1]
		R[1] = vAtoms[i].R[1]
		distance = math.sqrt((R[0]-centroid[0])**2 + (R[1]-centroid[1])**2 + (R[2]-centroid[2])**2)
		if(distance > longest_distance):
				longest_distance = distance

	return longest_distance 

def getCoordinateofHydrogenAtom (r_d, r_a):
	d = [0.0,0.0,0.0]
	d[0] = r_a[0] - r_d[0]	
	d[1] = r_a[1] - r_d[1]	
	d[2] = r_a[2] - r_d[2]
	distance = math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2) 
	hydrogen_atom_coord = [0.0,0.0,0.0]
	hydrogen_atom_coord[0] = d[0] / distance
	hydrogen_atom_coord[1] = d[1] / distance
	hydrogen_atom_coord[2] = d[2] / distance	
	
	hydrogen_atom_coord[0] = hydrogen_atom_coord[0] + r_d[0]
	hydrogen_atom_coord[1] = hydrogen_atom_coord[1] + r_d[1]
	hydrogen_atom_coord[2] = hydrogen_atom_coord[2] + r_d[2]
	
	return hydrogen_atom_coord

#########################################
#check physical characters of atom groups
#########################################

def isDonorGroup (atom):
	if (atom.atom_type.strip(' ').upper() == "NH1" or atom.atom_type.strip(' ').upper() == "NH2"):
		return True
	else:
		return False

def isAcceptorGroup (atom):
	if (atom.atom_type.strip(' ').upper() == "NR" or atom.atom_type.strip(' ').upper() == "O"):
		return True
	else:
		return False

def isCationGroup (atom):
	if (atom.atom_type.strip(' ').upper() == "NH3" or atom.atom_type.strip(' ').upper() == "NC2"):
		return True
	else:
		return False

def isAnionGroup (atom):
	if (atom.atom_type.strip(' ').upper() == "OC"):
		return True
	else:
		return False

def isRingGroup (atom):
	if (atom.atom_type.strip(' ').upper() == "CR1E" or atom.atom_type.strip(' ').upper() == "CR"):
		return True
	else:
		return False

def isHydrophobicGroup (atom):
	if (atom.atom_type.strip(' ').upper() == "CH2E" or atom.atom_type.strip(' ').upper() == "CH3E"):
		return True
	else:
		return False

def isHydroxylGroup (atom):
	if (atom.atom_type.strip(' ').upper() == "OH1"):
		return True
	else:
		return False


def setCFFToGrids (grid, receptor):
#Set probes to the grids according to the simple selection rules
	vCFFAtoms = list()
	distance_charge_from = 2.5               
	distance_charge_to = 4.5
	distance_hydrophobic_contact_from = 2.5
	distance_hydrophobic_contact_to = 4.5        
	distance = 0.0
	distance_nearest = 100
	#grid_atom, receptor_atom
	for i in range(len(grid)):
		#grid_atom = list(grid[i])
		#print "i is : ",i
		for j in range(len(receptor)):
			#receptor_atom = list(receptor[j])
			distance = getDistance(grid[i].R, receptor[j].R)
			#print distance
			if (distance < distance_nearest):
				distance_nearest = float(distance) 
				#print "distance_nearest is :", distance_nearest
				nearest_atom_index  = int(j)
				#print "nearest_atom_index is: ", nearest_atom_index
		#print receptor[nearest_atom_index].atom_type
		isDonor = isDonorGroup( receptor[nearest_atom_index] )
		isAcceptor = isAcceptorGroup( receptor[nearest_atom_index]  )
		isCation = isCationGroup( receptor[nearest_atom_index]  )
		isAnion = isAnionGroup( receptor[nearest_atom_index]  )
		isRing = isRingGroup( receptor[nearest_atom_index]  )
		isHydrophobe = isHydrophobicGroup( receptor[nearest_atom_index]  )
		isHydroxyl = isHydroxylGroup( receptor[nearest_atom_index]  )
		#print isDonor, isAcceptor, isCation,isAnion,isRing,isHydrophobe,isHydroxyl
		#print isDonor,isAcceptor,isCation,isAnion,isRing,isHydrophobe,isHydroxyl
		# *** IMPLICITMILLSDEAN COLOR TYPES ***
		# Type index: 1 donor
		# Type index: 2 acceptor
		# Type index: 3 cation
		# Type index: 4 anion
		# Type index: 5 rings
		# Type index: 6 hydrophobe
		# Type index: 7 hydroxyl
		if( ( isDonor or isAcceptor or isCation or isAnion or isHydroxyl ) and 	( ( distance_nearest >= distance_charge_from ) and ( distance_nearest < distance_charge_to ) ) ):
			if( isDonor ):
				vCFFAtoms.append( Atom( "2  ", "C", "GRI", "  1", grid[i].R ) )
				#print "vCFFAtoms[%d].atom_name is : %s"  % (len(vCFFAtoms),vCFFAtoms[len(vCFFAtoms)-1].atom_name)
				#print vCFFAtoms[len(vCFFAtoms)-1].R
				
			elif( isAcceptor ): 
				vCFFAtoms.append( Atom( "1  ", "C", "GRI", "  1", grid[i].R ) )
				#print "vCFFAtoms[%d].atom_name is : %s"  % (len(vCFFAtoms),vCFFAtoms[len(vCFFAtoms)-1].atom_name)
				#print vCFFAtoms[len(vCFFAtoms)-1].R
			elif( isCation ):					
				vCFFAtoms.append( Atom( "2  ", "C", "GRI", "  1", grid[i].R ) )
				#print "vCFFAtoms[%d].atom_name is : %s"  % (len(vCFFAtoms),vCFFAtoms[len(vCFFAtoms)-1].atom_name)
				#print vCFFAtoms[len(vCFFAtoms)-1].R
				vCFFAtoms.append( Atom( "4  ", "C", "GRI", "  1", grid[i].R ) )
				#print "vCFFAtoms[%d].atom_name is : %s"  % (len(vCFFAtoms),vCFFAtoms[len(vCFFAtoms)-1].atom_name)
				#print vCFFAtoms[len(vCFFAtoms)-1].R
			elif( isAnion ):
				vCFFAtoms.append( Atom( "1  ", "C", "GRI", "  1", grid[i].R ) )
				#print "vCFFAtoms[%d].atom_name is : %s"  % (len(vCFFAtoms),vCFFAtoms[len(vCFFAtoms)-1].atom_name)
				#print vCFFAtoms[len(vCFFAtoms)-1].R
				vCFFAtoms.append( Atom( "3  ", "C", "GRI", "  1", grid[i].R ) )
				#print "vCFFAtoms[%d].atom_name is : %s"  % (len(vCFFAtoms),vCFFAtoms[len(vCFFAtoms)-1].atom_name)
				#print vCFFAtoms[len(vCFFAtoms)-1].R
			elif( isHydroxyl ):
				vCFFAtoms.append( Atom( "7  ", "C", "GRI", "  1", grid[i].R ) )
				#print grid[i].atom_name
				#print grid[i].R
				#print "vCFFAtoms[%d].atom_name is : %s"  % (len(vCFFAtoms),vCFFAtoms[len(vCFFAtoms)-1].atom_name)
				#print vCFFAtoms[len(vCFFAtoms)-1].R
		elif( isRing and ( ( distance_nearest >= distance_hydrophobic_contact_from ) and ( distance_nearest < distance_hydrophobic_contact_to ) ) ):
			vCFFAtoms.append(Atom( "5  ", "C", "GRI", "  1", grid[i].R ) )
			#print "vCFFAtoms[%d].atom_name is : %s"  % (len(vCFFAtoms),vCFFAtoms[len(vCFFAtoms)-1].atom_name)
			#print vCFFAtoms[len(vCFFAtoms)-1].R
		elif( isHydrophobe and ( ( distance_nearest >= distance_hydrophobic_contact_from ) and ( distance_nearest < distance_hydrophobic_contact_to ) ) ):
			vCFFAtoms.append(Atom( "6  ", "C", "GRI", "  1", grid[i].R ) )
			#print "vCFFAtoms[%d].atom_name is : %s"  % (len(vCFFAtoms),vCFFAtoms[len(vCFFAtoms)-1].atom_name)
			#print vCFFAtoms[len(vCFFAtoms)-1].R
		distance_nearest = 100
			
	#print vCFFAtoms[5].atom_name
	#print len(vCFFAtoms)
	'''
	for i in range(len(vCFFAtoms)):
		print "i is :",i
		print vCFFAtoms[i].R
	'''
	return vCFFAtoms
	

def getUnitVectorFromIToJ( atom_i, atom_j ):
	length = getDistance( atom_i.R, atom_j.R )
	unit_vector = [0.0,0.0,0.0]
	unit_vector[0] = ( atom_j.R[0] - atom_i.R[0] )/length
	unit_vector[1] = ( atom_j.R[1] - atom_i.R[1] )/length
	unit_vector[2] = ( atom_j.R[2] - atom_i.R[2] )/length
	return unit_vector

def  setResidueKinds(top_file):
#Read residues names from CHARMM topology file
	line = ""
	number = 0
	mHash_Residue_names = list()
	fileinput = open(top_file,"r")
	lineall = fileinput.read()
	fileinput.close()
	for line in lineall.splitlines():
		if "RESI" in line:
			mHash_Residue_names.append(line[5:8])
			#mHash_Residue_names.put( line.substring( 5, 8 ), new Integer(number)
			number += 1
		if END in line:
			break
		
def setAtomTypeAndChargesUsingCharmmTopologyFile (top_file, rec_coord) :
#Read atom names and types from CHARMM topology file
	line = ""
	total_assigned_atoms = 0 
	###### top_file ########
	fileinput = open(top_file,"r")
	charge_dict = dict()
	atom_type_dict = dict()
	
	flag = True
	while True:
		top_line = fileinput.readline()
		if (top_line.strip().startswith("END")):
			break
		if top_line.startswith("RESI"):
			resi_name = top_line[5:8].strip().upper()
			charge_dict[resi_name] = dict()
			atom_type_dict[resi_name] = dict()
			while True:
				top_line = fileinput.readline()
				if top_line.startswith ("ATOM"):
					#atom_name = ""
					atom_name = top_line[5:8].strip().upper()
					charge_dict[resi_name][atom_name] = float(top_line[17:22].strip())
					atom_type_dict[resi_name][atom_name] = top_line[10:14].strip().upper()
				if top_line.strip() == "":
					break
			
	fileinput.close()
	for i in range(len(rec_coord)):
		residue = rec_coord[i].residue_name.strip()
		atom = rec_coord[i].atom_name.strip()
		if not residue in charge_dict or not atom in charge_dict[residue]:
			sys.stdout.write( "ERROR: " + " " + rec_coord[i].atom_name + " (" + rec_coord[i].residue_name + rec_coord[i].sequence_number +") has been not assigned." )
			continue
		else:
			rec_coord[i].atom_type = atom_type_dict[residue][atom]
			rec_coord[i].charge =charge_dict[residue][atom]
	########################
def getProbes (par_file):
#Initializing Probes
	sys.stdout.write( "\n> Initializing probe elements" )
	vProbes = list()
	line = ""
	fileinput = open(par_file,"r")
	line = fileinput.readline() #Skip comments
	line = fileinput.readline()
	while line :
		if( line.startswith( "*" ) ):
			break
		else:
			probe = Atom()
			probe.atom_type = line[0:4]
			probe.charge = float(line[10:15].strip(" "))
			vProbes.append( probe )
			line = fileinput.readline()
	fileinput.close()
	sys.stdout.write( "Total number of probes: " + len(vProbes))
	return vProbes

def getProbesForExplicitHydrogens (par_file):
	sys.stdout.write( "\n> Initializing probe elements for explicit hydrogen atoms" )
	line = ""
	fileinput = open(par_file,"r")
	line = fileinput.readline() #Skip comments
	line = fileinput.readline()
	vProbes = list()
	while line:
		if( line.startswith( "*" ) ):
			break
		else:
			probe = Atom()
			probe.atom_name = "H   "
			probe.atom_type = line[0:4]
			probe.residue_name = line[12:15]
			probe.atom_name = line[23:26]
			probe.charge = float(line[36:42].strip(" ") )	
			vProbes.append( probe )
			line = fileinput.readline()
	fileinput.close()
	sys.stdout.write( "Total number of probes: " + len(vProbes))
	return vProbes

