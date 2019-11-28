#!/usr/bin/env python
# 2016-1-6 wenyi Zhang
#GenerateBindingSiteNegativeImage
docstring='''
generate the negative iamge by the binding sites
cd the running file folder then run like following:
example: ./GenerateBindingSiteNegaviveImage.py receptor lig_center
inputfile: rec.pdb lig_charge.mol2 ligand_center.pdb(this file is the predicted binding sites)
'''

import sys,os,math,string,copy
from GenerateBindingSiteNegativeImage_MOD import *
#import pdb; pdb.set_trace()



#dir_target = dir_bin + sys.argv[1] + "/"
#rec_file = dir_target + "receptor.pdb"
#rec_file = dir_target + sys.argv[2]+".pdb"
rec_file = sys.argv[1] + ".pdb"
#predicted_sites_file = dir_target + "predicted_sites.pdb"
#predicted_sites_file = dir_target + "ligand_center.pdb"
#predicted_sites_file = dir_target + sys.argv[3] + ".pdb"
predicted_sites_file = sys.argv[2] + ".pdb"
#predicted_sites_num_lig_file = dir_target + "predicted_sites_num_ligands.txt"
#lig_file = dir_target + "lig.sdf"
#lig_file = dir_target + "lig_charge.mol2"
lig_file = "lig_charge.mol2"

mBoxCentroid = [0.0,0.0,0.0] # the coordinates of centroid of the box
mGridBoxSize = [0.0,0.0,0.0] 

mGridSize = 2.0  
mClashCutoffDist = 2.5
dist_incre = 1.0
num_ni_dissected = 3
#For enclosure calculation
mRay_length = 5.0
mEnclosure_ratio = 0.5    # default: 0.5



def set_grid_box (mBoxCentroid,mGridBoxSize):
# SET GRID BOX
#8 * 3 array, mGridBox0 = (a+x, b-y, c+z), mGridBox1 = (a+x, b-y, c-z), mGridBox2 = (a+x, b+y, c-z), mGridBox3 = (a+x, b+y, c+z), mGridBox4 = (a-x, b-y, c+z), mGridBox5 = (a-x, b-y, c-z), mGridBox6 = (a-x, b+y, c-z), mGridBox7 = (a-x, b+y, c+z),
	mGridBox = [ [ mBoxCentroid[0] + mGridBoxSize[0], mBoxCentroid[1] - mGridBoxSize[1], mBoxCentroid[2] + mGridBoxSize[2] ],[ mBoxCentroid[0] + mGridBoxSize[0], mBoxCentroid[1] - mGridBoxSize[1], mBoxCentroid[2] - mGridBoxSize[2]], [mBoxCentroid[0] + mGridBoxSize[0], mBoxCentroid[1] + mGridBoxSize[1], mBoxCentroid[2] - mGridBoxSize[2]], [mBoxCentroid[0] + mGridBoxSize[0], mBoxCentroid[1] + mGridBoxSize[1], mBoxCentroid[2] + mGridBoxSize[2]], [mBoxCentroid[0] - mGridBoxSize[0], mBoxCentroid[1] - mGridBoxSize[1], mBoxCentroid[2] + mGridBoxSize[2]], [mBoxCentroid[0] - mGridBoxSize[0], mBoxCentroid[1] - mGridBoxSize[1], mBoxCentroid[2] - mGridBoxSize[2]], [mBoxCentroid[0] - mGridBoxSize[0], mBoxCentroid[1] + mGridBoxSize[1], mBoxCentroid[2] - mGridBoxSize[2]], [mBoxCentroid[0] - mGridBoxSize[0], mBoxCentroid[1] + mGridBoxSize[1], mBoxCentroid[2] + mGridBoxSize[2]] ]
	return mGridBox

def set_grids (mGridBoxSize, mGridBox):	
#set grids
	grid_num_x = int( ( mGridBoxSize[0] * 2 ) / mGridSize )
	grid_num_y = int( ( mGridBoxSize[1] * 2 ) / mGridSize )
	grid_num_z = int( ( mGridBoxSize[2] * 2 ) / mGridSize )
	vGrids = list()
	probe_coord = [ mGridBox[0][0] - mGridSize/2.0, mGridBox[0][1] + mGridSize/2.0, mGridBox[0][2] - mGridSize/2.0 ]
	current_probe_coord = [0.0,0.0,0.0]
	for i in range (grid_num_x):
		current_probe_coord[0] = probe_coord[0] - mGridSize * i
		for j in range (grid_num_y):
			current_probe_coord[1] = probe_coord[1] + mGridSize * j
			for k in range (grid_num_z):
				current_probe_coord[2] = probe_coord[2] - mGridSize * k
				vGrids.append(Atom( "C  ", "C", "GRI", "  1", current_probe_coord ))	
	return vGrids

def zeros(N,M):
	'''make an NxM matrix filled with 0'''
	return [[0. for i in range(M)] for j in range(N)]

def filtering_grids (vGrids,vReceptor_filered, predicted_index,mBoxCentroid):
#FILTERING GRIDS
	#print "I am here"
	vTemp_Grids = list()
	#FILTERING 1: Removing grids too close to the receptor atoms
	mNumGrids = len(vGrids)	
	mNumReceptorAtoms = len(vReceptor_filered) 
	#print "mNumReceptorAtoms ",mNumReceptorAtoms
	for  i in range (mNumGrids-1, -1, -1):
		for j in range(0,mNumReceptorAtoms):
			distance = getDistance( vGrids[i].R, vReceptor_filered[j].R )
			if (distance < mClashCutoffDist ):
				del vGrids[i]
				break
	print "filter1: vGrid=",len(vGrids)
	# FILTERING 2: Removing grids too far from the receptor atoms
	array_mini_distance = [1000] * len(vGrids)
	cutoff_mini_distance = 4.5
	mNumGrids = len(vGrids)
	vTemp_Grids =  list()
	#print "mNumReceptorAtoms",mNumReceptorAtoms
	for i in range(mNumGrids):
		for j in range(mNumReceptorAtoms):
			distance =getDistance( vGrids[i].R, vReceptor_filered[j].R );
			#print i,distance
			if distance < array_mini_distance[i]:
				array_mini_distance[i] = float(distance)
				#print "array_mini_distance,",array_mini_distance[i]
	mNumGrids = len(vGrids)
		
	for i in range(mNumGrids):
		if( array_mini_distance[i] <= cutoff_mini_distance ):
			vTemp_Grids.append( vGrids[i] )
	vGrids = list()
	mNumGrids = len(vTemp_Grids)
	for i in range(mNumGrids):
		vGrids.append( vTemp_Grids[i] )
	vTemp_Grids = list()
	print "filter2: vGrid=",len(vGrids)
		
	# FILTERING 3: Removing separated grids from the main grid cluster
	cutoff_distance = math.sqrt( mGridSize*mGridSize + mGridSize*mGridSize + mGridSize*mGridSize )
	cutoff_num_grids = 3
	array_num_neighbors = [0] * mNumGrids
	for i in range (mNumGrids):
		for j in range(i+1,mNumGrids):
			distance = getDistance( vGrids[i].R, vGrids[j].R )
			if( distance <= cutoff_distance ):
				array_num_neighbors[i] += 1 
				array_num_neighbors[j] += 1
		
	vTemp_Grids = list()
	for i in range (mNumGrids):
		if( array_num_neighbors[i] >= cutoff_num_grids ):
			vTemp_Grids.append( vGrids[i] )
	mNumGrids = len(vTemp_Grids)
	vGrids = list()
		
	for i in range (mNumGrids):
		vGrids.append( vTemp_Grids[i] )
	print "filter3: vGrid=",len(vGrids)

	vTemp_Grids = list()
	#FILTERING 4: Calculating enclosure: removing highly solvent-exposed points by
	# Rotation and Coordinates setting
	mRx = zeros(3,3); mRy = zeros(3,3); mRz = zeros(3,3)
	#mRx = [[0.0] * 3 ] * 3; mRy = [[0.0] * 3 ] * 3; mRz = [[0.0] * 3 ] * 3
	mNumRays = 146; # ==2+8*18
	mRays = zeros(mNumRays,3)
	#mRays = [[0.0] * 3 ] * mNumRays
	temp = [0.0] * 3; R = [0.0] * 3; r = [0.0] * 3 
	ray_index = 0; R_dot_r = 0.0; d = 0.0 
	numStrike = 0 
	mNumGrids = len(vGrids)
	mNumReceptorAtoms = len(vReceptor_filered)		
	#print mNumGrids, mNumReceptorAtoms
	for i in range(mNumGrids-1,-1,-1):
		numStrike = 0
		atom = vGrids[i]
		#print atom.R
		mRays[0][0] = atom.R[0]; mRays[0][1] = atom.R[1]; mRays[0][2] = atom.R[2] + mRay_length
		mRays[1][0] = atom.R[0]; mRays[1][1] = atom.R[1]; mRays[1][2] = atom.R[2] - mRay_length
		temp[0] = mRays[0][0]; temp[1] = mRays[0][1]; temp[2] = mRays[0][2];
		horizotalTranslation( temp, atom.R )
		ray_index = 1
		for j in range(1,9):
			setRotationMatrix( mRx, mRy, mRz, 0, 20, 0 )
			setRotatedCoordinates( mRx, mRy, mRz, temp )
			for k in range(18):
				setRotationMatrix( mRx, mRy, mRz, 0, 0, 20 )
				setRotatedCoordinates( mRx, mRy, mRz, temp )
				ray_index += 1 
				mRays[ray_index][0] = temp[0] + atom.R[0];
				mRays[ray_index][1] = temp[1] + atom.R[1];
				mRays[ray_index][2] = temp[2] + atom.R[2];
		# mRays is a sphere around atom.R with radius being mRay_length
		#writeXYZinPDBFormat (mRays, open("mrays.pdb",'w'))

		#Calculating the fraction of radial rays that strike the receptor surface
		# for each of the 146 mRays points, draw a cylinder whose radius is 2.0. the center
		# axis of the cylinder is the vector from atom.R to mRays point. for each cylinder,
		# judge if there is any vReceptor_filtered point inside it.
		for j in range(mNumRays):
			R[0] = mRays[j][0] - atom.R[0]
			R[1] = mRays[j][1] - atom.R[1]
			R[2] = mRays[j][2] - atom.R[2]
			# R is a sphere whose center is (0,0,0) and radius is mRay_length
			for k in range(mNumReceptorAtoms):
				if( getDistance( atom.R, vReceptor_filered[k].R ) <= mRay_length ):
					r[0] = vReceptor_filered[k].R[0] - atom.R[0];r[1] = vReceptor_filered[k].R[1] - atom.R[1];r[2] = vReceptor_filered[k].R[2] - atom.R[2];
					R_dot_r = R[0]*r[0] + R[1]*r[1] + R[2]*r[2]
					if( R_dot_r >= 0 ):
						d = math.sqrt( ( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] ) - R_dot_r * R_dot_r / ( R[0]*R[0] + R[1]*R[1] + R[2]*R[2] ) )
						#print d
						if( d < 2.0 ):
							numStrike += 1
							break
		if( 1.*numStrike/mNumRays < mEnclosure_ratio ):
			#print "yes"
			del vGrids[i]
			mNumGrids = mNumGrids -1 
		
	print "filter4: vGrid=",len(vGrids)
	#print len(vGrids)
	#Dissecting Binding Sites 
	vNI_dissected = [[],[],[]]
	mNumGrids = len(vGrids)
	
	dist_init = longest_distance - dist_incre # radius-of-gyration of ligand - 1
	writeflag =  False
	#num_ni_dissected is 3 
	for s in range(num_ni_dissected):
		longest_distance_cutoff = dist_init + s * dist_incre
		for i in range(mNumGrids):
			#selected_atom = vGrids[i]
			distance = getDistance( vGrids[i].R, mBoxCentroid )
			#print distance
			if( distance <= longest_distance_cutoff ):
				vNI_dissected[s].append( vGrids[i] )
	print "grid num:" 
	print len(vNI_dissected[2])
	outputfilename = "negativeimage_model" + sys.argv[2][-1] + ".pdb"
	#outputfilename = "negativeimage_model" + ".pdb"
	if( len(vNI_dissected[2]) >=1): 
		#negativeimagepdbout = open(dir_target + "negativeimage_model.pdb","w") # outputfile
		#outputfilename=dir_target + "negativeimage_model" + sys.argv[3][-1] + "pdb"
		#outputfilename = "negativeimage_model" + ".pdb"
		negativeimagepdbout = open(outputfilename,"w")
		writeNegativeimageInPDBFormat(vNI_dissected[2],negativeimagepdbout)
		negativeimagepdbout.close()
		'''
		ligandclusterindex = int(predicted_index)-1
		print "ligandclusterindex is :",ligandclusterindex
		# write the template ligands
		template_ligands = open(dir_target + "negativeimage_ligandcluster.txt","w")
		ligandscluster = open(dir_target + "clusterbindingsite_ligands.txt","r")
		ligandsclusterreader = ligandscluster.read()
		ligandscluster.close()
		writeligand = False
		clusterarray = ligandsclusterreader.split("Cluster")
		template_ligands.write(str('\n'.join(clusterarray[predicted_index].splitlines()[1:]))+"\n")
		'''
		writeflag = True
	if ( len(vNI_dissected[2])==0):
		negativeimagepdbout = open(outputfilename,"w")
		#writeNegativeimageInPDBFormat(mBoxCentroid,negativeimagepdbout)
		negativeimagepdbout.write( "ATOM  "+'    '+'1'+" DUM  PBS  "+ \
			'    '+'1'+'    ')
		X = '%.3f' % float(mBoxCentroid[0])
		Y = '%.3f' % float(mBoxCentroid[1])
		Z = '%.3f' % float(mBoxCentroid[2])
		negativeimagepdbout.write(" " * (8-len(X)) + X + " " * (8-len(Y))+ Y + " " * (8-len(Z)) + Z + '\n')
		negativeimagepdbout.write("TER" + "\n")
		writeflag =  True
		
	return writeflag
	

if __name__ == "__main__":
# GET NEGATIVE IMAGES FOR EACH PREDICTED SITES
	#print sys.argv[1]
	mNumNI = 0
	vReceptor_filered = list()
	#print "haha"
	#print rec_file
	vReceptor = getHeavyAtomsFromFile(rec_file) 
	print "vREceptor num,",len(vReceptor)	
	#Binding site file: .pdb
	list_inputstream = open(predicted_sites_file,"r")
	#print predicted_sites_file
	#Number of ligands file: .txt
	
	#list_num_lig_inputstream = open(predicted_sites_num_lig_file, "r")
	#print predicted_sites_num_lig_file
	#out_shape = open(dir_target + "target_shape.sdf", "w")
	#out_cff = open( dir_target + "target_cff.sdf", "w")
	#longest_distance = getLongestDistanceFromCentroidForSDFFile( lig_file )
	longest_distance = getLongestDistanceFromCentroidForMol2File (lig_file)
	#print longest_distance
	
	mBoxSize = longest_distance + dist_incre 
	mGridBoxSize[0] = float(mBoxSize) 
	mGridBoxSize[1] = float(mBoxSize) 
	mGridBoxSize[2] = float(mBoxSize)
	#print  mGridBoxSize
	list_line = list_inputstream.readline()
	predicted_flag = False
	predicted_index = 0
	while list_line:
		#print list_line
		if (not list_line.strip()):
			continue
		#list_num_lig_line = list_num_lig_inputstream.readline()
		
		predicted_index += 1
		#array_num_lig = list_num_lig_line.split("\t")
		#num_ligands = int (array_num_lig[1].strip(" "))
		# step 1 BOX CENTROID
		
		mBoxCentroid[0] = float(list_line[30:38].strip(" "))
		#print list_line[38:46]
		mBoxCentroid[1] = float(list_line[38:46].strip(" "))
		#print list_line[46:54]
		mBoxCentroid[2] = float(list_line[46:54].strip(" "))
		#print list_line[46:54]
		####### SET GRID BOX ######
		mGridBox = set_grid_box (mBoxCentroid, mGridBoxSize)
		#print mGridBox
		###### set grids ######
		# the number of grids in x, y, and z axes
		vGrids = list()
		vGrids = set_grids (mGridBoxSize, mGridBox)	
			
		# step 2 FILTER RECEPTOR COORDINATES
		vReceptor_filered = list() ######### 
		mNumReceptorAtoms = len(vReceptor)
		mCutoff_receptor_atoms = mBoxSize + mRay_length
		
		#print "mCutoff_receptor_atoms",mCutoff_receptor_atoms
		#print "mNumReceptorAtoms,",mNumReceptorAtoms
		print mBoxCentroid
		for i in range(mNumReceptorAtoms):
			#print i,"--",mNumReceptorAtoms
			selected_atom = vReceptor[i];
			#print selected_atom.R
			distance = getDistance( selected_atom.R, mBoxCentroid )
			
			#print mCutoff_receptor_atoms
			if( distance <= mCutoff_receptor_atoms ):
				#print distance,"------",mCutoff_receptor_atoms
				vReceptor_filered.append( selected_atom )
				#print "step2,",len(vReceptor_filered);
		
		
		#step 3 FILTERING GRIDS
		writeflag = filtering_grids(vGrids,vReceptor_filered,predicted_index, mBoxCentroid)
		
		
		if(writeflag == True):
			#print predicted_index
			break	
		if(writeflag == False):
			#print predicted_index
			list_line = list_inputstream.readline()	
			
	list_inputstream.close()
	#list_num_lig_inputstream.close()
	
	


