#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string.h>

#include "translate_rotate.h"
#include "mol.h"
using namespace std;
void ShowMyvector(const vector<float> &cc)
{
   for(int i=0; i<cc.size(); ++i)
   {
      cout<<cc[i]<<'\t';
   }
   cout<<endl;
}


inline float deg2rad(float deg)
{
   //if(deg<0)
   //    return 2*PI-deg*PI/180;

   return deg*PI/180;
}

void translate_rotate(DOCKMol & mol) {
    vector<float> axisA(3, 0);
    vector<float> axisB(3, 0);
    vector<float> centers(3, 0);
    //centers[0]=centers[1]=centers[2]=0.0;
    for(int i;i<mol.num_atoms;i++)
    {
        centers[0]+=mol.x[i];
        centers[1]+=mol.y[i];
        centers[2]+=mol.z[i];
    }
    centers[0]=centers[0]/mol.num_atoms;
    centers[1]=centers[1]/mol.num_atoms;
    centers[2]=centers[2]/mol.num_atoms;
    
    axisA[0] = 1.0+centers[0]; axisA[1] = 1.0+centers[1]; axisA[2] = 1.0+centers[2];
    axisB[0] = 2.0+centers[0]; axisB[1] = 2.0+centers[1]; axisB[2] = 2.0+centers[2];
    float angle = 30.0;
    
    GroupRotation(axisA, axisB, angle, mol);
    
}
bool GroupRotation(const vector<float> &axisA, const vector<float> &axisB, float angle, DOCKMol & mol){
    if(axisA.size()!=3 || axisB.size()!=3)
   {
      cout<<"Error in GroupRotation()"<<endl;
      return false;
   }
   //cout<<"GroupRotation "<<angle<<endl;
   vector<float> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   //cout<<"deg2rad(angle): "<<deg2rad(angle)<<endl;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;  
   
   //test
   //for(int i=0; i<3; i++)
	//for(int j=0; j<3; j++)
	    //cout<<rotmtx[i][j]<<endl;
   float point_A[3];
   for(int i=0; i<mol.num_atoms; i++)
   {
      point_A[0]=mol.x[i]-axisA[0];
      point_A[1]=mol.y[i]-axisA[1];
      point_A[2]=mol.z[i]-axisA[2];
      
      //pointB[i]=axisA;
      mol.x[i]=axisA[0];
      mol.y[i]=axisA[1];
      mol.z[i]=axisA[2];
      mol.x[i]+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
      mol.y[i]+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
      mol.z[i]+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];
      /*
      point_A[0]=pointB[i][0]-axisA[0];
      point_A[1]=pointB[i][1]-axisA[1];
      point_A[2]=pointB[i][2]-axisA[2];
      point_A[3]=1;
      
      if(!MatrixTimesTransVector(rotmtx, point_A, pointB[i]))
      { 
         return false;
      }
      pointB[i].pop_back();
   
      pointB[i][0]+=axisA[0];
      pointB[i][1]+=axisA[1];
      pointB[i][2]+=axisA[2];
      */
   }
   return true;

}
/******************************************************************************
功能：用于大量空间点坐标的旋转变换
参数：axisA-axisB两个点构成空间旋转轴，angle是旋转角度，
     pointB是原来的空间点坐标，运算完毕后更新为新的坐标。
作者: CaoYang
时间：2008.6.14.
******************************************************************************/
bool GroupRotation1(const vector<float> &axisA, const vector<float> &axisB, float angle, vector<vector<float> > &pointB)
{
   if(axisA.size()!=3 || axisB.size()!=3)
   {
      cout<<"Error in GroupRotation()"<<endl;
      return false;
   }

   vector<float> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   cout<<"deg2rad(angle): "<<deg2rad(angle)<<endl;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
   //vector<float> point_A(4, 0);
   float point_A[3];
   for(int i=0; i<pointB.size(); ++i)
   {
      point_A[0]=pointB[i][0]-axisA[0];
      point_A[1]=pointB[i][1]-axisA[1];
      point_A[2]=pointB[i][2]-axisA[2];

      pointB[i]=axisA;
      pointB[i][0]+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
      pointB[i][1]+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
      pointB[i][2]+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];
      /*
      point_A[0]=pointB[i][0]-axisA[0];
      point_A[1]=pointB[i][1]-axisA[1];
      point_A[2]=pointB[i][2]-axisA[2];
      point_A[3]=1;
      
      if(!MatrixTimesTransVector(rotmtx, point_A, pointB[i]))
      { 
         return false;
      }
      pointB[i].pop_back();
   
      pointB[i][0]+=axisA[0];
      pointB[i][1]+=axisA[1];
      pointB[i][2]+=axisA[2];
      */
   }
   return true;
}

/************************Rotation Matrix****************************
* Give the Rotation Axis and the Rotation Angle
* Generate the Rotation Matrix
* Author: C.Y.
* Date 2006.12.6.
*************************End***************************************/
inline bool RotationMatrixB(const vector<float> &axis, float angle, matrix &romtx)
{
   if(axis.size()!=3) 
   {
      cout<<"Error 1 in RotationMatrixB()"<<endl;
      return false;
   }
   vector<float> ouc(3, 0);
   if(!norm(axis, ouc)) 
   {
      cout<<"Error 2 in RotationMatrixB()"<<endl;
      return false;
   }
   float c=cos(angle);
   float s=sin(angle);
   float t=1-c;
   SetMatrix(romtx, 4, 4);
   romtx[0][0]=t*ouc[0]*ouc[0]+c;romtx[0][1]=t*ouc[0]*ouc[1]+s*ouc[2];
   romtx[0][2]=t*ouc[0]*ouc[2]-s*ouc[1];romtx[0][3]=0;
   
   romtx[1][0]=t*ouc[1]*ouc[0]-s*ouc[2];romtx[1][1]=t*ouc[1]*ouc[1]+c;
   romtx[1][2]=t*ouc[1]*ouc[2]+s*ouc[0];romtx[1][3]=0;
   
   romtx[2][0]=t*ouc[2]*ouc[0]+s*ouc[1];romtx[2][1]=t*ouc[2]*ouc[1]-s*ouc[0];
   romtx[2][2]=t*ouc[2]*ouc[2]+c;romtx[2][3]=0;
   
   romtx[3][0]=0;romtx[3][1]=0;romtx[3][2]=0;romtx[3][3]=1;
   
   return true;
}



//Method to initialize a matrix
inline void SetMatrix(matrix &sm, int m, int n)
{
   vector<float> tmp(n, 0);
   sm.assign(m, tmp);
   //tmp.~vector<float>();
}

inline void subtract(const vector<float> &c1, const vector<float> &c2, vector<float> &cc)
{
   cc[0] = c1[0] - c2[0];
   cc[1] = c1[1] - c2[1];
   cc[2] = c1[2] - c2[2];
}

inline void crossproduct(const vector<float> &c1, const vector<float> &c2, vector<float> &cc)
{
   cc[0] = c1[1] * c2[2] - c1[2] * c2[1];
   cc[1] = c1[2] * c2[0] - c1[0] * c2[2];
   cc[2] = c1[0] * c2[1] - c2[0] * c1[1];
}

inline bool norm(const vector<float> &c, vector<float> &cc)
{
   float len = c[0]*c[0] + c[1]*c[1] +c[2]*c[2];
   if(len<Extra) 
   {
      cout<<"Error in norm()! length~=0!\n";
      ShowMyvector(c);
      ShowMyvector(cc);
      return false;
   }
   len = 1/sqrt(len);
   cc[0] = c[0]*len;
   cc[1] = c[1]*len;
   cc[2] = c[2]*len;
   return true;
}

inline float innerproduct(const vector<float> &c1, const vector<float> &c2)
{
   return c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2];
}


/************************Points to Dihedral************************/
//If use inline here there will be a error when combiling codes.
//Return the angle of c1-c2-c3-c4. Unit: radian
//points c1-c2-c3-c4 should be in two different pane.
//or there may be faults.
//by Yang Cao
//
//Change by chengxin: if they are on the same plane, return 2*PI
float Points2Dihedral(const vector<float> &c1, const vector<float> &c2, const vector<float> &c3, const vector<float> &c4)
{
   vector<float> vector1(3,0), vector2(3,0), vector3(3,0);
   subtract(c1, c2, vector1);
   subtract(c2, c3, vector2);
   subtract(c3, c4, vector3);
   
   vector<float> v1(3,0), v2(3,0);
   crossproduct(vector2, vector1, v1);
   crossproduct(vector3, vector2, v2);
   
   vector<float> v3(3,0), v4(3,0);
   if(!norm (v1, v3)) {cout<<"Error in Points2Dihedral 1\n"<<endl; return 2*PI;}//exit(1);}
   if(!norm (v2, v4)) {cout<<"Error in Points2Dihedral 2\n"<<endl; return 2*PI;}//exit(1);}
   
   float dihedral = innerproduct(v3, v4);
   
   if (dihedral>1 && dihedral<1+Extra)
   {
      //cout<<"dihedral "<<dihedral<<" "<<acos(dihedral)<<"\n";
      dihedral=1;
   }
   else if(dihedral<-1 && dihedral>-1-Extra)
   {
      dihedral=-1;
   }
   else if(dihedral>1+Extra || dihedral<-1-Extra)
   {
      cout<<"Error, float Points2Dihedral()\n";
      exit(0);
   }
   
   vector<float> v5(3,0);
   crossproduct(v4, v3, v5);
   float direction = innerproduct(v5, vector2);
   
   if (direction>0)
    {
       return  (acos(dihedral)*180.0)/PI;
    }  
    else
    {  
       return (-acos(dihedral)*180.0)/PI;
    }
   
}


bool GroupRotation(const vector<float> &axisA, const vector<float> &axisB, float angle, DOCKMol & mol, vector<bool> &flexibleatoms){
    if(axisA.size()!=3 || axisB.size()!=3)
   {
      cout<<"Error in GroupRotation()"<<endl;
      return false;
   }
   //cout<<"GroupRotation "<<angle<<endl;
   vector<float> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   //cout<<"deg2rad(angle): "<<deg2rad(angle)<<endl;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;  
   
   //test
   //for(int i=0; i<3; i++)
	//for(int j=0; j<3; j++)
	    //cout<<rotmtx[i][j]<<endl;
   float point_A[3];
   for(int i=0; i<mol.num_atoms; i++)
   {
      //cout<<flexibleatoms[i]<<" ";
      if(flexibleatoms[i]==0)
          continue;
      
      point_A[0]=mol.x[i]-axisA[0];
      point_A[1]=mol.y[i]-axisA[1];
      point_A[2]=mol.z[i]-axisA[2];
      
      //pointB[i]=axisA;
      mol.x[i]=axisA[0];
      mol.y[i]=axisA[1];
      mol.z[i]=axisA[2];
      mol.x[i]+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
      mol.y[i]+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
      mol.z[i]+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];
      /*
      point_A[0]=pointB[i][0]-axisA[0];
      point_A[1]=pointB[i][1]-axisA[1];
      point_A[2]=pointB[i][2]-axisA[2];
      point_A[3]=1;
      
      if(!MatrixTimesTransVector(rotmtx, point_A, pointB[i]))
      { 
         return false;
      }
      pointB[i].pop_back();
   
      pointB[i][0]+=axisA[0];
      pointB[i][1]+=axisA[1];
      pointB[i][2]+=axisA[2];
      */
   }
   return true;

}
