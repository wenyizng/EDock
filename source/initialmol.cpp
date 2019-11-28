#include <sstream> 
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <random>

#include "energy.h"
#include "MC_simulation.h"
#include "mol.h"
#include "translate_rotate.h"
#include "output.h"
#include "generaterandnum.h"
bool intialmol2(vector<vector<float> > &ini_x, vector<vector<float> > &ini_y, vector<vector<float> > &ini_z, int num_molvec,DOCKMol mol,float spherecenterx,float spherecentery,float spherecenterz,int initialnum,mt19937 &gen)
{
    //string output="initial_con.pdb";
    //ofstream ini_out(output.c_str());
    //if( num_molvec ==0)
    {
       DOCKMol tempmol;
       copy_molecule(tempmol,mol);
       cout<<"I am molvec.size()==0"<<endl;
       float ligcenterx,ligcentery,ligcenterz;
       ligcenterx=ligcentery=ligcenterz=0;
       for(int i=0;i<mol.num_atoms;i++)
       {
           ligcenterx +=mol.x[i];
           ligcentery +=mol.y[i];
           ligcenterz +=mol.z[i];
       }
       ligcenterx=ligcenterx/mol.num_atoms;
       ligcentery=ligcentery/mol.num_atoms;
       ligcenterz=ligcenterz/mol.num_atoms;
       cout<<"mol center"<<ligcenterx<<" "<<ligcentery<<" "<<ligcenterz<<" "<<endl;
       double t[3];
       t[0]=spherecenterx-ligcenterx;
       t[1]=spherecentery-ligcentery;
       t[2]=spherecenterz-ligcenterz;
       
       for(int atom=0;atom<mol.num_atoms;atom++)
       {
           mol.x[atom] +=t[0];
           mol.y[atom] +=t[1];
           mol.z[atom] +=t[2];
           
       }
       
       for(int atom=0;atom<mol.num_atoms;atom++)
       {
            ini_x[0][atom] = mol.x[atom];
            ini_y[0][atom] = mol.y[atom];
            ini_z[0][atom] = mol.z[atom];
            tempmol.x[atom] = mol.x[atom];
            tempmol.y[atom] = mol.y[atom];
            tempmol.z[atom] = mol.z[atom]; 
            
       }
       //outpdb(tempmol,ini_out,0);
       vector<float> axis1(3, 0);
       //axis1[0]=ligcenterx; axis1[1]=ligcentery; axis1[2]=ligcenterz;
       axis1[0]=spherecenterx; axis1[1]=spherecentery; axis1[2]=spherecenterz;
       for(int i=1;i<initialnum-num_molvec;i++)
       {
            double a,b,c;
            
            vector<float> axis2(3, 0);
            
            //cout<<axis1[0]<<"\t"<<axis1[1]<<"\t"<<axis1[2]<<endl;
            a = rand0to1(gen);
            b = rand0to1(gen);
            //cout<<"a"<<a<<" "<<"b"<<b<<" "<<"c"<<c<<endl;
            double theta = 2 * 3.1415926 * a;
            double phi = acos(2*b-1.0);
    
            axis2[0] = sin(phi) * cos(theta) + ligcenterx;
            axis2[1] = sin(phi) * sin(theta) + ligcentery;
            axis2[2]= cos(phi) + ligcenterz;
            //c=randgeneration(gen);
            //float angle = c*30.0;
            float angle = (360.0/(initialnum*1.0))*i;
            //cout<<angle<<endl;
            //cout<<"axis2:"<<axis2[0]<<" "<<axis2[1]<<" "<<axis2[2]<<endl;
            //float r=(axis2[0]-axis1[0])*(axis2[0]-axis1[0]) + (axis2[1]-axis1[1])*(axis2[1]-axis1[1])+ (axis2[2]-axis1[2])*(axis2[2]-axis1[2]);

            GroupRotation(axis2, axis1,angle,mol); // random rotate the mol   
            vector <float>(axis2).swap(axis2);
            
            for(int atom=0;atom<mol.num_atoms;atom++)
            {
                ini_x[i][atom] = mol.x[atom];
                ini_y[i][atom] = mol.y[atom];
                ini_z[i][atom] = mol.z[atom];
                tempmol.x[atom] = mol.x[atom];
                tempmol.y[atom] = mol.y[atom];
                tempmol.z[atom] = mol.z[atom]; 
            }
            //outpdb(tempmol,ini_out,i);
           
       } 
       //ini_out.close();
       /*
       //MC part 
       vector<DOCKMol> MCrotatemols;
       vector <float> tmpenergy;
       vector<vector<float> > tmpx,tmpy,tmpz;
       vector<float> replica_MC_energy;
       if(MC(amber, mol, c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,1,tmpx,tmpy,tmpz,tmpenergy,500.0,1.0,gen,replica_MC_energy,torsions,flexible_flag))
       {
           cout<<"MC can get how many conformaitons for the REMC!"<<endl;
           cout<<tmpx.size()<<endl;
           bool MC_flag=0;
           if(tmpx.size()>=initialnum && MC_flag==0)
           {
	       for(int i=0;i<initialnum;i++)
               {
                   DOCKMol copymol;
                   copy_molecule(copymol, mol);
                   
                   for(int atomi=0;atomi<mol.num_atoms;atomi++)
                   {
                       copymol.x[atomi]=tmpx[i][atomi];
                       copymol.y[atomi]=tmpy[i][atomi];   
                       copymol.z[atomi]=tmpz[i][atomi]; 
                   }
                   copymol.current_score=tmpenergy[i];   
                   molvec.push_back(copymol);
               }
               MC_flag=1;
           }
           if(tmpx.size()<initialnum && MC_flag==0)
           {
	       for(int i=0;i<tmpx.size();i++)
               {
                   DOCKMol copymol;
                   copy_molecule(copymol, mol);                  
                   for(int atomi=0;atomi<mol.num_atoms;atomi++)
                   {
                       copymol.x[atomi]=tmpx[i][atomi];
                       copymol.y[atomi]=tmpy[i][atomi];   
                       copymol.z[atomi]=tmpz[i][atomi]; 
                   } 
                   copymol.current_score=tmpenergy[i];   
                   molvec.push_back(copymol);
              }
	      for (int i=0;i<initialnum-tmpx.size();i++)
              {
                  DOCKMol copymol;
                  copy_molecule(copymol, mol); 
                  for(int atomi=0;atomi<mol.num_atoms;atomi++)
                  {
                      copymol.x[atomi]=tmpx[i%tmpx.size()][atomi];
                      copymol.y[atomi]=tmpy[i%tmpx.size()][atomi];
                      copymol.z[atomi]=tmpz[i%tmpx.size()][atomi];
                  }
                  copymol.current_score=tmpenergy[i%tmpx.size()];  
                  molvec.push_back(copymol);
              }	       

               MC_flag=1;  
          }
          if(tmpx.size()==0 && MC_flag==0)
          {
               cout<<"the initial conformation is bad, can't get the docking result!"<<endl;
               //finial output mol
               for(int i=0;i<initialnum;i++)
                   molvec.push_back(mol);
               //return false;
               return true;
          }
    }
    */
    return true;
    }
    
}
