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

//bool MC(AMBER_TYPER & amber, DOCKMol &mol,Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,int cutoff,vector<vector <float> > &x,vector<vector <float> > &y,vector<vector <float> > &z,vector <float> &tmpenergy, float temprature,float angel,mt19937 &gen,vector<int> &MCsteps)
bool MC(AMBER_TYPER & amber, DOCKMol &mol,Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,float cutoff,vector<vector <float> > &x,vector<vector <float> > &y,vector<vector <float> > &z,vector <float> &tmpenergy, float temprature,mt19937 &gen,vector<float> &replica_MC_energy,vector<TORSION> &torsions, bool flexible_flag,vector<float> sphx, vector<float> sphy, vector<float> sphz)
{

    DOCKMol original,endmol;
    copy_molecule(original,mol); // input is mol, copy it to the original
    int maxsimulations=501; // max simulations times
    //cout<<"start:temprature is:"<<temprature<<"\t"<<"method is:"<<method<<"\t"<<mol.current_score<<endl;
    cout<<"start:temprature is:"<<temprature<<"\t"<<mol.current_score<<endl;
    int tcout=0;
    copy_molecule(endmol,mol);
    int totaltimes=10000;
    bool flexible_option;
    if(flexible_flag==0) 
        flexible_option=0;
    
    DOCKMol tmpmol;
    copy_molecule(tmpmol, mol);
    while(maxsimulations > 1)
    {
	
	//cout<<"temprotatemols.size()\t"<<temprotatemols.size()<<endl;
	//cout<<temprotatemols[temprotatemols.size()-1].current_score<<endl;
	//DOCKMol endmol;
	//cout<<"MC"<<maxsimulations<<endl;
	/**************************generatationmol**********************************/
        //add the flexible informaitons
        //option1, when iterations return 20, will do one time flexible docking
        //where do the flexible movement
        if(flexible_flag == 1 && maxsimulations%10==0 )
        {
        //flexible_option=1;
            for(int flex_num=0;flex_num<1;flex_num++) // 10 times flexible movement
            {
                swapmol(tmpmol,endmol);//copy the current mol to the tmpmol, and flexible tmpmol
                if(flexible_mol(tmpmol,gen,torsions)) // it can be flexible movement, if not flexible movement, then the endmol not change
                {
                     if(energy(amber, tmpmol,c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,sphx,sphy,sphz)) // sometime, the flexible mol can not be calulated the energy, so ther is a justify wether it can be calculated enenrgy.
                     {
                         swapmol(endmol,tmpmol);
                         //cout<<"mc flexible part:"<<endmol.internal_energy<<" "<<endmol.current_score<<endl;
                         //go to the MC to justy whether the flexible can be acceptable
                         int id = 1;		
	                 double detaenergy = endmol.current_score - original.current_score;
		         //cout<<"dateenergy is:"<<detaenergy<<endl;
                         if(detaenergy>0)
	                 {
	                     double r = rand0to1(gen);
                             double condation = exp((-detaenergy)/(temprature));
                             //cout<<r<<"\t"<<condation<<endl;
		             //cout<< "r and condation    "<<r<<"\t"<<condation<<endl;
		             if(r>condation)
		             {
			         id = 3;
			         //cout<<"before:"<<endmol.current_score<<endl;
			         swapmol(endmol,original);
                                 //replica_MC_energy.push_back(endmol.current_score);
			         //cout<<"after:"<<endmol.current_score<<endl;
		              }
		             //else
                             //{
                             //    cout<<r<<"\t"<<condation<<"\t"<<endmol.current_score<<endl;
                             //} 
	                 }
                         //accept
                         if(id == 1)
		         {
		             /**************add the endmol x y z coordate*************************/
		             vector <float> tmp(endmol.num_atoms,0.0);
		             x.push_back(tmp);
		             y.push_back(tmp);
                             z.push_back(tmp);
		             for(int i=0;i<endmol.num_atoms;i++)
		             {
			         x[x.size()-1][i]=endmol.x[i];
			         y[y.size()-1][i]=endmol.y[i];
			         z[z.size()-1][i]=endmol.z[i];
		             }
		             tmpenergy.push_back(endmol.current_score);
		             //int tempstep=20000-maxsimulations;
                             //MCsteps.push_back(tempstep);
		             swapmol(original, endmol);//put endmol in original
		             replica_MC_energy.push_back(endmol.current_score);
                             cout<<endmol.current_score<<endl;
		             /*
		             for(int j=0;j<temprotatemols[temprotatemols.size()-1].num_atoms;j++)
	                     {
	                        cout<< "ATOM  " << setw(5) << j+1 << " " << setw(4)<<temprotatemols[temprotatemols.size()-1].atom_names[j] << " " << setw(3)<<"ALA" <<"  "<<setw(4) << j+1 <<"   " ;
	                        cout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<temprotatemols[temprotatemols.size()-1].x[j]<<setw(8)<<temprotatemols[temprotatemols.size()-1].y[j]<<setw(8)<<temprotatemols[temprotatemols.size()-1].z[j]<<"\n";
	                     }
		             cout<<"\n"<<endl;
		            */	
		          }

                     //***********************************************************//
                     }   
                }   
            }    
         }
        int method = 1;
        if(maxsimulations%20==0)
	{
            method=2;
            //cout<<"^^^^^^^^^^^^^I will turn upside down"<<endl;
        }
            
	
	if(generatationmol(original,endmol,method,gen))
	{ 
            
            	 	
	    maxsimulations--;
	    totaltimes--;
	    if(totaltimes<0) // if 100*1000 times also can finish, break;
		break;
	    if(energy(amber, endmol,c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,sphx,sphy,sphz))
	    {
                
		int id = 1;
		
	        double detaenergy = endmol.current_score -
                                       original.current_score;
		//cout<<"dateenergy is:"<<detaenergy<<endl;
		
	        if(detaenergy>0)
	        {
	            double r = rand0to1(gen);
                    double condation = exp((-detaenergy)/(temprature));
                    //cout<<r<<"\t"<<condation<<endl;
		    //cout<< "r and condation    "<<r<<"\t"<<condation<<endl;
		    if(r>condation)
		    {
			id = 3;
			//cout<<"before:"<<endmol.current_score<<endl;
			swapmol(endmol,original);
                        replica_MC_energy.push_back(endmol.current_score);
			//cout<<"after:"<<endmol.current_score<<endl;
		    }
		    //else
                    //{
                    //    cout<<r<<"\t"<<condation<<"\t"<<endmol.current_score<<endl;
                    //} 
    		    
	        }
		//else if(detaenergy<0)
                //{
                //    cout<<"0.00000"<<"\t"<<"1.0000"<<"\t"<<endmol.current_score<<endl;
                //}
               
	    //accept
		if(id == 1)
		{
                    if(method==2)
                    {
                        cout<<"turn upside down "<<endmol.current_score<<"     accept"<<endl;
                    }
			
		    /**************add the endmol x y z coordate*************************/
		    vector <float> tmp(endmol.num_atoms,0.0);
		    x.push_back(tmp);
		    y.push_back(tmp);
                    z.push_back(tmp);
		    for(int i=0;i<endmol.num_atoms;i++)
		    {
			x[x.size()-1][i]=endmol.x[i];
			y[y.size()-1][i]=endmol.y[i];
			z[z.size()-1][i]=endmol.z[i];
		    }
		    tmpenergy.push_back(endmol.current_score);
		    //int tempstep=20000-maxsimulations;
                    //MCsteps.push_back(tempstep);
		    swapmol(original, endmol);//put endmol in original
		    replica_MC_energy.push_back(endmol.current_score);
		    /*
		    for(int j=0;j<temprotatemols[temprotatemols.size()-1].num_atoms;j++)
	            {
	                cout<< "ATOM  " << setw(5) << j+1 << " " << setw(4)<<temprotatemols[temprotatemols.size()-1].atom_names[j] << " " << setw(3)<<"ALA" <<"  "<<setw(4) << j+1 <<"   " ;
	                cout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<temprotatemols[temprotatemols.size()-1].x[j]<<setw(8)<<temprotatemols[temprotatemols.size()-1].y[j]<<setw(8)<<temprotatemols[temprotatemols.size()-1].z[j]<<"\n";
	            }
		    cout<<"\n"<<endl;
		    */	
		}
	    }
	    else
            {
                maxsimulations++;
            }
	}	
    }
    cout<<tmpenergy.size()<<endl;
    if(tmpenergy.size()>0)
        cout<<"end:acceptance "<<tmpenergy.size()<<"\t"<<tmpenergy.size()*1.0/500.0<<" final energy:"<<tmpenergy[tmpenergy.size()-1]<<endl;
    if(tmpenergy.size()==0)
   {cout<<"end:acceptance 0"<<" final energy:"<<mol.current_score<<"\n I want check the program, xsize is "<<x.size()<<endl;
	return false;
   }
    if(x.size() == 0)
	return false;
    //cout<<"temprotatemols.size():"<<temprotatemols.size()<<endl;
    //cout<<"tcout"<<"\t"<<tcout<<endl;
    else
    {
          return true;
    }
   
}


vector<int> sort_index(RotatemolVec temprotatemols)
{
    vector<int> original_index;
    vector<float> molscores;
    for(int i =0;i<temprotatemols.size()-1;i++)
    {
        original_index.push_back(i);
        molscores.push_back(temprotatemols[i].endmol.current_score);
    }
	


    for(int i=0;i<molscores.size();i++)
    {
	for(int j=i+1;j<molscores.size();j++)
	{
	     int temp;
	     float tempscore;
	     if(molscores[j] <molscores[i])
	     {
		 tempscore = molscores[j];
		 molscores[j] = molscores[i];
		 molscores[i] = tempscore;
	         temp=original_index[j];
		 original_index[j]=original_index[i];
		 original_index[i]=temp; 
	     }
	     
	}        
		
    }
    //for(int i=0;i<original_index.size();i++)
    //{
    //    cout<<i<<":"<<original_index[i]<<temprotatemols[original_index[i]].endmol.current_score<<endl;
    //}
    vector <float>(molscores).swap(molscores);
    return original_index;
}

bool generatationmol(DOCKMol & original,DOCKMol & endmol,int method,mt19937 &gen)
{
    
    swapmol(endmol,original); //put original in endmol;
    double a,b,c,t[3],angle;
    vector<float> centers(3, 0);
    for(int i=0;i<endmol.num_atoms;i++)
    {
        centers[0]+=endmol.x[i];
        centers[1]+=endmol.y[i];
        centers[2]+=endmol.z[i];
    }
    centers[0]=centers[0]/endmol.num_atoms;
    centers[1]=centers[1]/endmol.num_atoms;
    centers[2]=centers[2]/endmol.num_atoms;

    float tspacing, rspacing;
    if (method==1)
    {
        tspacing = 0.2; // translate spacing
        rspacing = 3.0;
        c=randgeneration(gen);
        angle = c* rspacing;
        
	
    }
    if (method==2)
    {
        tspacing = 0.0;
        rspacing = largeangle(gen);
        angle=rspacing;
	//cout<<"turn angle " <<angle<<endl;
    
    }
    
    //translation operator
    
    vector<float> axis1(3, 0);
    vector<float> axis2(3, 0);
    for(int i=0;i<3;i++)
    {
        //t[i] = ((rand()*1.0/(RAND_MAX*1.0))*2.0-1) * spacing;
        t[i] = randgeneration(gen) * tspacing;
        //cout<<t[i]<<endl;
    } 
    //cout<<t[0]<<" "<<t[1]<<" "<<t[2]<<endl;      
    for(int atom=0;atom<endmol.num_atoms;atom++)
    {
         endmol.x[atom] +=t[0];
         endmol.y[atom] +=t[1];
         endmol.z[atom] +=t[2];
    }
    
    //rotation operator
    axis1[0]=centers[0]; axis1[1]=centers[1]; axis1[2]=centers[2];
    a = rand0to1(gen);
    b = rand0to1(gen);
    double theta = 2 * 3.1415926 * a;
    double phi = acos(2*b-1.0);
    
    axis2[0] = sin(phi) * cos(theta) + centers[0];
    axis2[1] = sin(phi) * sin(theta) + centers[1];
    axis2[2]= cos(phi) + centers[2];
    
    //cout<<"angle is:"<<angle<<endl;
    //cout<<temprotatemols[0].startmol.current_score<<endl;
    GroupRotation(axis1, axis2, angle, endmol); // random rotate the mol    
    vector <float>(axis1).swap(axis1);
    vector <float>(axis2).swap(axis2);
    vector <float>(centers).swap(centers);
    
    //old 
    /*	
    for(int i=0;i<3;i++)
    {
        //a = (rand()*1.0/(RAND_MAX*1.0))*2.0-1;
        //b = (rand()*1.0/(RAND_MAX*1.0))*2.0-1;
        a=randgeneration(gen);
        b=randgeneration(gen);
        axis1[i] = a+centers[i];
        axis2[i] = b+centers[i];  
    }
    */
    //new
    return true;
}

bool flexible_mol(DOCKMol & mol,mt19937 &gen,vector<TORSION> &torsions)
{
    vector<bool> flexibleatoms(mol.num_atoms,0);
    
    int bondindex=randint(gen,torsions.size());
    //cout<<"bondindex is"<<bondindex<<endl;    
    for(int j=0;j<torsions[bondindex].Latoms.size();j++)   
    {   
        flexibleatoms[torsions[bondindex].Latoms[j]]=1;
    }
    vector<float> axisA(3,0);
    vector<float> axisB(3,0);
    float flexible_angle=10.0;
        
    float angle=flexible_angle*randgeneration(gen); //randnum from-1 to 1
    axisA[0]=mol.x[torsions[bondindex].atom2];axisA[1]=mol.y[torsions[bondindex].atom2];
    axisA[2]=mol.z[torsions[bondindex].atom2];
    axisB[0]=mol.x[torsions[bondindex].atom3];axisB[1]=mol.y[torsions[bondindex].atom3];
    axisB[2]=mol.z[torsions[bondindex].atom3];
    GroupRotation(axisA, axisB, angle, mol, flexibleatoms);
    vector <float>(axisA).swap(axisA);
    vector <float>(axisB).swap(axisB);
    vector <bool>(flexibleatoms).swap(flexibleatoms); 
    return true;
}

void swapmol(DOCKMol & mol1, DOCKMol & mol2)//put mol2 in mol1
{
    for(int i=0;i<mol1.num_atoms;i++)
    {
	mol1.x[i]=mol2.x[i];
	mol1.y[i]=mol2.y[i];
	mol1.z[i]=mol2.z[i];	
    }
    mol1.current_score=mol2.current_score;
    mol1.internal_energy=mol2.internal_energy;
    mol1.intral_energy=mol2.intral_energy;
}












