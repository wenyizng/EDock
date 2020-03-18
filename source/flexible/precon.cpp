#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string.h>

#include "precon.h"
#include "mol.h"

bool readconformations(const char *filein, vector<vector<float> > &x,vector<vector<float> > &y,vector<vector<float> > &z, vector<float> &energy,DOCKMol &mol)
{

    std :: ifstream ifs(filein, ios::in);
    string line;
    /********************************/
    int natoms,nbonds,nresidues;
    //DOCKMol mol;
    char  tmp1[100],tmp2[100],tmp3[100],tmp4[100],subst_name[100];
    //vector<vector<float> > x,y,z;
    
    vector <float> tmp;
    float tempenergy;
    //char tmpstr1[100],tmpstr2[100];
    int molindex;
    /********************************/

    bool rflag=0;
    bool firstmol=1;
    

    while(getline(ifs,line))
    {
	if(line.substr(0,10)=="##########")
	{
            rflag=1;
            //cout<<line<<endl;
            sscanf(line.c_str(), "########## Name:\t%d\t%f",&molindex,&tempenergy);
            energy.push_back(tempenergy);
            
            //add the energy part
	    continue;
	    
	}    


	if(rflag==1 && firstmol==1)
	{
	    if(line.substr(0,17)=="@<TRIPOS>MOLECULE")
	    {
                
                getline(ifs,line);//CP6.pdb
                mol.title=line;
		getline(ifs,line);//30 31 1 0 0
                mol.mol_info_line = line;
		sscanf(line.c_str(), "%d %d %d", &natoms, &nbonds, &nresidues);
                mol.allocate_arrays(natoms, nbonds, nresidues);
                //the first the x y z setting
                tmp.assign(mol.num_atoms,0.0);
		x.push_back(tmp);
		y.push_back(tmp);
                z.push_back(tmp);
                //mol informaitons
                getline(ifs,line);
                mol.comment1 = line;
                getline(ifs,line);
                mol.comment2 = line;
                getline(ifs,line);
                mol.energy=line;
                getline(ifs,line);
                mol.comment3=line;
                continue;
            }    
	    if(line.substr(0,14)=="@<TRIPOS>ATOM")
            {
                for (int i = 0; i < mol.num_atoms; i++) 
                {
                    if (!getline(ifs, line)) 
                    {
                        cout <<
                        "ERROR:  Atom information missing from ligand file.  Program will terminate."
                        << endl;
                        return false;
                    }
                    sscanf(line.c_str(), "%s %s %f %f %f %s %s %s %f", tmp4,tmp1, &mol.x[i], &mol.y[i],
                           &mol.z[i], tmp2, tmp3,subst_name, &mol.charges[i]);
                    mol.atom_names[i] = tmp1;
                    mol.atom_types[i] = tmp2;
	            mol.atom_number[i] = tmp4;
	            mol.atom_residue_numbers[i] = tmp3;
                    mol.subst_names[i] = subst_name;
                    x[x.size()-1][i]=mol.x[i];
                    y[y.size()-1][i]=mol.y[i];
                    z[z.size()-1][i]=mol.z[i];
                }
                /*
                for(int j=0;j<mol.num_atoms;j++)
                {
                     cout<< "ATOM  " << setw(5) << j+1 << " " << setw(4)<<mol.atom_names[j] << " " << setw(3)<<"ALA" <<"  "<<setw(4) << j+1 <<"   " ;
	             cout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<mol.x[j]<<setw(8)<<mol.y[j]<<setw(8)<<mol.z[j]<<"\n";    
                }
                */
                for (;;) 
                {
                    if (!getline(ifs,line)) 
                    {
                         cout <<
                         "ERROR: @<TRIPOS>BOND indicator missing from ligand file.  Program will terminate."
                         << endl;
                         return false;
                    }

                   if (line.substr(0,13) == "@<TRIPOS>BOND")
                        break;
                }
                
                for (int i = 0; i < mol.num_bonds; i++) 
                {
                    if (!getline(ifs,line)) 
                    {
                        cout <<
                        "ERROR: Bond information missing from ligand file.  Program will terminate."
                        << endl;
                        return false;
                    }
                    sscanf(line.c_str(), "%*d %d %d %s", &mol.bonds_origin_atom[i],
                           &mol.bonds_target_atom[i], tmp1);

                   // adjust bond atom #'s to start at 0
                   mol.bonds_origin_atom[i]--; //
	
	           //cout <<mol.bonds_origin_atom[i]<<endl;
                   mol.bonds_target_atom[i]--;
                   mol.bond_types[i] = tmp1;
                }
                mol.id_ring_atoms_bonds();
                firstmol=0;
		rflag=0;
                //cout<<"the first mol end"<<endl;
                continue;               
	    } //end of read the first mol
        }
        if(rflag==1)//from the second
        {
            if(line.substr(0,14)=="@<TRIPOS>ATOM")
            {
                
                x.push_back(tmp);
		y.push_back(tmp);
                z.push_back(tmp);
                for (int i = 0; i < mol.num_atoms; i++) 
                {
                    if (!getline(ifs, line)) 
                    {
                        cout <<
                        "ERROR:  Atom information missing from ligand file.  Program will terminate."
                        << endl;
                        return false;
                    }
                    sscanf(line.c_str(), "%s %s %f %f %f %s %s %s ", tmp4,tmp1, &x[x.size()-1][i], 
                           &y[y.size()-1][i],&z[z.size()-1][i], tmp2, tmp3,subst_name);
                    
               }
               rflag=0;
           }         
       }	
		
    }
/*
    for(int i=0;i<x.size();i++)
    {
        cout<<i<<"***************"<<endl;
        for(int j=0;j<mol.num_atoms;j++)
        {
            cout<<x[i][j]<<"\t"<<y[i][j]<<"\t"<<z[i][j]<<endl;
        }
    }
*/
   
    return true;
}
