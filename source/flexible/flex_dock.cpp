 #include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string.h>

#include "flex_dock.h"
#include "mol.h"
#include "translate_rotate.h"

using namespace std;
/******************************************************/
void
id_torsions(DOCKMol & mol, FLOATVec & vertex, vector<TORSION> &torsions)
{
    int             i,
                    j,
                    max_central,
                    central;
    int             a1,
                    a2,
                    a3,
                    a4;
    int             nbr;
    vector < int   >nbrs;
    TORSION         tmp_torsion;
    
    BREADTH_SEARCH  bfs;
    
    //torsions.clear();

    // loop over bonds- add flex bonds to torsion list

    for (i = 0; i < mol.num_bonds; i++) {
        if (mol.bond_active_flags[i]) {
            if (mol.bond_is_rotor(i)) {
                torsions.push_back(tmp_torsion);
                torsions[torsions.size() - 1].bond_num = i;
                vertex.push_back(0.0);
            }
        }
    }
    //vector<vector <int> > Latoms(torsions.size());//add by wenyizng
    //vector<vector <int> > Ratoms(torsions.size());//add by wenyizng
    //cout<<"torsions.size()"<<torsions.size()<<endl;
    //cout<<"Latoms.size()"<<Latoms.size()<<endl;
    //cout<<"Ratoms.size()"<<Ratoms.size()<<endl;
    
    // ID the inter-segment rot-bonds
    for (i = 0; i < torsions.size(); i++) {

        a2 = mol.bonds_origin_atom[torsions[i].bond_num];
        a3 = mol.bonds_target_atom[torsions[i].bond_num];

        max_central = -1;
        // nbrs = mol.get_atom_neighbors(a2);
        nbrs = mol.neighbor_list[a2];
	//cout<<"nbrs.size()"<<nbrs.size()<<endl;
        for (j = 0; j < nbrs.size(); j++) {
            nbr = nbrs[j];
		
            if (nbr != a3) {
                vector<int> tmp_Latoms;
                central = bfs.get_search_radius(mol, nbr, a2,tmp_Latoms);//central is int
                
                if (central > max_central) {
                    a1 = nbr;
                    max_central = central;
                    //cout<<"left:"<<central<<"tmp_Latoms.size()"<<tmp_Latoms.size()<<endl;
                     //because the a4 has been changed,so should remove the atoms in Latoms[i] firstly, and then add the new atoms
                    torsions[i].Latoms.clear();
                    torsions[i].Latoms.push_back(a1);
                    for(int m=0;m<tmp_Latoms.size();m++)
                    {
                        torsions[i].Latoms.push_back(tmp_Latoms[m]);
                        //cout<<tmp_Latoms[m]<<" ";
                    }
                    //cout<<"\n";
                        
                        
                }
            }
        }

        max_central = -1;
        // nbrs = mol.get_atom_neighbors(a3);
        nbrs = mol.neighbor_list[a3];

        for (j = 0; j < nbrs.size(); j++) {
            nbr = nbrs[j];
            if (nbr != a2) {
                vector<int> tmp_Ratoms;
                central = bfs.get_search_radius(mol, nbr, a3,tmp_Ratoms);
                if (central > max_central) {
                    a4 = nbr;
                    max_central = central;
                    //cout<<"right:"<<central<<"tmp_Ratoms.size()"<<tmp_Ratoms.size()<<endl;
                    //because the a4 has been changed,so should remove the atoms in Ratoms[i] firstly, and then add the new atoms
                    torsions[i].Ratoms.clear();
                    torsions[i].Ratoms.push_back(a4);
                    for(int m=0;m<tmp_Ratoms.size();m++)
                    {
                        torsions[i].Ratoms.push_back(tmp_Ratoms[m]);
                        //cout<<tmp_Ratoms[m]<<" ";
                    }
                    //cout<<"\n";
                        
                }
            }
        }

        torsions[i].atom1 = a1;
        torsions[i].atom2 = a2;
        torsions[i].atom3 = a3;
        torsions[i].atom4 = a4;
	//cout<<"1:"<<torsions[i].atom1<<"cor "<<mol.x[a1]<<","<<mol.y[a1]<<","<<mol.z[a1]<<" 2:"<<torsions[i].atom2<<" 3:"<<torsions[i].atom3<<" 4:"<<torsions[i].atom4<<" bond_num:"<<torsions[i].bond_num<<endl;
        //cout<<torsions[i].bond_num<<" bond:"<<torsions[i].atom1<<" "<<torsions[i].atom2<<" "<<torsions[i].atom3<<" "<<torsions[i].atom4<<" "<<endl;
        //cout<<"left atoms:"<<endl;
        //for(int m=0;m<torsions[i].Latoms.size();m++)
            //cout<<torsions[i].Latoms[m]<<" ";
        //cout<<endl;
        //cout<<"right atoms:"<<endl;
        //for(int m=0;m<torsions[i].Ratoms.size();m++)
           // cout<<torsions[i].Ratoms[m]<<" ";
        //cout<<endl;
        
        //compulate the 2Dihedral for a1 a2 a3 a4
        vector<float> c1(3,0);
        vector<float> c2(3,0);
        vector<float> c3(3,0);
        vector<float> c4(3,0);
        c1[0]=mol.x[torsions[i].atom1]; c1[1]=mol.y[torsions[i].atom1]; c1[2]=mol.z[torsions[i].atom1];
        c2[0]=mol.x[torsions[i].atom2]; c2[1]=mol.y[torsions[i].atom2]; c2[2]=mol.z[torsions[i].atom2];
        c3[0]=mol.x[torsions[i].atom3]; c3[1]=mol.y[torsions[i].atom3]; c3[2]=mol.z[torsions[i].atom3];
        c4[0]=mol.x[torsions[i].atom4]; c4[1]=mol.y[torsions[i].atom4]; c4[2]=mol.z[torsions[i].atom4];
        //cout<<"2dihedralis:"<<Points2Dihedral(c1,c2,c3,c4)<<endl;
        torsions[i].ihedral=Points2Dihedral(c1,c2,c3,c4);
	if((torsions[i].ihedral-(-180.0))<=1e-4)
		torsions[i].ihedral=180.0;
        c1.clear();c2.clear();c3.clear();c4.clear();
       
    }

}

vector<string> split(string str, string pattern)
{
	vector<string> ret;
	if (pattern.empty()) return ret;
	size_t start = 0, index = str.find_first_of(pattern, 0);
	while (index != str.npos)
	{	
		if (start != index)
		{
			ret.push_back(str.substr(start, index - start));
		}
			
			
		start = index + 1;
		index = str.find_first_of(pattern, start);
		
	}
	
	if (!str.substr(start).empty())
		ret.push_back(str.substr(start));
	return ret;
}

