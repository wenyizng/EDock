#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string.h>

#include "MathTools.h"
#include "mol.h"

/* read the ligand informaiton about the title, mol_info_line, comment1, comment2, comment3,     
    score_text_data,energy, simplex_text, mol_data, num_atoms, num_bonds, num_residues, and so on.
    and compulate the ring information*/
bool readlig(DOCKMol & mol, const char *file){
 
    char            line[1000];
    int             count;
    int             i;
    int             n1;
    bool            atom_line;
    char            typ[100],
                    col[100];

    char            tmp1[100],
                    tmp2[100],
		    jwu_tmp4[100],
		    jwu_tmp3[100],
                    subst_name[100];
    int             natoms,
                    nbonds,
                    nresidues; //defined by jwu
    std :: string   l1,
                    l2,
                    l3,
                    l4,
                    l5,
                    l6;
    float           f1,
                    f2,
                    f3,
                    f4,
                    f5;

    atom_line = false;
    std :: ifstream ifs(file, ios::in); //read file
    // read forward until the tripos molecule tag is reached
    for (;;) {
	
        if (! ifs.getline(line, 1000)) {
            cout << endl << "ERROR: Ligand file empty.  Program will terminate." << endl;
            return false;
        }

        if (!strncmp(line, "@<TRIPOS>MOLECULE", 17))
            break;
         
    }
    for (count = 0;; count++) {

        if (!ifs.getline(line, 1000)) {
            cout << endl <<
                "ERROR:  Ligand file empty.  Program will terminate." << endl;
            return false;
        }

        if (!strncmp(line, "@<TRIPOS>ATOM", 13)) {
            atom_line = true;
            break;
        }
        // assign the first 5 header lines to the proper fields
        switch (count) {

        case 0:
            l1 = line;
            break;

        case 1:
            l2 = line;
            break;

        case 2:
            l3 = line;
            break;

        case 3:
            l4 = line;
            break;

        case 4:
            l5 = line;
            break;

        case 5:
            l6 = line;
            break;
        }

    }
    if (!atom_line) {
        cout <<
            "ERROR: @<TRIPOS>ATOM indicator missing from ligand file.  Program will terminate."
            << endl;
        return false;
    }
    sscanf(l2.c_str(), "%d %d %d", &natoms, &nbonds, &nresidues);
    mol.allocate_arrays(natoms, nbonds, nresidues);

    mol.title = l1;
    mol.mol_info_line = l2;
    mol.comment1 = l3;
    mol.comment2 = l4;
    mol.energy = l5;
    mol.comment3 = l6;
    for (i = 0; i < mol.num_atoms; i++) {
        if (!ifs.getline(line, 1000)) {
            cout <<
                "ERROR:  Atom information missing from ligand file.  Program will terminate."
                << endl;
            return false;
        }
// 1     N       15.2586  -59.3416   35.3528 N.4     1 PRO1  -0.2020

        sscanf(line, "%s %s %f %f %f %s %s %s %f", jwu_tmp4,tmp1, &mol.x[i], &mol.y[i],
               &mol.z[i], tmp2, jwu_tmp3,subst_name, &mol.charges[i]);
	//jwu_tmp4:atom number, tmp1:atom name, x,y,z,tmp2,subst_name,
        mol.atom_names[i] = tmp1;
  //jwu
        mol.atom_types[i] = tmp2;
	mol.atom_number[i] = jwu_tmp4;
	mol.atom_residue_numbers[i] = jwu_tmp3;
        mol.subst_names[i] = subst_name;

    }
   for (;;) {
        if (!ifs.getline(line, 1000)) {
            cout <<
                "ERROR: @<TRIPOS>BOND indicator missing from ligand file.  Program will terminate."
                << endl;
            return false;
        }

        if (!strncmp(line, "@<TRIPOS>BOND", 13))
            break;
    }
  for (i = 0; i < mol.num_bonds; i++) {
        if (!ifs.getline(line, 1000)) {
            cout <<
                "ERROR: Bond information missing from ligand file.  Program will terminate."
                << endl;
            return false;
        }
       sscanf(line, "%*d %d %d %s", &mol.bonds_origin_atom[i],
               &mol.bonds_target_atom[i], tmp1);

        // adjust bond atom #'s to start at 0
        mol.bonds_origin_atom[i]--; //
	
	//cout <<mol.bonds_origin_atom[i]<<endl;
        mol.bonds_target_atom[i]--;
        mol.bond_types[i] = tmp1;

    }
  mol.id_ring_atoms_bonds();
 

  return true;
} 

void read_binding_sites(const char *file, vector<vector <int> > &binding_sites)
{
    
    
    ifstream fp(file, ios::in);
    string line;
    if (!fp)
    {
        cout << endl <<
            "ERROR:  coach_clsuter file empty.  Program will terminate." << endl;
        exit(1);
    }

    while(!fp.eof())
    {

        getline(fp,line);
        
        if(line.substr(0,4) == "Site")
        {
            //cout<<line<<"zwy"<<endl;
            vector <int> tmpint;
            binding_sites.push_back(tmpint);
            getline(fp,line);
            //cout<<line<<"(((((((((((((((((("<<endl;
            char *cstr = new char[line.length() + 1];
            strcpy(cstr, line.c_str());
            const char *sep = ","; 
            char *p;
            p = strtok(cstr, sep);
            while(p)
           {
               //cout<<"site "<<atoi(p)<<endl;
               binding_sites[binding_sites.size()-1].push_back(atoi(p));
               //cout<<"binding_sites.size()"<<binding_sites.size()<<endl;
               p = strtok(NULL, sep);
           }
           //cout<<"\n";
           delete [] cstr;
        }

        
    }
}

void read_dat(const char *file, vector<vector <float> > &x,vector<vector <float> > &y, vector<vector <float> > &z,vector <vector <int> > &atomindex)
{
    
    
    //std :: ifstream ifs(file, ios::in); //read file
    ifstream fp(file, ios::in);
    string line;
    if (!fp)
    {
        cout << endl <<"ERROR:  "<<file<<" file empty." << endl;
        exit(0);
    }

    while(!fp.eof())
    {

        getline(fp,line);
        //cout<<line;
        
        if(line.substr(0,4) == "####")
        {
            //cout<<line<<"zwy"<<endl;
            vector <float> tempf;
            vector <int> bindsite_temp;
            x.push_back(tempf); y.push_back(tempf); z.push_back(tempf);
            atomindex.push_back(bindsite_temp);
            while(1)
            {
               getline(fp,line);
               
               if(line.substr(0,3) == "TRE")
                break;
               //cout<<line<<endl;
               int residue;
               float xt,yt,zt;
               sscanf(line.c_str(), "%d %f %f %f", &residue,&xt,&yt,&zt);
               //cout<<residue<<"&&"<<xt<<"&&"<<yt<<"&&"<<zt<<endl;
               x[x.size()-1].push_back(xt);
               y[y.size()-1].push_back(yt);
               z[z.size()-1].push_back(zt);
               atomindex[atomindex.size()-1].push_back(residue);

            }   
        }
    }
    /*
    for(int i=0;i<atomindex.size();i++)
    {
        for(int j=0;j<atomindex[i].size();j++)
        {
            cout<<atomindex[i][j]<<"**"<<x[i][j]<<"**"<<y[i][j]<<"**"<<z[i][j]<<endl;
        }
        cout<<"TRE"<<endl;
    }
    */
}
