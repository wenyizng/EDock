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
#include "output.h"

void outlig(DOCKMol & mol){
    cout<<"****************test out******************"<<endl;
    cout<<"title is:\n"<<mol.title<<endl;
    cout<<"mol_info_line is: \n"<<mol.mol_info_line<<endl;
    cout<<"comment1 is: \n"<<mol.comment1<<endl;
    cout<<"comment2 is: \n"<<mol.comment2<<endl;
    cout<<"comment3 is: \n"<<mol.comment3<<endl;
    cout<<"score_text_data is: \n"<<mol.score_text_data<<endl;
    cout<<"energy is: \n"<<mol.energy<<endl;
    cout<<"simplex_text is: \n"<<mol.simplex_text<<endl;
    cout<<"mol_data is: \n"<<mol.mol_data<<endl;
    cout<<"num_atoms is: \n"<<mol.num_atoms<<endl;
    cout<<"num_bonds is: \n"<<mol.num_bonds<<endl;
    cout<<"num_residues is: \n"<<mol.num_residues<<endl;

    for (int i=0;i<mol.num_atoms;i++){
        cout <<mol.atom_names[i]<<"    "<<mol.atom_types[i]<<"    "<<mol.atom_number[i]<<"   "<<mol.atom_residue_numbers[i]<<"    "<<mol.subst_names[i]<<"    "<<mol.x[i]<<"    "<<mol.y[i]<<"     "<<
        mol.z[i]<<"    "<<mol.charges[i]<< "    " << mol.atom_active_flags[i] << endl;
    }

    cout<<"ring info:"<<endl;
    for(int m=0;m<mol.num_atoms;m++){
       cout<<mol.atom_ring_flags[m]<<" ";
    }
    cout<<"\n";
    for(int m=0;m<mol.num_bonds;m++){
       cout<<mol.bond_ring_flags[m]<<" ";
    }
    cout<<"\n";

    cout<<"****************test end******************"<<endl;

}

void outpdb(DOCKMol & mol, ofstream & fp,const int index){
    //ofstream fp(fpname);
    fp << "MODEL     "<<setw(4)<<index<<endl;
    //cout <<"mol.num_atoms=" <<mol.num_atoms <<endl;

    for(int i=0;i<mol.num_atoms;i++){
	//cout<<"i="<<i;
	//cout<<"\tatom_names[i]="<<mol.atom_names[i];
	//cout<<"\tx[i]="<<mol.x[i] <<endl;

	fp << "ATOM  " << setw(5) << i+1 << " " << setw(4)<<mol.atom_names[i] << " " << setw(3)<<"LIG" <<"  "<<setw(4) << i+1 <<"   " ;
        fp<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<mol.x[i]<<setw(8)<<mol.y[i]<<setw(8)<<mol.z[i]<<"\n";
    }
    fp<<"ENDMDL\n";
}

/***********************************************************************/
bool Write_Mol2(DOCKMol & mol, ofstream & ofs, string s)
{
    char            line[1000];
    int             i;
    vector < int   >renumber;
    int             current_atom,
                    current_bond;

    // init atom/bond renumbering data
    renumber.resize(mol.num_atoms + 1);
    current_atom = 1;
    current_bond = 1;

    // write out header information
    ofs << "@<TRIPOS>MOLECULE" << endl;

    if (mol.title.empty())
        ofs << "*****" << endl;
    
    else
    {
        //ofs << mol.title << endl;
        ofs<<s<<"\t"<<endl;  
    }
        
        

    sprintf(line, " %d %d 1 0 0", mol.num_active_atoms, mol.num_active_bonds);
    ofs << line << endl;

    ofs << mol.comment1 << endl;
    ofs << mol.comment2 << endl;
    ofs << mol.energy << endl;
    ofs << mol.comment3 << endl;

    // write out atom lines
    ofs << "@<TRIPOS>ATOM" << endl;

    for (i = 0; i < mol.num_atoms; i++) {
        if (mol.atom_active_flags[i]) {

            sprintf(line,
                    "%7d%1s%-6s%12.4f%10.4f%10.4f%1s%-5s%4s%1s %-8s%10.4f",
                    current_atom, "", mol.atom_names[i].c_str(), mol.x[i],
                    mol.y[i], mol.z[i], "", mol.atom_types[i].c_str(), "1", "",
                    mol.subst_names[i].c_str(), mol.charges[i]);

            ofs << line << endl;

            renumber[i] = current_atom;
            current_atom++;
        }
    }

    // write out bond lines
    ofs << "@<TRIPOS>BOND" << endl;

    for (i = 0; i < mol.num_bonds; i++) {
        if (mol.bond_active_flags[i]) {

            sprintf(line, "%6d%6d%6d%3s%2s", current_bond,
                    renumber[mol.bonds_origin_atom[i]],
                    renumber[mol.bonds_target_atom[i]], "",
                    mol.bond_types[i].c_str());

            ofs << line << endl;
            current_bond++;
        }
    }

    ofs << "@<TRIPOS>SUBSTRUCTURE" << endl;
    sprintf(line,
            "     1 %-4s        1 TEMP              0 ****  ****    0 ROOT",
            mol.subst_names[0].c_str());
    ofs << line << endl << endl;


    return true;
}

bool Write_Mol2_2(DOCKMol &mol, vector<float> &x,vector<float> &y, vector<float> &z,ofstream & ofs)
{
    char            line[1000];
    int             i;
    vector < int   >renumber;
    int             current_atom,
                    current_bond;

    // init atom/bond renumbering data
    renumber.resize(mol.num_atoms + 1);
    current_atom = 1;
    current_bond = 1;

    // write out header information
    ofs << "@<TRIPOS>MOLECULE" << endl;

    if (mol.title.empty())
        ofs << "*****" << endl;
    else
        ofs << mol.title << endl;

    sprintf(line, " %d %d 1 0 0", mol.num_active_atoms, mol.num_active_bonds);
    ofs << line << endl;

    ofs << mol.comment1 << endl;
    ofs << mol.comment2 << endl;
    ofs << mol.energy << endl;
    ofs << mol.comment3 << endl;

    // write out atom lines
    ofs << "@<TRIPOS>ATOM" << endl;

    for (i = 0; i < mol.num_atoms; i++) {
        if (mol.atom_active_flags[i]) {

            sprintf(line,
                    "%7d%1s%-6s%12.4f%10.4f%10.4f%1s%-5s%4s%1s %-8s%10.4f",
                    current_atom, "", mol.atom_names[i].c_str(), x[i],
                    y[i], z[i], "", mol.atom_types[i].c_str(), "1", "",
                    mol.subst_names[i].c_str(), mol.charges[i]);

            ofs << line << endl;

            renumber[i] = current_atom;
            current_atom++;
        }
    }

    // write out bond lines
    ofs << "@<TRIPOS>BOND" << endl;

    for (i = 0; i < mol.num_bonds; i++) {
        if (mol.bond_active_flags[i]) {

            sprintf(line, "%6d%6d%6d%3s%2s", current_bond,
                    renumber[mol.bonds_origin_atom[i]],
                    renumber[mol.bonds_target_atom[i]], "",
                    mol.bond_types[i].c_str());

            ofs << line << endl;
            current_bond++;
        }
    }

    ofs << "@<TRIPOS>SUBSTRUCTURE" << endl;
    sprintf(line,
            "     1 %-4s        1 TEMP              0 ****  ****    0 ROOT",
            mol.subst_names[0].c_str());
    ofs << line << endl << endl;


    return true;
}

