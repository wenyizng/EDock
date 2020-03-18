#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string.h>

#include "mol.h"
//#include "MathTools.hpp"
using namespace std;
/*********************************************/
bool
DOCKMol::atoms_are_one_three(int a1, int a2)
{
    vector < int   >nbrs1,
                    nbrs2;
    int             i,
                    j;

    nbrs1 = get_atom_neighbors(a1);
    nbrs2 = get_atom_neighbors(a2);

    for (i = 0; i < nbrs1.size(); i++) {
        for (j = 0; j < nbrs2.size(); j++) {
            if (nbrs1[i] == nbrs2[j])
                return true;
        }
    }

    return false;
}
/*********************************************/
bool
DOCKMol::atoms_are_one_four(int a1, int a2)
{
    vector < int   >nbrs1,
                    nbrs2;
    int             i,
                    j;

    nbrs1 = get_atom_neighbors(a1);
    nbrs2 = get_atom_neighbors(a2);

    for (i = 0; i < nbrs1.size(); i++) {
        for (j = 0; j < nbrs2.size(); j++) {
            if (get_bond(nbrs1[i], nbrs2[j]) != -1)
                return true;
        }
    }

    return false;
}
void
DOCKMol::allocate_arrays(int natoms, int nbonds, int nresidues) //added jwu
{
    int             i;

    //clear_molecule();

    num_atoms    = natoms;
    num_bonds    = nbonds;
    num_residues = nresidues; // number in current mol

//    num_residues_in_rec = 1000; // number in current rec


    atom_data = new string[num_atoms];

    x = new float[num_atoms];
    y = new float[num_atoms];
    z = new float[num_atoms];

    charges = new float[num_atoms];
    atom_types = new string[num_atoms];
// hbond
    flag_acceptor    = new bool[num_atoms];
    flag_donator     = new bool[num_atoms];
    acc_heavy_atomid = new int[num_atoms];

//jwu code
    atom_number=new string[num_atoms];
    atom_residue_numbers=new string[num_atoms];
    atom_names = new string[num_atoms];
    subst_names = new string[num_atoms];
    atom_color = new string[num_atoms]; // kxr205
    atom_psol = new float[num_atoms];   // kxr205
    atom_apsol = new float[num_atoms];  // kxr205

    // DTM - 11-12-08 - allocate the atom_segment_ids array
    atom_segment_ids = new int[num_atoms];

    // added for inter molecular H-Bonds in xlogp
    // commented out by sudipto and trent
    /*number_of_H_Donors = 0;
    number_of_H_Acceptors = 0;
    H_bond_donor = new int[num_atoms];
    H_bond_acceptor = new int[num_atoms]; */

    atom_ring_flags = new bool[num_atoms];
    bond_ring_flags = new bool[num_bonds];

    bonds_origin_atom = new int[num_bonds];
    bonds_target_atom = new int[num_bonds];
    bond_types = new string[num_bonds];

    atom_active_flags = new bool[num_atoms];
    bond_active_flags = new bool[num_bonds];

    num_active_atoms = num_atoms;
    num_active_bonds = num_bonds;

    amber_at_bump_id = new int[num_atoms];
    amber_at_heavy_flag = new int[num_atoms];
    amber_at_id = new int[num_atoms];
    amber_at_valence = new int[num_atoms];
    amber_at_radius = new float[num_atoms];
    amber_at_well_depth = new float[num_atoms];

    chem_types = new string[num_atoms];

    amber_bt_id = new int[num_bonds];
    amber_bt_minimize = new int[num_bonds];
    amber_bt_torsion_total = new int[num_bonds];
    amber_bt_torsions = new FLOATVec[num_bonds];

    gb_hawkins_radius = new float[num_atoms];
    gb_hawkins_scale = new float[num_atoms];
    
    // footprints.
    //footprints.clear(); // performed in clear_molecule() 

    neighbor_list = new INTVec[num_atoms];
    // ie_neighbor_list = new bool[num_atoms*num_atoms];

    atom_child_list = new INTVec[num_bonds * 2];

    arrays_allocated = true;

    // assign init values
    for (i = 0; i < num_atoms; i++) {

        atom_data[i] = "";
        atom_names[i] = "";

        //footprint
        atom_residue_numbers[i] ="";
        atom_names[i] = "";

        subst_names[i] = "";
        atom_types[i] = "";
        atom_color[i] = "null"; // kxr205
        atom_psol[i] = 0.0;     // kxr205
        atom_apsol[i] = 0.0;
        charges[i] = 0.0;

        // DTM - 11-12-08 - initialize the atom_segment_ids array
       	atom_segment_ids[i] = -1;

// hbond
        flag_acceptor[i]    = false;
        flag_donator[i]     = false;
        acc_heavy_atomid[i] = 0;

        atom_ring_flags[i] = false;
        atom_active_flags[i] = true;

        amber_at_bump_id[i] = 0;
        amber_at_heavy_flag[i] = 0;
        amber_at_id[i] = 0;
        amber_at_valence[i] = 0;
        amber_at_radius[i] = 0.0;
        amber_at_well_depth[i] = 0.0;

        chem_types[i] = "";

        gb_hawkins_radius[i] = 0.0;
        gb_hawkins_scale[i] = 0.0;

    }

    for (i = 0; i < num_bonds; i++) {

        bonds_origin_atom[i] = 0;
        bonds_target_atom[i] = 0;
        bond_types[i] = "";

        bond_ring_flags[i] = false;
        bond_active_flags[i] = true;

        amber_bt_id[i] = 0;
        amber_bt_minimize[i] = 0;
        amber_bt_torsion_total[i] = 0;
    }
    // test_child_list = new bool[2*num_bonds*num_atoms];

}

void
DOCKMol::id_ring_atoms_bonds()
{
    vector < bool > atoms_visited, bonds_visited;
    vector < int   >atom_path,
                    bond_path;
    int             i;

    atoms_visited.clear();
    bonds_visited.clear();

    atoms_visited.resize(num_atoms, false);
    bonds_visited.resize(num_bonds, false);
    
    atom_path.clear();
    bond_path.clear();

    for (i = 0; i < num_atoms; i++)
        atoms_visited[i] = false;

    for (i = 0; i < num_bonds; i++)
        bonds_visited[i] = false;

    // loop over all atoms and find rings
    
    for (i = 0; i < num_atoms; i++){
	//cout<<i<<":"<<atoms_visited[i]<<"   ";
        if (!atoms_visited[i]){
        
        find_rings(atom_path, bond_path, atoms_visited, bonds_visited, i);
	    //find every atom's ring info
         
	}
            
       
       }
    //cout<<"print_vector"<<endl;
    //print_vector(bond_path);

}

/*************************************/
void
DOCKMol::find_rings(vector < int >apath, vector < int >bpath,
                    vector < bool > &atoms, vector < bool > &bonds, int atnum)
{
    
    
    
    int             i,
                    j,
                    nbr_bond;
    vector < int   >nbrs;
    //cout<<atnum<<":"<<atoms[atnum]<<endl;
    
    if (atoms[atnum]) {

        i = apath.size() - 1;
        j = bpath.size() - 1;
        //cout<<"atom_path is :"<<apath.size()<<"****"<<"bond_path is:"<<bpath.size()<<endl;
        while ((i >= 0) && (j >= 0)) {  // ///// changed from > to >= fixes the 
                                        // first atom bug
            atom_ring_flags[apath[i--]] = true;
            bond_ring_flags[bpath[j--]] = true;

            if (i == -1)        // added fix to address if i = 0
                break;
      
            else if (apath[i] == atnum)
                break;

        }
 
    //cout<<"atom_ring_flags is"<<endl;
    //for (int m;m<num_atoms;m++)
       //cout<<m+1<<":"<<atom_ring_flags[m]<<" ";
    } else {

        atoms[atnum] = true;
        nbrs = get_atom_neighbors(atnum);

	//cout<<"print nbrs"<<endl;
	//cout<<atnum<<":";
	//print_vector(nbrs);
	//cout<<nbrs.size()<<endl;

        for (i = 0; i < nbrs.size(); i++) {
	    //cout << "neighber---"<<nbrs[i]<<":"<<atnum<<endl;
            nbr_bond = get_bond(nbrs[i], atnum);
	    
	    //cout<<"nbr_bond:"<<nbr_bond<<endl;
	    

            if (!bonds[nbr_bond]) {

                bonds[nbr_bond] = true;

                apath.push_back(nbrs[i]);
                bpath.push_back(nbr_bond);

                find_rings(apath, bpath, atoms, bonds, nbrs[i]);
                //cout<<"atom_path is :"<<apath.size()<<"****"<<"bond_path is:"<<bpath.size()<<endl;
                apath.pop_back();
                bpath.pop_back();
            }

        }

    }

}

vector < int >
DOCKMol::get_atom_neighbors(int index) //find the atom ID adjacent with index
{
    vector < int   >nbrs;
    int             i;

    nbrs.clear();

    // loop over bonds and find any that include the index atom
    for (i = 0; i < num_bonds; i++) {

        if (bonds_origin_atom[i] == index)
            nbrs.push_back(bonds_target_atom[i]);
	    //print_vector(nbrs);
        if (bonds_target_atom[i] == index)
            nbrs.push_back(bonds_origin_atom[i]);

    }
    //
    return nbrs;
}

int
DOCKMol::get_bond(int a1, int a2) // whether  a1---a2 has a bond, if have, return the bond index, otherwise, return -1 
{
    int             bond_id;
    int             i;

    bond_id = -1;

    for (i = 0; i < num_bonds; i++) {

        if ((bonds_origin_atom[i] == a1) && (bonds_target_atom[i] == a2))
            bond_id = i;

        if ((bonds_origin_atom[i] == a2) && (bonds_target_atom[i] == a1))
            bond_id = i;

    }

    return bond_id;
}


vector < int   >
DOCKMol::get_atom_children(int a1, int a2) // find the the atoms which is adjacent a2
{
    vector < int   >children;
    bool           *visited;
    vector < int   >nbrs,
                    new_nbrs;
    int             i;

    if (get_bond(a1, a2) != -1) {
        visited = new bool[num_atoms];
        memset(visited, 0, num_atoms * sizeof(bool));
        // visited.resize(num_atoms, false);

        visited[a1] = true;
        visited[a2] = true;

        nbrs = neighbor_list[a2];//a2 neighbors 
	//cout << "nbrs" <<endl;
	//print_vector(nbrs);
        while (nbrs.size() > 0) {

            if (!visited[nbrs[nbrs.size() - 1]]) { //from the last one

                children.push_back(nbrs[nbrs.size() - 1]);
                new_nbrs.clear();
                // new_nbrs = get_atom_neighbors(nbrs[nbrs.size()-1]);
                new_nbrs = neighbor_list[nbrs[nbrs.size() - 1]];
                visited[nbrs[nbrs.size() - 1]] = true;
                nbrs.pop_back();
                for (i = 0; i < new_nbrs.size(); i++)
                    nbrs.push_back(new_nbrs[i]);

            } else {
                nbrs.pop_back();
            }

        }

        delete[]visited;
    }
    //cout<<"children"<<endl;
    //print_vector(children);
    return children;
    
}


/*************************************/
void
DOCKMol::clear_molecule()
{
    int             i;

    amber_at_assigned = false;
    amber_bt_assigned = false;
    
    if (arrays_allocated) {
  
        //cout<<"arrays_allocated is:"<<arrays_allocated<<endl;
        //cout<<"num_bonds is:"<<num_bonds<<endl;
        //cout<<"amber_bt_torsions size is:"<<amber_bt_torsions[0].size()<<endl;
        for (i = 0; i < num_bonds; i++)
            amber_bt_torsions[i].clear();

        for (i = 0; i < num_atoms; i++)
            neighbor_list[i].clear();

        for (i = 0; i < 2 * num_bonds; i++)
            atom_child_list[i].clear();

        delete[]atom_data;
        delete[]x;
        delete[]y;
        delete[]z;

        // hbond
        delete[]flag_acceptor;
        delete[]flag_donator;
        delete[]acc_heavy_atomid;
  
        // Arrays obtained from within xlogp.cpp
        // Commented out since we are not defining this any more: sudipto and trent
        /*delete[]H_bond_donor;
        delete[]H_bond_acceptor; */

 
        delete[]charges;
        delete[]atom_types;
        delete[]atom_names;

        //footprint info.
        delete[]atom_number;
        delete[]atom_residue_numbers;
        delete[]subst_names; // residue name

        delete[]atom_color;     // kxr205
        delete[]atom_psol;      // kxr205
        delete[]atom_apsol;

	// DTM - 11-12-08 - delete the atom_segment_ids array
	delete[]atom_segment_ids;
	
        delete[]bonds_origin_atom;
        delete[]bonds_target_atom;
        delete[]bond_types;

        delete[]atom_ring_flags;
        delete[]bond_ring_flags;

        delete[]atom_active_flags;
        delete[]bond_active_flags;

        delete[]amber_at_id;
        delete[]amber_at_radius;
        delete[]amber_at_well_depth;
        delete[]amber_at_heavy_flag;
        delete[]amber_at_valence;
        delete[]amber_at_bump_id;

        delete[]chem_types;

        delete[]amber_bt_id;
        delete[]amber_bt_minimize;
        delete[]amber_bt_torsion_total;
        delete[]amber_bt_torsions;

        delete[]gb_hawkins_radius;
        delete[]gb_hawkins_scale;

        // footprints
        //footprints.clear(); 

        delete[]neighbor_list;
        // delete [] ie_neighbor_list;
        // delete [] test_child_list;
        delete[]atom_child_list;

        arrays_allocated = false;
    }
    // clear MOL2 header info
    title = "";
    mol_info_line = "";
    comment1 = "";
    comment2 = "";
    comment3 = "";
    simplex_text = "";
    energy = "";
    mol_data = "";

    // clear scalar data
    num_atoms = 0;
    num_bonds = 0;
    num_residues = 0;
    //num_residues_in_rec = 0;
    num_active_atoms = 0;
    num_active_bonds = 0;
    score_text_data = "";
    current_score = 0.0;
    internal_energy = 0.0;
    intral_energy = 0.0;
    current_data = "ERROR: Conformation could not be scored by DOCK.\n";
    primary_data = "ERROR: Conformation could not be scored by DOCK.\n";
    hbond_text_data = "";
    amber_score_ligand_id = "";
    grid_num = 0; // This defaults to the first grid read in
    total_dsol = 0.0;

}

/***********************************************************************/
void copy_molecule(DOCKMol & target, DOCKMol & original)
{
    int             i;
    
    target.allocate_arrays(original.num_atoms, original.num_bonds, original.num_residues);
    
    // copy scalar data
    target.title = original.title;
    target.mol_info_line = original.mol_info_line;
    target.comment1 = original.comment1;
    target.comment2 = original.comment2;
    target.comment3 = original.comment3;
    target.simplex_text = original.simplex_text;
    target.energy = original.energy;

    target.mol_data = original.mol_data;

    target.num_atoms = original.num_atoms;
    target.num_bonds = original.num_bonds;
    target.num_residues = original.num_residues; // this is the number in this molicule

    target.num_active_atoms = original.num_active_atoms;
    target.num_active_bonds = original.num_active_bonds;
    target.score_text_data = original.score_text_data;

    target.amber_at_assigned = original.amber_at_assigned;
    target.amber_bt_assigned = original.amber_bt_assigned;
    target.chem_types_assigned = original.chem_types_assigned;

    target.current_score = original.current_score;
    target.internal_energy = original.internal_energy;
    target.intral_energy = original.intral_energy;
    target.current_data = original.current_data;
    target.primary_data = original.primary_data;
    target.hbond_text_data = original.hbond_text_data;

    target.grid_num = original.grid_num;
    target.amber_score_ligand_id = original.amber_score_ligand_id;
    target.total_dsol = original.total_dsol;
    //  xlogp code : commented out sudipto and trent
    /*target.number_of_H_Donors = original.number_of_H_Donors;
    target.number_of_H_Acceptors = original.number_of_H_Acceptors; */
    
    // copy arrays
    for (i = 0; i < original.num_atoms; i++) {

        target.atom_data[i] = original.atom_data[i];

        target.x[i] = original.x[i];
        target.y[i] = original.y[i];
        target.z[i] = original.z[i];

        //hbond
        target.flag_acceptor[i] = original.flag_acceptor[i]; // this is a h-bond acceptor
        target.flag_donator[i]  = original.flag_donator[i]; // this is true for polar hydrogen
        target.acc_heavy_atomid[i] = original.acc_heavy_atomid[i]; // atom id of acceptor that is connected to polar h
                           // equal to zero if not polar hydrogen

        target.charges[i] = original.charges[i];
        target.atom_types[i] = original.atom_types[i];
        target.atom_names[i] = original.atom_names[i];

        //footprint info.
        target.atom_number[i] = original.atom_number[i];
        target.atom_residue_numbers[i] = original.atom_residue_numbers[i];
        target.subst_names[i] = original.subst_names[i];

        target.atom_color[i] = original.atom_color[i];  // kxr205
        target.atom_psol[i] = original.atom_psol[i];    // kxr205
        target.atom_apsol[i] = original.atom_apsol[i];  // kxr205

        //  xlogp
        // commmented out: sudipto and trent 
        /* target.H_bond_donor[i]    = original.H_bond_donor[i];
        target.H_bond_acceptor[i] = original.H_bond_acceptor[i]; */


        target.atom_ring_flags[i] = original.atom_ring_flags[i];
        target.atom_active_flags[i] = original.atom_active_flags[i];

        target.amber_at_id[i] = original.amber_at_id[i];
        target.amber_at_radius[i] = original.amber_at_radius[i];
        target.amber_at_well_depth[i] = original.amber_at_well_depth[i];
        target.amber_at_heavy_flag[i] = original.amber_at_heavy_flag[i];
        target.amber_at_valence[i] = original.amber_at_valence[i];
        target.amber_at_bump_id[i] = original.amber_at_bump_id[i];

        target.chem_types[i] = original.chem_types[i];

        target.gb_hawkins_radius[i] = original.gb_hawkins_radius[i];
        target.gb_hawkins_scale[i] = original.gb_hawkins_scale[i];

        target.neighbor_list[i] = original.neighbor_list[i];

        // DTM - 11-12-08 - copy the new atom_segment_ids array
        target.atom_segment_ids[i] = original.atom_segment_ids[i];
    }
    
    
    for (i = 0; i < original.num_bonds; i++) {

        target.bonds_origin_atom[i] = original.bonds_origin_atom[i];
        target.bonds_target_atom[i] = original.bonds_target_atom[i];
        target.bond_types[i] = original.bond_types[i];

        target.bond_ring_flags[i] = original.bond_ring_flags[i];
        target.bond_active_flags[i] = original.bond_active_flags[i];

        target.amber_bt_id[i] = original.amber_bt_id[i];
        target.amber_bt_minimize[i] = original.amber_bt_minimize[i];
        target.amber_bt_torsion_total[i] = original.amber_bt_torsion_total[i];
        target.amber_bt_torsions[i] = original.amber_bt_torsions[i];

    }
    /*
    // footprints
    target.footprints.clear();
    FOOTPRINT_ELEMENT temp_footprint_ele;
    int num_residues_in_rec = original.footprints.size();
    // loops over all residues and appends the values of the orginal mol 
    // on to the new mol
    for (i = 0; i < num_residues_in_rec; i++) {
        temp_footprint_ele.resname = original.footprints[i].resname;
        temp_footprint_ele.resid   = original.footprints[i].resid;
        temp_footprint_ele.vdw     = original.footprints[i].vdw;
        temp_footprint_ele.es      = original.footprints[i].es;
        temp_footprint_ele.hb      = original.footprints[i].hb;
        target.footprints.push_back(temp_footprint_ele);
    }
    */
    // test_child_list = new bool[2*num_bonds*num_atoms];

    // for(i=0;i<original.num_atoms*original.num_atoms;i++)
    // target.ie_neighbor_list[i] = original.ie_neighbor_list[i];

    for (i = 0; i < 2 * original.num_bonds; i++)
        target.atom_child_list[i] = original.atom_child_list[i];
    
}


/*********************************************/
void
DOCKMol::translate_mol(const DOCKVector & vec)
{
    int             i;

    for (i = 0; i < num_atoms; i++) {
        x[i] += vec.x;
        y[i] += vec.y;
        z[i] += vec.z;
    }

}

/*********************************************/
void
DOCKMol::rotate_mol(double mat[9])
{
    int             i;
    double          nx,
                    ny,
                    nz;

    for (i = 0; i < num_atoms; i++) {
        nx = x[i];
        ny = y[i];
        nz = z[i];

        x[i] = mat[0] * nx + mat[1] * ny + mat[2] * nz;
        y[i] = mat[3] * nx + mat[4] * ny + mat[5] * nz;
        z[i] = mat[6] * nx + mat[7] * ny + mat[8] * nz;
    }

}

/*********************************************/
void
DOCKMol::rotate_mol(double mat[3][3])
{
    double          new_mat[9];
    int             i,
                    j,
                    k;

    k = 0;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            new_mat[k++] = mat[i][j];

    rotate_mol(new_mat);

}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
DOCKVector::operator=(const DOCKVector & v)
{

    x = v.x;
    y = v.y;
    z = v.z;

    return *this;

}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::operator+=(const DOCKVector & v)
{

    x += v.x;
    y += v.y;
    z += v.z;

    return *this;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::operator+=(const float &f)
{

    x += f;
    y += f;
    z += f;

    return *this;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::operator-=(const DOCKVector & v)
{

    x -= v.x;
    y -= v.y;
    z -= v.z;

    return *this;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::operator-=(const float &f)
{

    x -= f;
    y -= f;
    z -= f;

    return *this;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::operator*=(const float &f)
{

    x *= f;
    y *= f;
    z *= f;

    return *this;
}

//bool   DOCKVector operator<(const DOCKMol &mol1, const DOCKMol &mol2){
//    return (mol1.current_score < mol2.current_score);
//}
// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::operator/=(const float &f)
{

    x /= f;
    y /= f;
    z /= f;

    return *this;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::normalize_vector()
{
    float
        len;

    len = length();

    if (len != 0) {
        x = x / len;
        y = y / len;
        z = z / len;
    }

    return *this;
}

// +++++++++++++++++++++++++++++++++++++++++
float
DOCKVector::length() const
{
    float           len;

    len = sqrt(x * x + y * y + z * z);

    return len;
}

// +++++++++++++++++++++++++++++++++++++++++
float
DOCKVector::squared_length() const
{
    float           len;

    len = x * x + y * y + z * z;

    return len;
}

// +++++++++++++++++++++++++++++++++++++++++
float
DOCKVector::squared_dist(const DOCKVector & v) const
{
    float           dist;

    dist =
        (x - v.x) * (x - v.x) + (y - v.y) * (y - v.y) + (z - v.z) * (z - v.z);

    return dist;
}

// +++++++++++++++++++++++++++++++++++++++++
int
operator==(const DOCKVector & v1, const DOCKVector & v2)
{

    if ((v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z))
        return (true);
    else
        return (false);

}

// +++++++++++++++++++++++++++++++++++++++++
int
operator!=(const DOCKVector & v1, const DOCKVector & v2)
{

    if ((v1.x != v2.x) || (v1.y != v2.y) || (v1.z != v2.z))
        return (true);
    else
        return (false);

}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator+(const DOCKVector & v1, const DOCKVector & v2)
{
    DOCKVector      vec;

    vec.x = v1.x + v2.x;
    vec.y = v1.y + v2.y;
    vec.z = v1.z + v2.z;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator+(const float &f, const DOCKVector & v2)
{
    DOCKVector      vec;

    vec.x = f + v2.x;
    vec.y = f + v2.y;
    vec.z = f + v2.z;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator+(const DOCKVector & v1, const float &f)
{
    DOCKVector      vec;

    vec.x = v1.x + f;
    vec.y = v1.y + f;
    vec.z = v1.z + f;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator-(const DOCKVector & v1, const DOCKVector & v2)
{
    DOCKVector      vec;

    vec.x = v1.x - v2.x;
    vec.y = v1.y - v2.y;
    vec.z = v1.z - v2.z;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator-(const DOCKVector & v)
{
    DOCKVector      vec;

    vec.x = -v.x;
    vec.y = -v.y;
    vec.z = -v.z;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator-(const DOCKVector & v1, const float &f)
{
    DOCKVector      vec;

    vec.x = v1.x - f;
    vec.y = v1.y - f;
    vec.z = v1.z - f;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator*(const DOCKVector & v1, const DOCKVector & v2)
{
    DOCKVector      vec;

    vec.x = v1.x * v2.x;
    vec.y = v1.y * v2.y;
    vec.z = v1.z * v2.z;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator*(const float &f, const DOCKVector & v2)
{
    DOCKVector      vec;

    vec.x = f * v2.x;
    vec.y = f * v2.y;
    vec.z = f * v2.z;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator*(const DOCKVector & v1, const float &f)
{
    DOCKVector      vec;

    vec.x = v1.x * f;
    vec.y = v1.y * f;
    vec.z = v1.z * f;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator/(const DOCKVector & v1, const float &f)
{
    DOCKVector      vec;

    vec.x = v1.x / f;
    vec.y = v1.y / f;
    vec.z = v1.z / f;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
float
dot_prod(const DOCKVector & v1, const DOCKVector & v2)
{
    float           dp;

    dp = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

    return dp;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
cross_prod(const DOCKVector & v1, const DOCKVector & v2)
{
    DOCKVector      vec;

    vec.x = v1.y * v2.z - v1.z * v2.y;
    vec.y = -v1.x * v2.z + v1.z * v2.x;
    vec.z = v1.x * v2.y - v1.y * v2.x;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
float
get_vector_angle(const DOCKVector & v1, const DOCKVector & v2)
{
    float           mag,
                    prod;
    float           result;

    mag = v1.length() * v2.length();
    prod = dot_prod(v1, v2) / mag;

    if (prod < -0.999999)
        prod = -0.9999999f;

    if (prod > 0.9999999)
        prod = 0.9999999f;

    if (prod > 1.0)
        prod = 1.0f;

    //result = (acos(prod) / PI) * 180;
    result = (acos(prod) / 3.141592653589793) * 180;
    return result;
}

// +++++++++++++++++++++++++++++++++++++++++
float
get_torsion_angle(DOCKVector & v1, DOCKVector & v2, DOCKVector & v3,
                  DOCKVector & v4)
{
    float           torsion;
    DOCKVector      b1,
                    b2,
                    b3,
                    c1,
                    c2,
                    c3;

    b1 = v1 - v2;
    b2 = v2 - v3;
    b3 = v3 - v4;

    c1 = cross_prod(b1, b2);
    c2 = cross_prod(b2, b3);
    c3 = cross_prod(c1, c2);

    if (c1.length() * c2.length() < 0.001) {
        torsion = 0.0;
    } else {
        torsion = get_vector_angle(c1, c2);
        if (dot_prod(b2, c3) > 0.0)
            torsion *= -1.0;
    }

    return (torsion);
}

//+++++++++++++++++++++++++++++++++++++++++
int
BREADTH_SEARCH::get_search_radius(DOCKMol & mol, int root_atom, int avoid_atom,vector<int> & rotateatoms)
{
    int             root;
    int             x,
                    max_radius;
    int             nbr_atom;
    int             i;
    vector < int   >atom_nbrs;

    atoms.clear();
    nbrs.clear();
    nbrs_next.clear();

    atoms.resize(mol.num_atoms, -1);
    nbrs.push_back(root_atom);
    atoms[root_atom] = 0;
    atoms[avoid_atom] = -2;

    while (nbrs.size() > 0) {
        root = nbrs[nbrs.size() - 1];
        nbrs.pop_back();

        // atom_nbrs = mol.get_atom_neighbors(root);
        atom_nbrs = mol.neighbor_list[root];
	
        for (i = 0; i < atom_nbrs.size(); i++) {
            nbr_atom = atom_nbrs[i];

            if (atoms[nbr_atom] == -1) {
                atoms[nbr_atom] = atoms[root] + 1;
                nbrs.push_back(nbr_atom);
                rotateatoms.push_back(nbr_atom);
            }
        }

        if ((nbrs.size() == 0) && (nbrs_next.size() > 0)) {
            nbrs = nbrs_next;
            nbrs_next.clear();
        }

    }

    max_radius = 0;
    for (x = 0; x < atoms.size(); x++) {
        if (atoms[x] > max_radius)
            max_radius = atoms[x];
    }

    return max_radius;
}

/*************************************/
void
DOCKMol::prepare_molecule()
{
    int             i,
                    j,
                    idx;

    // pre-cache neighbor list
    for (i = 0; i < num_atoms; i++)
        neighbor_list[i] = get_atom_neighbors(i);

        /**
	// pre-cache ie_neighbor_list
	for(i=0;i<num_atoms;i++) {
		for(j=0;j<num_atoms;j++) {
			if((get_bond(i, j)==-1)&&(!atoms_are_one_three(i, j))&&(!atoms_are_one_four(i, j))&&(i!=j))
				ie_neighbor_list[i*num_atoms + j] = false;
			else
				ie_neighbor_list[i*num_atoms + j] = true;
		}
	}
	**/

    // pre-cache atom child list
    for (i = 0; i < num_atoms; i++) {
        for (j = 0; j < num_atoms; j++) {
            if (i != j) {
                idx = get_bond(i, j);
                if (idx != -1) {
                    if (i < j)
                        idx = 2 * idx;
                    else
                        idx = 2 * idx + 1;

                    atom_child_list[idx] = get_atom_children(i, j);
                }
            }
        }
    }


}

