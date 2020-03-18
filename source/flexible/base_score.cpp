#include <iostream>
#include "energy.h"
#include "base_score.h"
#include "mol.h"

using namespace std;


// +++++++++++++++++++++++++++++++++++++++++
// static member initializers

const string Base_Score::DELIMITER    = "########## ";
// 20 is magic
const int    Base_Score::FLOAT_WIDTH  = 20;
// length of "Score:" plus length of longest score name.
const int    Base_Score::STRING_WIDTH = 7 + 7;


// +++++++++++++++++++++++++++++++++++++++++
Base_Score::Base_Score()
{
    use_internal_energy = false;
    vdwA = NULL;
    vdwB = NULL;
    ie_vdwA = NULL;
    ie_vdwB = NULL;
    nb_int.clear();
    rep_radius_scale = 1.0;
    // consider scaling the repulsive energy during ligand growth
    // add the question tree later?
}

// +++++++++++++++++++++++++++++++++++++++++
Base_Score::~Base_Score()
{
    delete[]vdwA;
    delete[]vdwB;
    delete[]ie_vdwA;
    delete[]ie_vdwB;
    nb_int.clear();
}

// +++++++++++++++++++++++++++++++++++++++++
// Calculate Internal Energy of the ligand kxr 010506
// DTM - 11-12-08
// I changed this to only compute interactions between atoms of different segments.
// Also to skip any atoms that are inactive, so
// this can be used during flexible growth.

// Trent & Sudipto 2009-02-12 -- This function is all so used for 
// single point and rigid minimazation. (tortions move dering fin min.)
float
Base_Score::compute_ligand_internal_energy(DOCKMol & mol)

{
 
    float           distancesq;
    unsigned int    i;
    unsigned int    a1,a2;
    float           ligand_internal_energy = 0.0;
    float           int_vdw_rep = 0.0;
    //cout<<"use_internal_energy is "<<use_internal_energy<<endl; 
    if (use_internal_energy) {
        for (i = 0; i < nb_int.size(); i++) {
            //cout << "nb_int.size() is :"<< nb_int.size()<<endl;
            // nb_int is a neighbor list of all non-bonded pair of ligand atoms
               a1 = nb_int[i].first;
	       //cout << "a1 is:" << a1 << ","; 
               a2 = nb_int[i].second;
               //cout << "a2 is:" << a2 << endl;
               //cout << mol.atom_active_flags[a1]<<"," << mol.atom_active_flags[a2]<<endl;
               if ( mol.atom_active_flags[a1] && mol.atom_active_flags[a2] )
                {
                        distancesq = pow((mol.x[a1] - mol.x[a2]), 2)
                                   + pow((mol.y[a1] - mol.y[a2]), 2) 
                                   + pow((mol.z[a1] - mol.z[a2]), 2);
                        if(distancesq == 0) cout<< a1 << "," << a2<<endl;
                        //cout<< "distancesq is:"<<distancesq<<endl;
 
 			//ie_rep_exp = 12
                        //cout<< ie_vdwA[a1] << "," << ie_vdwA[a2] << endl;
                        int_vdw_rep += (ie_vdwA[a1]*ie_vdwA[a2]) / pow(distancesq, float(ie_rep_exp/2.0));
                        //cout << "int_vdw_rep is:" << int_vdw_rep << endl;
 
                }
        }
    }
    ligand_internal_energy = int_vdw_rep;
    mol.internal_energy = ligand_internal_energy;
    //cout<< "mol.internal_energy is:"<<mol.internal_energy<<endl;
    //cout<<"mol.internal_energy="<<ligand_internal_energy<<endl;
    return ligand_internal_energy;
}


// Bad design warning: sudipto
// This function stores data for each specific ligand inside class base_score
// this is bad because this is the parent class of all scoring functions, 
// and should be a stateless function wrt each ligand being scored 
// However, since the simplex minimizer needs to call this function as well,
// but there is no way out but to store the precalculated VDW A&B as well as 
// non-bonded atoms list in here. Let me know you if find a better place
// to store the precalculated data.
// +++++++++++++++++++++++++++++++++++++++++
void            
Base_Score::initialize_internal_energy(DOCKMol & mol)
{
   
    //ie_att_exp, ie_rep_exp, ie_diel must be set before calling this
    //clear the ie_vdw arrays
    method = 0;
    use_internal_energy = 1;
    ie_att_exp = 6;
    ie_rep_exp = 12;
    ie_diel = 2.0;
    delete[]ie_vdwA;
    ie_vdwA = NULL;
    delete[]ie_vdwB;
    ie_vdwB = NULL;

    ie_vdwA = new float[mol.num_atoms];
    ie_vdwB = new float[mol.num_atoms];

    // calculate vdwA and vdwB terms
    for (int i = 0; i < mol.num_atoms; i++) {
        ie_vdwA[i] = sqrt(mol.amber_at_well_depth[i] *
                 (ie_att_exp / (ie_rep_exp - ie_att_exp)) *
                 pow((2 * mol.amber_at_radius[i]), ie_rep_exp));

        ie_vdwB[i] = sqrt(mol.amber_at_well_depth[i] *
                 (ie_rep_exp / (ie_rep_exp - ie_att_exp)) *
                 pow((2 * mol.amber_at_radius[i]), ie_att_exp));
    }

   // neighbor list for internal energy 
   nb_int.clear();  //clear list of non-bonded interactions
   int a1,a2;      // atom number counters
   INTPair temp;   // to push_back into array nb_int 
  
   // this double loop makes sure we never check the same atom against itself 
   for (a1 = 0; a1 < mol.num_atoms - 1; a1++) {
       for (a2 = a1 + 1; a2 < mol.num_atoms; a2++) {
       
          // nonbondpair is used in the if statement below.
          // if neither method is set then all are false

          bool nonbondpair = false;  
          //cout<<"method is:"<<method<<endl;
          switch (method) {
       
          case 0: //Rigid Docking, Single point calcs, Minimization without
                  // orientation -- perform all atom internal.
              // for rigid, segments are not assinged, all segment_id's are -1
              nonbondpair = ((mol.get_bond(a1, a2) == -1) // a1 & a2 are non-bonded
                   && (!mol.atoms_are_one_three(a1, a2)) // not 1-3
                   && (!mol.atoms_are_one_four(a1, a2)) // not 1-4
                   );
              break;
       
          case 1: //Flex Anchor and Grow -- perform inter-segment internal
              nonbondpair = ( (mol.atom_segment_ids[a1] != mol.atom_segment_ids[a2])
                              //not within same segment
                   && (mol.get_bond(a1, a2) == -1) // a1 & a2 are non-bonded
                   && (!mol.atoms_are_one_three(a1, a2)) // not 1-3
                   && (!mol.atoms_are_one_four(a1, a2)) // not 1-4
                   );

              break;
       
          // case 2: //HDB
              // break;
          }

          if ( nonbondpair) // different criteria for different methods 
           {
                temp.first = a1;
                temp.second = a2;
                nb_int.push_back(temp);
		//cout << nb_int[nb_int.size()-1].first<< "  "<< nb_int[nb_int.size()-1].second<<endl;
           }
       }   // for loop a1
   }  // for loop a2

}


// +++++++++++++++++++++++++++++++++++++++++
void
Base_Score::init_vdw_energy(AMBER_TYPER & typer, float att_exp, float rep_exp)
{

    //att_exp and rep_exp are read from the energy grid file
    delete[]vdwA;
    vdwA = NULL;

    delete[]vdwB;
    vdwB = NULL;

    vdwA = new float[typer.atom_typer.types.size()];
    vdwB = new float[typer.atom_typer.types.size()];
    //cout<<"typer.atom_typer.types.size() is:"<<typer.atom_typer.types.size()<<endl;
    for (int i = 0; i < typer.atom_typer.types.size(); i++) {
        vdwA[i] = sqrt(typer.atom_typer.types[i].well_depth *
                 (att_exp / (rep_exp - att_exp)) *
                 pow((2 * rep_radius_scale * typer.atom_typer.types[i].radius), rep_exp));
	//cout<<"vdwA["<<i<<"]:"<<vdwA[i]<<endl;

        vdwB[i] = sqrt(typer.atom_typer.types[i].well_depth *
                 (rep_exp / (rep_exp - att_exp)) *
                 pow((2 * typer.atom_typer.types[i].radius), att_exp));
       //cout<<"vdwB["<<i<<"]:"<<vdwB[i]<<endl;
    }
}
