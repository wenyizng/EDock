#ifndef ENERGY_CPP
#define ENERGY_CPP 1

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include <assert.h>

#include "mol.h"
#include "energy.h"
#include "base_score.h"

#define MIN_FLOAT -1e+15
#define PI 3.1415926535897932384626433f
#define NINT(x) (int) ((x) > 0 ? ((x) + 0.5) : ((x) - 0.5))
#define INTFLOOR(x) (int) (floor(x + 0.00001))
#define MAX(x, y) ((x) > (y) ? (x) : (y))


char           *
white_line(char *line)
{
    for (unsigned int i = 0; i < strlen(line); i++)
        if (isspace(line[i]))
            line[i] = ' '; 

    return line;
}


int
assign_node(ATOM_TYPE_NODE & node, int include)
{
    char            temp[6];
    char           *branch;

    node.next_total = 0;
    strcpy(temp, strtok(NULL, " "));

    if (isdigit(temp[0])) {
        node.multiplicity = atoi(temp);

        if (node.multiplicity < 0) {
            cout << "Cannot specify negative multiplicity in definition.\n";
            return 0;
        }

        strcpy(node.type, strtok(NULL, " "));

    } else {
        node.multiplicity = 0;
        strcpy(node.type, temp);
    }

    node.include = include;

    while (branch = strtok(NULL, " ")) {

        if ((!strncmp(branch, "(", 1)) || (!strncmp(branch, "[", 1))) {

            if (node.next_total >= 6) {
                cout <<
                    "Cannot exceed 6 substituents for each atom in definition.\n";
                return 0;
            }

            ATOM_TYPE_NODE  tmp_node;
            node.next.push_back(tmp_node);

            if ((!strncmp(branch, "(", 1))
                && (!assign_node(node.next[node.next_total], 1)))
                return 0;

            if ((!strncmp(branch, "[", 1))
                && (!assign_node(node.next[node.next_total], 0)))
                return 0;

            node.next_total++;
        }

        else if ((!strncmp(branch, ")", 1)) || (!strncmp(branch, "]", 1)))
            return 1;

        else
            return 0;
    }

    return 1;
}

int
count_atom_neighbors(DOCKMol & mol, int atom_num)
{
    int             count;

    count = mol.get_atom_neighbors(atom_num).size();

    return count;
}

/*
 * =================================================================== 
 */
int
check_type(const char *candidate, const char *reference)
{
    if ((strstr(candidate, reference)) || (reference[0] == '*'))
        return 1;
    else
        return 0;
}

/*
 * =================================================================== 
 */

/*
 * =================================================================== 
 */
int
check_bonded_atoms(DOCKMol & mol, int current_atom, int previous_atom,
                   ATOM_TYPE_NODE & node)
{
    int             i,
                    j,
                    next_atom;
    int             match,
                    match_count;
    vector < int   >nbrs;

    match_count = 0;

    // loop over neighbors of current atom
    nbrs = mol.get_atom_neighbors(current_atom);
    for (i = 0; i < nbrs.size(); i++) {

        next_atom = nbrs[i];
        match = 0;

        if ((next_atom != previous_atom)
            && (match =
                check_type(mol.atom_types[next_atom].c_str(), node.type))) {

            for (j = 0; j < node.next_total; j++) {
                if (!check_bonded_atoms
                    (mol, next_atom, current_atom, node.next[j]))
                    match = 0;
            }

        }

        if (match)
            match_count++;
    }

    if (node.multiplicity) {
        if (node.multiplicity == match_count)
            match = 1;
        else
            match = 0;
    } else {
        if (match_count)
            match = 1;
        else
            match = 0;
    }

    if (match == node.include)
        return 1;
    else
        return 0;

}

/*
 * =================================================================== 
 */


/*
 * =================================================================== 
 */
int
check_atom(DOCKMol & mol, int current_atom, ATOM_TYPE_NODE & node)
{
    int             match = 0;
    int             i;

    if (match = check_type(mol.atom_types[current_atom].c_str(), node.type)) {

        for (i = 0; i < node.next_total; i++) {

            if (!check_bonded_atoms
                (mol, current_atom, current_atom, node.next[i]))
                match = 0;

        }

    }

    return match;
}


/*
 * =================================================================== 
 */
bool
count_bond_neighbors(DOCKMol & mol, int bond_num)
{
    int             c1,
                    c2;

    c1 = mol.get_atom_neighbors(mol.bonds_origin_atom[bond_num]).size();
    c2 = mol.get_atom_neighbors(mol.bonds_target_atom[bond_num]).size();

    if ((c1 <= 1) || (c2 <= 1))
        return false;
    else
        return true;
}

/*
 * =================================================================== 
 */


/***********************************/
void
ATOM_TYPER::get_vdw_labels(string fname)
{
    FILE           *ifp;
    char            line[100],
                    model[100];
    ATOM_TYPE       tmp_type;
    ATOM_TYPE_NODE  tmp_node;
    int             i;

    ifp = fopen(fname.c_str(), "r");

    if (ifp == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    while (fgets(line, 100, ifp) != NULL) {

        if (!strncmp(line, "name", 4)) {        // read in name field
            types.push_back(tmp_type);

            if (sscanf(line, "%*s %s", types[types.size() - 1].name) < 1) {
                cout << "Incomplete vdw member declaration.\n";
                exit(0);
            }

        } else if (!strncmp(line, "atom_model", 10)) {  // read in atom model
                                                        // field
            if (sscanf(line, "%*s %s", model) != 1) {
                cout << "Incomplete atom_model specification.\n";
                exit(0);
            }

            types[types.size() - 1].atom_model = tolower(model[0]);

            if ((types[types.size() - 1].atom_model != 'a')
                && (types[types.size() - 1].atom_model != 'u')
                && (types[types.size() - 1].atom_model != 'e')) {
                cout <<
                    "Atom_model specification restricted to ALL, UNITED, or EITHER.\n";
                exit(0);
            }
        } else if (!strncmp(line, "heavy_flag", 8)) {   // read in heavy flag
                                                        // field
            if (sscanf(line, "%*s %d", &types[types.size() - 1].heavy_flag) !=
                1) {
                cout << "Incomplete heavy_flag specification.\n";
                exit(0);
            }
        } else if (!strncmp(line, "radius", 6)) {       // read in radius field
            if (sscanf(line, "%*s %f", &types[types.size() - 1].radius) != 1) {
                cout << "Incomplete radius specification.\n";
                exit(0);
            }
        } else if (!strncmp(line, "well_depth", 10)) {  // read in well depth
                                                        // field
            if (sscanf(line, "%*s %f", &types[types.size() - 1].well_depth) !=
                1) {
                cout << "Incomplete well_depth specification.\n";
                exit(0);
            }
        } else if (!strncmp(line, "valence", 7)) {      // read in valence
                                                        // field
            if (sscanf(line, "%*s %d", &types[types.size() - 1].valence) != 1) {
                cout << "Incomplete valence specification.\n";
                exit(0);
            }
        } else if (!strncmp(line, "gbradii", 7)) {
            if (sscanf(line, "%*s %f", &types[types.size() - 1].gbradius) != 1) {
                cout << "Incomplete GBRadius specification.\n";
                exit(0);
            }
        } else if (!strncmp(line, "gbscale", 7)) {
            if (sscanf(line, "%*s %f", &types[types.size() - 1].gbscale) != 1) {
                cout << "Incomplete GBScale specification.\n";
                exit(0);
            }
        } else if (!strncmp(line, "definition", 10)) {  // read in definition
                                                        // fields
            strtok(white_line(line), " ");
            types[types.size() - 1].definitions.push_back(tmp_node);

            if (!assign_node
                (types[types.size() - 1].
                 definitions[types[types.size() - 1].definitions.size() - 1],
                 true)) {
                cout << "Error assigning vdw member definitions.\n";
                exit(0);
            }
        }

    }                           // end While

    // check to see if any vdw values have been read in
    if (types.size() < 1) {
        cout << "ERROR:  VDW parameter file empty." << endl;
        exit(0);
    }
    // check to see if all values for vdw parameters have been read in
    for (i = 0; i < types.size() - 1; i++) {
        if (types[i].atom_model == '\0') {
            cout << "ERROR:  No atom_model assigned for " << types[i].
                name << " in VDW parameter file." << endl;
            exit(0);
        }
        if (types[i].radius == MIN_FLOAT) {
            cout << "ERROR:  No radius assigned for " << types[i].
                name << " in VDW parameter file." << endl;
            exit(0);
        }
        if (types[i].well_depth == MIN_FLOAT) {
            cout << "ERROR:  No well_depth assigned for " << types[i].
                name << " in VDW parameter file." << endl;
            exit(0);
        }
        if (types[i].heavy_flag == INT_MIN) {
            cout << "ERROR:  No heavy_flag assigned for " << types[i].
                name << " in VDW parameter file." << endl;
            exit(0);
        }
        if (types[i].valence == INT_MIN) {
            cout << "ERROR:  No valence assigned for " << types[i].
                name << " in VDW parameter file." << endl;
            exit(0);
        }
       
            if (types[i].gbradius == MIN_FLOAT) {
                cout << "ERROR:  No gbradii assigned for " << types[i].
                    name << " in VDW parameter file." << endl;
                exit(0);
            }

            if (types[i].gbscale == MIN_FLOAT) {
                cout << "ERROR:  No gbscale assigned for " << types[i].
                    name << " in VDW parameter file." << endl;
                exit(0);
            }
        
        if (types[i].definitions.size() < 1) {
            cout << "ERROR:  No definitions assigned for " << types[i].
                name << " in VDW parameter file." << endl;
            exit(0);
        }

    }

    for (i = 0; i < types.size() - 2; i++) {    // calculate bump_id values
        if (types[i].heavy_flag)
            types[i].bump_id = NINT(10.0 * types[i].radius);
        else
            types[i].bump_id = 0;
    }

    fclose(ifp);
}

void ATOM_TYPER :: print_vector(vector<ATOM_TYPE>& types) {
    
    for (int i=0;i<types.size();i++){
        cout<<"types["<<i<<"]"<<".name is "<< types[i].name<<endl;
        cout<<"types["<<i<<"]"<<".atom_model is "<< types[i].atom_model<<endl;
        cout<<"types["<<i<<"]"<<".radius is "<< types[i].radius<<endl;
        cout<<"types["<<i<<"]"<<".well_depth is "<< types[i].well_depth<<endl;
        cout<<"types["<<i<<"]"<<".heavy_flag is "<< types[i].heavy_flag<<endl;
        cout<<"types["<<i<<"]"<<".valence is "<< types[i].heavy_flag<<endl;
        cout<<"types["<<i<<"]"<<".bump_id is "<< types[i].bump_id<<endl;
        cout<<"types["<<i<<"]"<<".gbradius is "<< types[i].gbradius<<endl;
        cout<<"types["<<i<<"]"<<".gbscale is "<< types[i].gbscale<<endl;
        for(int j=0;j<types[i].definitions.size();j++){
            cout<<"types["<<i<<"]"<<".definitions.type is "<< types[i].definitions[j].type<<endl;
            cout<<"types["<<i<<"]"<<".definitions.include is "<< types[i].definitions[j].include<<endl;
            cout<<"types["<<i<<"]"<<".definitions.next_total is "<< types[i].definitions[j].next_total<<endl;
            cout<<"types["<<i<<"]"<<".definitions.vector_atom is "<< types[i].definitions    
                                                                     [j].vector_atom<<endl;
            cout<<"types["<<i<<"]"<<".definitions.multiplicity is "<< types[i].definitions
                                                                     [j].multiplicity<<endl;
            cout<<"types["<<i<<"]"<<".definitions.weight is "<< types[i].definitions[j].weight<<endl;
            
        }
             
    }
}

/***********************************/
void
BOND_TYPER::get_flex_labels(string fname)
{
    FILE           *ifp;
    char            line[100];
    BOND_TYPE       tmp_flex;
    int             i;
    int             definition_count;

    ifp = fopen(fname.c_str(), "r");

    if (ifp == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    while (fgets(line, 100, ifp)) {

        if (!strncmp(line, "name", 4)) {
            definition_count = 0;

            types.push_back(tmp_flex);

            if (types.size() > 1) {
            }

            if (sscanf(line, "%*s %s", types[types.size() - 1].name) < 1) {
                cout << "Incomplete Flex Definition Failure.\n";
                exit(0);
            }

            for (i = 0; i < strlen(types[types.size() - 1].name); i++)
                types[types.size() - 1].name[i] =
                    (char) tolower(types[types.size() - 1].name[i]);

            types[types.size() - 1].drive_id = -1;
            types[types.size() - 1].minimize = -1;

        }                       // End "name" if

        else if (!strncmp(line, "drive_id", 6))
            sscanf(line, "%*s %d", &types[types.size() - 1].drive_id);

        else if (!strncmp(line, "minimize", 8))
            sscanf(line, "%*s %d", &types[types.size() - 1].minimize);

        else if (!strncmp(line, "definition", 10)) {

            strtok(white_line(line), " ");
            assign_node(types[types.size() - 1].definition[definition_count],
                        1);
            definition_count++;

        }                       // End "definition" if
    }
    fclose(ifp);
}


/***********************************/
void
BOND_TYPER::get_flex_search(string fname)
{
    FILE           *ifp;
    char            line[100];
    char           *token;
    int             i;

    int             drive_id;
    int             torsion_total;
    FLOATVec        torsions;

    ifp = fopen(fname.c_str(), "r");

    if (ifp == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    for (i = 0; i < types.size(); i++)
        types[i].torsions.clear();

    while (fgets(line, 100, ifp)) {

        token = strtok(white_line(line), " ");

        if (!strcmp(token, "drive_id")) {

            torsions.clear();

            if (token = strtok(NULL, " "))
                drive_id = atoi(token);

            else {
                cout <<
                    "ERROR get_flex_search: Search_id value not specified in ";
                cout << fname << endl;
                exit(0);
            }

            if (!fgets(line, 100, ifp)
                || !(token = strtok(white_line(line), " "))
                || strcmp(token, "positions")) {
                cout <<
                    "ERROR get_flex_search: Positions field doesn't follow Id in ";
                cout << fname << endl;
                exit(0);
            }

            if (token = strtok(NULL, " "))
                torsion_total = atoi(token);
            else {
                cout <<
                    "ERROR get_flex_search: Postions value not specified in ";
                cout << fname << endl;
                exit(0);
            }

            if (!fgets(line, 100, ifp)
                || !(token = strtok(white_line(line), " "))
                || strcmp(token, "torsions")) {
                cout <<
                    "ERROR get_flex_search: Torsions doesn't follow Positions in ";
                cout << fname << endl;
                exit(0);
            }

            for (i = 0; i < torsion_total; i++) {
                if (token = strtok(NULL, " ")) {
                    torsions.push_back(atof(token));
                } else {
                    cout <<
                        "ERROR get_flex_search: Insufficient number of torsions in ";
                    cout << fname << endl;
                    exit(0);
                }

            }

            for (i = 0; i < types.size(); i++) {
                if (types[i].drive_id == drive_id) {
                    types[i].torsion_total = torsion_total;
                    types[i].torsions = torsions;
                }
            }

        }                       // End if drive_id

    }                           // End while


    for (i = 0; i < types.size(); i++) {
        if (types[i].torsion_total < 1) {
            cout << "ERROR get_flex_search: Missing torsion parameters in ";
            cout << fname << endl;
            exit(0);
        }
    }

    fclose(ifp);


}

void BOND_TYPER:: print_vector() {
    for (int i=0;i<types.size();i++){
        cout<<i<<":"<<endl;
        cout<<"name is:"<<types[i].name<<endl;
        cout<<"drive_id is:"<<types[i].drive_id<<endl;
        cout<<"minimize is:"<<types[i].minimize<<endl;
        cout<<"torsion_total is:"<<types[i].torsion_total<<endl;
        for (int j=0;j<types[i].torsions.size();j++){
	    cout<<"torsions["<<j<<"]"<<types[i].torsions[j]<<endl;
        }
    }
}

/***********************************/
void
BOND_TYPER::apply_flex_labels(DOCKMol & mol)
{
    int             i,
                    j;

    flex_ids.clear();
    total_torsions = 0;

    for (i = 0; i < mol.num_bonds; i++) {

        flex_ids.push_back(0);
        flex_ids[i] = -1;

        if (mol.bond_ring_flags[i])
            continue;

        if (!count_bond_neighbors(mol, i))
            continue;

        for (j = 0; j < types.size(); j++) {    // Loop over bond type (1 to
                                                // size)

            if (check_atom
                (mol, mol.bonds_origin_atom[i], types[j].definition[0])
                && check_atom(mol, mol.bonds_target_atom[i],
                              types[j].definition[1])) {
                flex_ids[i] = j;
            } else
                if (check_atom
                    (mol, mol.bonds_origin_atom[i], types[j].definition[1])
                    && check_atom(mol, mol.bonds_target_atom[i],
                                  types[j].definition[0])) {
                flex_ids[i] = j;
            }

        }                       // End loop over bond types

        if (flex_ids[i] != -1)
            total_torsions++;
    }
  //print_vector(flex_ids);

}


/***********************************/
// H-bond add
void
AMBER_TYPER::assign_hbond_labels( DOCKMol & mol )
{
 // assign hbond labels
 for (int i=1;i<mol.num_atoms;i++){
      stringstream ss;
      stringstream ss1;
      char atom_type_char[10];
      ss << mol.atom_types[i]; 
      ss >> atom_type_char;
      // see if the atom is a acceptor
      if (    atom_type_char[0] == 'O' 
          ||  atom_type_char[0] == 'N' 
          ||  atom_type_char[0] == 'S'
          ||  atom_type_char[0] == 'F'
          || (atom_type_char[0] == 'C' && atom_type_char[1] == 'l')
          || (atom_type_char[0] == 'C' && atom_type_char[1] == 'L'))
      {
         mol.flag_acceptor[i] = true;
      }
      // see if atom is a donator
      // atoms must be an h and must be conected to a acceptor through a bond
      if(   atom_type_char[0] == 'H')
      {
            for (int j=0;j<mol.num_bonds;j++)
            {
                 if( i == mol.bonds_origin_atom[j])// atom i (which is an h) the start a bond
                 {
                    ss1 << mol.atom_types[mol.bonds_target_atom[j]]; 
                    ss1 >> atom_type_char;
                    if (    atom_type_char[0] == 'O'
                        ||  atom_type_char[0] == 'N'
                        ||  atom_type_char[0] == 'S'
                        ||  atom_type_char[0] == 'F'
                        || (atom_type_char[0] == 'C' && atom_type_char[1] == 'l')
                        || (atom_type_char[0] == 'C' && atom_type_char[1] == 'L'))
                    {
                           mol.flag_donator[i] = true;
                           mol.acc_heavy_atomid[i] = mol.bonds_target_atom[j];
                    }
                 }

                 if( i == mol.bonds_target_atom[j])// atom i (which is an h) the terminates a bond
                 {
                    ss1 << mol.atom_types[mol.bonds_origin_atom[j]]; 
                    ss1 >> atom_type_char;
                    if (    atom_type_char[0] == 'O'
                        ||  atom_type_char[0] == 'N'
                        ||  atom_type_char[0] == 'S'
                        ||  atom_type_char[0] == 'F'
                        || (atom_type_char[0] == 'C' && atom_type_char[1] == 'l')
                        || (atom_type_char[0] == 'C' && atom_type_char[1] == 'L'))
                    {
                           mol.flag_donator[i] = true;
                           mol.acc_heavy_atomid[i] = mol.bonds_origin_atom[j];
                    }
                 }

            }
      }
 }

}
/***********************************/
int
ATOM_TYPER::assign_vdw_labels(DOCKMol & mol, int atom_model)
{
    int             i,
                    j,
                    k;
    int             vdw_assigned = false;
    vector < int   >nbrs;

    atom_types.clear();

    // loop over all atoms in the mol
    for (i = 0; i < mol.num_atoms; i++) {

        // by default, assign dummy type
        atom_types.push_back(0);

        // loop over types read in from vdw defn file
        for (j = 0; j < types.size(); j++) {

            if ((atom_model == 'a') && (types[j].atom_model == 'u'))
                continue;

            if ((atom_model == 'u') && (types[j].atom_model == 'a'))
                continue;

            for (k = 0; k < types[j].definitions.size(); k++) {

                if (check_atom(mol, i, types[j].definitions[k])) {
                    atom_types[i] = j;
                    vdw_assigned = true;
                }               // end k loop

            }

        }                       // end j loop

        // check for vdw label assignment
        if (!vdw_assigned) {
            cout << "WARNING assign_vdw_labels: No vdw parameters for ";
            cout << i << " " << " " << mol.atom_types[i] << endl;
            // return false;
        }
        // check that no atom's valence is violated
        if (count_atom_neighbors(mol, i) > types[atom_types[i]].valence) {

            // loop over neighbor atoms
            nbrs = mol.get_atom_neighbors(i);
            k = 0;

            for (j = 0; j < nbrs.size(); j++)
                if ((strcmp(mol.atom_types[nbrs[j]].c_str(), "LP"))
                    && (strcmp(mol.atom_types[nbrs[j]].c_str(), "Du")))
                    k++;

            if (k > types[atom_types[i]].valence) {
                cout << "WARNING assign_vdw_labels: Atom valence violated for ";
                cout << mol.title << " ";
                cout << "atom number: " << i << endl;

                // return false;
            }

        }
        // transfer partial charges for united models
        if ((atom_model == 'u')
            && (fabs(types[atom_types[i]].well_depth) < 0.00001)) {
            if (count_atom_neighbors(mol, i) == 1) {
                nbrs = mol.get_atom_neighbors(i);
                mol.charges[nbrs[0]] = mol.charges[nbrs[0]] + mol.charges[i];
                mol.charges[i] = 0.0;
            } else {
                cout <<
                    "WARNING assign_vdw_labels: Unable to transfer partial charge away from ";
                cout << mol.title << " ";
                cout << "atom number: " << i << endl;

                mol.charges[i] = 0.0;
            }
        }
    }                           // end i loop

    return true;

}
/***********************************/
float AMBER_TYPER::getMW(DOCKMol & mol)
{
     //atomic weights from General Chemistry 3rd Edition Darrell D. Ebbing
     float mw = 0.0;
     char sybyl_type[10];

     for (int i = 0; i < mol.num_atoms; i++) {

        strcpy(sybyl_type, mol.atom_types[i].c_str());
        string element = strtok(sybyl_type,".");

        if (element == "O") mw += 15.9994;
        else if (element == "N") mw += 14.00674;
        else if (element == "C") mw += 12.011;
        else if (element == "F") mw += 18.9984032;
        else if (element == "Cl") mw += 35.4527;
        else if (element == "Br") mw += 79.904;
        else if (element == "I") mw += 126.90447;
        else if (element == "H" ) mw += 1.00794;
        else if (element == "B")  mw += 10.811;
        else if (element == "S" ) mw += 32.066;
        else if (element == "P")  mw += 30.973762;
        else if (element == "Li") mw += 6.941;
        else if (element == "Na") mw += 22.98968;
        else if (element == "Mg") mw += 24.3050;
        else if (element == "Al") mw += 26.981539;
        else if (element == "Si") mw += 28.0855;
        else if (element == "K") mw += 39.0983;
        else if (element == "Ca") mw += 40.078;
        else if (element == "Cr") mw += 51.9961;
        else if (element == "Mn") mw += 54.93805;
        else if (element == "Fe") mw += 55.847;
        else if (element == "Co") mw += 58.93320;
        else if (element == "Cu") mw += 63.546;
        else if (element == "Zn") mw += 65.39;
        else if (element == "Se") mw += 78.96;
        else if (element == "Mo") mw += 95.94;
        else if (element == "Sn") mw += 118.710;
        else if (element == "LP")  mw += 0.0;
        else cout << "Element " << element << " not found in MW code\n";
    }
    return mw;
}

/***********************************/
void
AMBER_TYPER::prepare_molecule(DOCKMol & mol )
{
    
    int i;

    // assign hbond acc and donor labels
    assign_hbond_labels(mol);
    
    // use the amber atom typers to assign atom and bond types
    
        atom_typer.assign_vdw_labels(mol, atom_model);
        //cout<< "atom_model is:" << atom_model << endl; 
        bond_typer.apply_flex_labels(mol);// find the rotate bond
	
   
    // copy atom types to molecule
        for (i = 0; i < mol.num_atoms; i++) {
            mol.amber_at_id[i] = atom_typer.atom_types[i];
            mol.amber_at_radius[i] =
                atom_typer.types[atom_typer.atom_types[i]].radius;
		//cout<<atom_typer.types[atom_typer.atom_types[i]].radius<<endl;
            mol.amber_at_well_depth[i] =
                atom_typer.types[atom_typer.atom_types[i]].well_depth;
		//cout<< atom_typer.types[atom_typer.atom_types[i]].well_depth<<endl;
            mol.amber_at_heavy_flag[i] =
                atom_typer.types[atom_typer.atom_types[i]].heavy_flag;
		//cout<<mol.amber_at_heavy_flag[i]<<endl;
            mol.amber_at_valence[i] =
                atom_typer.types[atom_typer.atom_types[i]].valence;
		//cout<<atom_typer.types[atom_typer.atom_types[i]].valence<<endl;
            mol.amber_at_bump_id[i] =
                atom_typer.types[atom_typer.atom_types[i]].bump_id;
		//cout<<atom_typer.types[atom_typer.atom_types[i]].bump_id<<endl;
            mol.gb_hawkins_radius[i] =
                atom_typer.types[atom_typer.atom_types[i]].gbradius;
		//cout<<atom_typer.types[atom_typer.atom_types[i]].gbradius<<endl;
            mol.gb_hawkins_scale[i] =
                atom_typer.types[atom_typer.atom_types[i]].gbscale;
		//cout<<atom_typer.types[atom_typer.atom_types[i]].gbscale<<endl;
        }
        mol.amber_at_assigned = true;
        // copy bond types to molecule
	//cout<<"mol.num_bonds is:" <<mol.num_bonds<<endl;
        for (i = 0; i < mol.num_bonds; i++) {
            
            mol.amber_bt_id[i] = bond_typer.flex_ids[i];
	    //cout<<i<<":"<<bond_typer.flex_ids[i]<<":"<<endl; 
            if (bond_typer.flex_ids[i] != -1) {
		//cout<<"bond_num:"<<i<<endl;
                //cout<<"drive_id:"<<bond_typer.types[bond_typer.flex_ids[i]].drive_id<<endl;
                mol.amber_bt_minimize[i] =
                    bond_typer.types[bond_typer.flex_ids[i]].minimize;
		    //cout<<"minimize:"<<bond_typer.types[bond_typer.flex_ids[i]].minimize<<endl;
                mol.amber_bt_torsion_total[i] =
                    bond_typer.types[bond_typer.flex_ids[i]].torsion_total;
		    //cout<<"types:"<<bond_typer.types[bond_typer.flex_ids[i]].torsion_total<<endl;
                mol.amber_bt_torsions[i] =
                    bond_typer.types[bond_typer.flex_ids[i]].torsions;
		    //cout<<"torsions:"<<endl;
		    //for(int j=0;j<mol.amber_bt_torsions[i].size();j++)
			//cout<<mol.amber_bt_torsions[i][j]<<endl;
		    //print_vector(mol.amber_bt_torsions[i]);

            } else {
                mol.amber_bt_minimize[i] = 0;
                mol.amber_bt_torsion_total[i] = 0;
                mol.amber_bt_torsions[i].clear();
            }

        }
        mol.amber_bt_assigned = true;

        mol.rot_bonds=0;  //defined in dockmol.cpp as part of dbfilter code
        for (i = 0; i < mol.num_bonds; i++)
            if (mol.bond_is_rotor(i)) mol.rot_bonds++;
	//cout<<"mol.rot_bonds is:"<<mol.rot_bonds<<endl;
        mol.formal_charge = 0.0;  //defined in dockmol.cpp as part of dbfilter code
        for (i = 0; i < mol.num_atoms; i++)
            mol.formal_charge += mol.charges[i];

        mol.heavy_atoms = 0;
        // Count the number of heavy atoms
        for (i = 0; i < mol.num_atoms; i++)
            if (mol.amber_at_heavy_flag[i]) mol.heavy_atoms++;

        mol.mol_wt = getMW(mol);  //defined in dockmol.cpp as part of dbfilter code
	//cout << "mol.mol_wt is :"<<endl;
	//cout<<mol.mol_wt<<endl;
        //cout<<"mol.rot_bonds is : "<<mol.rot_bonds<<endl;
   
}


void
AMBER_TYPER::initialize(string filevdw, string fileflex, string fileflex_drive_tbl)
{
     
	//cout <<vdw_defn_file.c_str()<<","<<read_gb_parm<<endl;
        //vdw_AMBER_parm99.defn,0

        atom_typer.get_vdw_labels(filevdw.c_str());
	//atom_typer.print_vector(atom_typer.types);
	
        bond_typer.get_flex_labels(fileflex.c_str());
	
	//cout <<flex_drive_tbl<<endl;
	//flex_drive.tbl
        bond_typer.get_flex_search(fileflex_drive_tbl.c_str());
        //bond_typer.print_vector();
   
    

}


/************************************************/
/*+++++++++++++++++++++ Energy_Grid +++++++++++++*/
// +++++++++++++++++++++++++++++++++++++++++++++++
// class Energy_Grid

// static member initializers
bool Energy_Grid :: got_the_grid = false;

// +++++++++++++++++++++++++++++++++++++++++
// Default constructor is private and a no op.
Energy_Grid :: Energy_Grid()
{
  avdw = NULL;
  bvdw = NULL;
  es   = NULL;
}

// +++++++++++++++++++++++++++++++++++++++++
// Default distructor is public.
Energy_Grid::~Energy_Grid()
{
    this->clear_grid();
}


// +++++++++++++++++++++++++++++++++++++++++
void Energy_Grid::clear_grid()
{
    delete[] avdw;
    delete[] bvdw;
    delete[] es;
}
// +++++++++++++++++++++++++++++++++++++++++
//reads header information from .bmp file
void
Base_Grid::read_header(string filename)
{
    string fname = filename + ".bmp";
    FILE* grid_in;
    grid_in = fopen( fname.c_str(), "rb");
    if (grid_in == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    cout << " Reading the grid box quantifiers from " << fname << endl;
    fread(&size, sizeof(int), 1, grid_in);
    fread(&spacing, sizeof(float), 1, grid_in);
    fread(origin, sizeof(float), 3, grid_in);
    fread(span, sizeof(int), 3, grid_in);
    //cout<<"size is:"<<size<<endl;
    //cout<<"spacing is:" <<spacing<<endl;
    //for(int i=0;i<3;i++){
     //   cout<<"origin["<<i<<"]:"<<origin[i]<<endl;
     //   cout<<"span["<<i<<"]:"<<span[i]<<endl;
   // }
    
    cout << " Done reading the grid box quantifiers." << endl;
    fclose(grid_in);
}

// +++++++++++++++++++++++++++++++++++++++++
//computes the minimum and maximum xyz coordinates for the dimensions
//of the grid, and collects them as the corners of the grid box
void
Base_Grid::calc_corner_coords()
{
    float dx = (span[0] - 1) * spacing;
    float dy = (span[1] - 1) * spacing;
    float dz = (span[2] - 1) * spacing;

    x_min = origin[0];
    x_max = origin[0] + dx;
    y_min = origin[1];
    y_max = origin[1] + dy;
    z_min = origin[2];
    z_max = origin[2] + dz;

    corners[0].assign_vals(x_min, y_min, z_min);
    corners[1].assign_vals(x_max, y_min, z_min);
    corners[2].assign_vals(x_max, y_min, z_max);
    corners[3].assign_vals(x_min, y_min, z_max);
    corners[4].assign_vals(x_min, y_max, z_min);
    corners[5].assign_vals(x_max, y_max, z_min);
    corners[6].assign_vals(x_max, y_max, z_max);
    corners[7].assign_vals(x_min, y_max, z_max);
}


// +++++++++++++++++++++++++++++++++++++++++
void
Energy_Grid::read_energy_grid(string filename)
{
    //cout<<"I am in read_energy_grid"<<endl;
    string fname = filename + ".nrg";
    //string fname = filename + ".bmp";
    FILE* grid_in;
    grid_in = fopen(fname.c_str(), "rb");
    if (grid_in == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    cout << " Reading the energy grid from " << fname << endl;

    fread(&size, sizeof(int), 1, grid_in);
    //cout<<"grid size is:"<<size<<endl;
    fread(&atom_model, sizeof(int), 1, grid_in);
    fread(&att_exp, sizeof(int), 1, grid_in);
    fread(&rep_exp, sizeof(int), 1, grid_in);
    
    //cout<<"size is:"<<size<<endl;
    avdw = new float[size];
    bvdw = new float[size];
    es = new float[size];

    fread( bvdw, sizeof(float), size, grid_in);
    fread( avdw, sizeof(float), size, grid_in);
    fread( es, sizeof(float), size, grid_in);
    /*
    for (int i=0;i<size;i++){
        cout<<"bvdw["<<i<<"]:"<<bvdw[i]<<endl;
        cout<<"avdw["<<i<<"]:"<<avdw[i]<<endl;
        cout<<"es["<<i<<"]:"<<es[i]<<endl;
    }
    */
    cout << " Done reading the energy grid." << endl;
    fclose(grid_in);
    
    read_header(filename);

    calc_corner_coords();
}




void
Energy_Grid :: get_instance(string filename)
{
    //cout<<"got_the_grid is:"<<got_the_grid<<endl;
    if ( ! got_the_grid ) {
        assert( 0 != filename.size() );
        read_energy_grid( filename );
        // if no errors then
        got_the_grid = true;
       
    }
}


// +++++++++++++++++++++++++++++++++++++++++
void
Energy_Score::initialize(AMBER_TYPER & typer)
{

        //cout << "Initializing Grid Score Routines..." << endl;
        energy_grid = new Energy_Grid(); 
        energy_grid->get_instance(grid_file_name);
	
	cout << "grid_file_name is : " << grid_file_name << endl; 
	//cout << "grid size is:"<<sizeof(energy_grid)<<endl;
        init_vdw_energy(typer, energy_grid->att_exp, energy_grid->rep_exp);//grid_bsp_cluster
     
}

// +++++++++++++++++++++++++++++++++++++++++
//determines whether atom is inside the boundaries of the grid box
//if not, return false
bool
Base_Grid::is_inside_grid_box(float x, float y, float z)
{
    //cout<<"x_min is:"<<x_min<<","<<"x_max is:"<<x_max<<endl;
    //cout<<"y_min is:"<<y_min<<","<<"y_max is:"<<y_max<<endl;
    //cout<<"z_min is:"<<z_min<<","<<"z_max is:"<<z_max<<endl;
    //cout<< "spacing is:"<<spacing<<endl;
    //cout<< "x is:"<<x<<","<<"y is:"<<y<<","<<"z is:"<<z<<endl;
    bool reflag = (x > x_min+spacing && x < x_max-spacing &&
              y > y_min+spacing && y < y_max-spacing &&
              z > z_min+spacing && z < z_max-spacing) ;
    //cout<< "reflag is:"<<reflag <<endl;
    return   (x > x_min+spacing && x < x_max-spacing &&
              y > y_min+spacing && y < y_max-spacing &&
              z > z_min+spacing && z < z_max-spacing) ;
    // Since we are using bicubic interpolation in the grid energy compuation,
    // any (x,y,z) point needs to be between x_min+spacing and x_max-spacing
    // in order to return valid scores for all 8 points being interpolated
    // For comparing floats, we need to account for round off error
}

// +++++++++++++++++++++++++++++++++++++++++
//computes index of array where value is stored
int 
Base_Grid::find_grid_index(int x, int y, int z)
{
    return span[0] * span[1] * z + span[0] * y + x;
}

// +++++++++++++++++++++++++++++++++++++++++
//finds the indices of the cubes that are neighboring
//the atom of interest
void
Base_Grid::find_grid_neighbors(float x, float y, float z)
{
    float x_int = x - origin[0];
    float y_int = y - origin[1];
    float z_int = z - origin[2];

    int x_nearest = NINT(x_int / spacing);
    int x_below = INTFLOOR(x_int / spacing);
    int x_above = x_below + 1;
    if (x_nearest >= span[0]) {
        if (x_below >= span[0])
            x_nearest = span[0] - 1;
        else
            x_nearest = x_below;
    }
    if (x_nearest < 0) {
        if (x_above < 0)
            x_nearest = 0;
        else
            x_nearest = x_above;
    }

    int y_nearest = NINT(y_int / spacing);
    int y_below = INTFLOOR(y_int / spacing);
    int y_above = y_below + 1;
    if (y_nearest >= span[1]) {
        if (y_below >= span[1])
            y_nearest = span[1] - 1;
        else
            y_nearest = y_below;
    }
    if (y_nearest < 0) {
        if (y_above < 0)
            y_nearest = 0;
        else
            y_nearest = y_above;
    }

    int z_nearest = NINT(z_int / spacing);
    int z_below = INTFLOOR(z_int / spacing);
    int z_above = z_below + 1;
    if (z_nearest >= span[2]) {
        if (z_below >= span[2])
            z_nearest = span[2] - 1;
        else
            z_nearest = z_below;
    }
    if (z_nearest < 0) {
        if (z_above < 0)
            z_nearest = 0;
        else
            z_nearest = z_above;
    }

    nearest_neighbor = find_grid_index(x_nearest, y_nearest, z_nearest);

    neighbors[0] = find_grid_index(x_above, y_above, z_above);
    neighbors[1] = find_grid_index(x_above, y_above, z_below);
    neighbors[2] = find_grid_index(x_above, y_below, z_above);
    neighbors[3] = find_grid_index(x_below, y_above, z_above);
    neighbors[4] = find_grid_index(x_above, y_below, z_below);
    neighbors[5] = find_grid_index(x_below, y_above, z_below);
    neighbors[6] = find_grid_index(x_below, y_below, z_above);
    neighbors[7] = find_grid_index(x_below, y_below, z_below);

    float corrected_coords[3];
    corrected_coords[0] = x - origin[0];
    corrected_coords[1] = y - origin[1];
    corrected_coords[2] = z - origin[2];

    cube_coords[0] = corrected_coords[0] / spacing - (float) (INTFLOOR(corrected_coords[0] / spacing)); 
    cube_coords[1] = corrected_coords[1] / spacing - (float) (INTFLOOR(corrected_coords[1] / spacing));
    cube_coords[2] = corrected_coords[2] / spacing - (float) (INTFLOOR(corrected_coords[2] / spacing)); 
}

// +++++++++++++++++++++++++++++++++++++++++
//function for interpolating values from a set of 8 grid points
float
Base_Grid::interpolate(float *grid)
{
    float           a1,
                    a2,
                    a3,
                    a4,
                    a5,
                    a6,
                    a7,
                    a8;
    float           value;
    int             out_of_bounds,
                    i;

    out_of_bounds = 0;
    for (i = 0; i < 8; i++)
        if ((neighbors[i] > size) || (neighbors[i] < 0))
            out_of_bounds = 1;
    //cout<<"out_of_bounds is:"<<out_of_bounds<<endl;
    //cout<<"grid[neighbors[7]] is:"<<grid[neighbors[7]]<<endl;
	
    if (out_of_bounds == 0) {
        a8 = grid[neighbors[7]];
        a7 = grid[neighbors[6]] - a8;
        a6 = grid[neighbors[5]] - a8;
        a5 = grid[neighbors[4]] - a8;
        a4 = grid[neighbors[3]] - a8 - a7 - a6;
        a3 = grid[neighbors[2]] - a8 - a7 - a5;
        a2 = grid[neighbors[1]] - a8 - a6 - a5;
        a1 = grid[neighbors[0]] - a8 - a7 - a6 - a5 - a4 - a3 - a2;
        value =
            a1 * cube_coords[0] * cube_coords[1] * cube_coords[2] +
            a2 * cube_coords[0] * cube_coords[1] +
            a3 * cube_coords[0] * cube_coords[2] +
            a4 * cube_coords[1] * cube_coords[2] + a5 * cube_coords[0] +
            a6 * cube_coords[1] + a7 * cube_coords[2] + a8;
	//cout<<"a1 is:"<<a1<<endl;
        return value;
    } else {
        return grid[nearest_neighbor];
    }
}

// +++++++++++++++++++++++++++++++++++++++++
bool
Energy_Score::compute_score(DOCKMol & mol,float cutoff)
{
    int atom;
    
    // check to see if molecule is inside grid box
    for (atom = 0; atom < mol.num_atoms; atom++) {
        //cout<<"mol.atom_active_flags[atom] is:"<<mol.atom_active_flags[atom]<<endl;
        bool testflag;
    // is an active atom and is inside the grid
	
        if (mol.atom_active_flags[atom] && !(energy_grid->is_inside_grid_box(mol.x[atom], mol.y[atom], mol.z[atom]))) {
            //cout<<"testflag is:"<<testflag<<endl;
            mol.current_score = -MIN_FLOAT;  // arbitrarily large score
            mol.current_data = "ERROR:  Conformation could not be scored."
                "\nConformation not completely within grid box.\n";
            //cout<<mol.current_data<<endl;
            return false;
         }
    }

    float es_val = 0.0;
    float vdw_val = 0.0;
    float total = 0.0;
    float vdw_val_tmp=0;
    float es_val_tmp=0;
    vdw_scale = es_scale = 1; // add by zwy
    for (atom = 0; atom < mol.num_atoms; atom++) {
        //cout<<"mol.atom_active_flags[atom] is:"<<mol.atom_active_flags[atom]<<endl;
        if (mol.atom_active_flags[atom]) {

           energy_grid->find_grid_neighbors(mol.x[atom], mol.y[atom], mol.z[atom]);
           //cout<<"vdw_scale: " <<vdw_scale<<endl;
            
	    //cout<<"vdwA[mol.amber_at_id[atom] is:"<<vdwA[mol.amber_at_id[atom]]<<endl;
	    //cout<<"vdwB[mol.amber_at_id[atom] is:"<<vdwB[mol.amber_at_id[atom]]<<endl;
	   float potentialA = energy_grid->interpolate(energy_grid->avdw);
	   float potentialB = energy_grid->interpolate(energy_grid->bvdw);
	   //cout<<"potentialA is:"<<potentialA<<endl;
	   //cout<<"potentialB is:"<<potentialB<<endl;
            vdw_val_tmp =
                ((vdwA[mol.amber_at_id[atom]] * energy_grid->interpolate(energy_grid->avdw)) -
                (vdwB[mol.amber_at_id[atom]] * energy_grid->interpolate(energy_grid->bvdw))) *
                vdw_scale;
	    
            //if(vdw_val_tmp > cutoff)
            //{
            //	vdw_val_tmp = cutoff;
            //}
            vdw_val += vdw_val_tmp;
	    //cout<<"mol.charges is:" << mol.charges[atom]<<endl;

            es_val_tmp = mol.charges[atom] * energy_grid->interpolate(energy_grid->es) * es_scale;
	    //float temp = energy_grid->interpolate(energy_grid->es);
            //cout<<"temp is:"<<temp<<endl;
             
            es_val += es_val_tmp;
            //cout<<"vdw_val_tmp\t"<<vdw_val_tmp<<'\t';
            //cout<<"es_val_tmp\t"<<es_val_tmp<<endl;
        }
    }

    total = vdw_val + es_val;

    vdw_component = vdw_val;
    es_component = es_val;

    //mol.current_score = total;
    mol.intral_energy=total;
    
    //mol.current_data = output_score_summary(mol);
    //cout<<"mol.current_score="<<mol.current_score<<endl;

    return true;
}
float contact_energy (DOCKMol & mol, vector<float> sphx, vector<float> sphy, vector<float> sphz)
{
    float contact_energy = 0.0;
    
    for(int i = 0 ;i<sphx.size();i++)
    {
         float mindis = 1000.0;
         for(int j = 0;j<mol.num_atoms;j++)
         {
             float dij =sqrt((sphx[i]-mol.x[j])*(sphx[i]-mol.x[j]) + (sphy[i]-mol.y[j])*(sphy[i]-mol.y[j]) + (sphz[i]-mol.z[j])*(sphz[i]-mol.z[j]));
             //cout<<dij<<endl;
             if (mindis > dij)
                  mindis = dij; 
	     //cout<<"mindis "<<mindis<<endl;
         }
         contact_energy = contact_energy + mindis;

    }
    return contact_energy;
}
/***********************************************/

bool energy (AMBER_TYPER & amber, DOCKMol & mol,Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file, float cutoff,vector<float> sphx, vector<float> sphy, vector<float> sphz) {
   
    Base_Score score;
    float energy_inter_score= 0.0;
    float internal_energy = 0.0 ;
    
   
    
    score.initialize_internal_energy(mol);
    //cout<<"score.initialize_internal_energy finish!"<<endl;
    score.compute_ligand_internal_energy (mol);
    //cout<<"score.compute_ligand_internal_energy finish!"<<endl;
    // read the receptor and calculate the ligand and receptor energy
    if(c_nrg.compute_score(mol,cutoff)){
	//cout<<"c_nrg.compute_score finish"<<endl;
        //energy_inter_score = mol.current_score;
        energy_inter_score = mol.intral_energy;
        internal_energy = mol.internal_energy;
        //mol.current_score = mol.intral_energy + mol.internal_energy;
	float contact_score = contact_energy(mol, sphx,sphy,sphz);
	//cout<<"contact_score "<<contact_score<<endl;
        mol.current_score=mol.intral_energy + contact_score;
        //cout<< "energy_inter_score is:" <<energy_inter_score<<endl;
        //cout<< "internal_energy is:" <<internal_energy<<endl;
        //return energy_inter_score, internal_energy;
        return true;
    }
    else
        return false;
    
}

#endif
