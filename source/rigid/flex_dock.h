#include <map>
#include <string>
#include <utility>  // pair
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mol.h"
#include "translate_rotate.h"

#if !defined(FLEX_DOCK_H)
    #define FLEX_DOCK_H 1
// +++++++++++++++++++++++++++++++++++++++++
class           TORSION {

  public:

    int             atom1;      // first atom bound to atom 2 (to define angle)
    int             atom2;      // first atom of rotatable bond
    int             atom3;      // second atom of rotatable bond
    int             atom4;      // first atom bound to atom 3 (to define angle)
    int             bond_num;   // bond number in Mol
    vector <int>    Latoms;     //add by wenyizng
    vector <int>    Ratoms;     //add by wenyizng
    float ihedral;              //Points2Dihedral

};
void
id_torsions(DOCKMol & mol, FLOATVec & vertex, vector<TORSION> &torsions);
#endif 

