#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string.h>


#if !defined(MOL_H)
    #define MOL_H 1



using namespace std;

typedef std::vector < float >TOR_LIST;
typedef std::vector < int >INTVec;
typedef std::vector < float >FLOATVec;

// +++++++++++++++++++++++++++++++++++++++++
class           DOCKVector {

  public:
    float           x;
    float           y;
    float           z;

    DOCKVector      operator=(const DOCKVector &);
                    DOCKVector & operator+=(const DOCKVector &);
                    DOCKVector & operator+=(const float &);
                    DOCKVector & operator-=(const DOCKVector &);
                    DOCKVector & operator-=(const float &);
                    DOCKVector & operator*=(const float &);
                    DOCKVector & operator/=(const float &);
                    DOCKVector & normalize_vector();


    float           length() const;
    float           squared_length() const;
    float           squared_dist(const DOCKVector &) const;

    friend int      operator==(const DOCKVector &, const DOCKVector &);
    friend int      operator!=(const DOCKVector &, const DOCKVector &);

    friend DOCKVector operator+(const DOCKVector &, const DOCKVector &);
    friend DOCKVector operator+(const float &, const DOCKVector &);
    friend DOCKVector operator+(const DOCKVector &, const float &);
    friend DOCKVector operator-(const DOCKVector &, const DOCKVector &);
    friend DOCKVector operator-(const DOCKVector &);
    friend DOCKVector operator-(const DOCKVector &, const float &);
    friend DOCKVector operator*(const DOCKVector &, const DOCKVector &);
    friend DOCKVector operator*(const float &, const DOCKVector &);
    friend DOCKVector operator*(const DOCKVector &, const float &);
    friend DOCKVector operator/(const DOCKVector &, const float &);
    //bool operator<(DOCKMol &, DOCKMol &);

    friend float    dot_prod(const DOCKVector &, const DOCKVector &);
    friend DOCKVector cross_prod(const DOCKVector &, const DOCKVector &);
    friend float    get_vector_angle(const DOCKVector &, const DOCKVector &);
    friend float    get_torsion_angle(DOCKVector &, DOCKVector &, DOCKVector &,
                                      DOCKVector &);
};

class DOCKMol {

  public:

    // Molecule General Info
    std::string     title;
    std::string     mol_info_line;
    std::string     comment1;
    std::string     comment2;
    std::string     comment3;
    std::string     score_text_data;
    std::string     energy;
    std::string     simplex_text;
    std::string     mol_data;   // misc data

    // Molecule Size and Property Info
    int             num_atoms;
    int             num_bonds;
    int             num_residues; // this is the number of residues in this molecule 
                                  // this is not the number in the receptor unless 
                                  // the molecule is the receptor
    float           total_dsol; // kxr
    

    bool            arrays_allocated;

    // atom information
    std::string    *atom_data;  // misc data

    float          *x;
    float          *y;
    float          *z;


    bool *flag_acceptor; // this is a h-bond acceptor
    bool *flag_donator; // this is true for polar hydrogen
    int  *acc_heavy_atomid; // atom id of acceptor that is connected to polar h
                           // equal to zero if not polar hydrogen

    float          *charges;
    std::string    *atom_types;
    std::string    *atom_names;


    //footprint info jwu
    std::string    *atom_number;
    std::string    *atom_residue_numbers;
    std::string    *subst_names; // residue name


    std::string    *atom_color; // kxr205
    float          *atom_psol;  // kxr205
    float          *atom_apsol; // kxr205

    // DTM - 11-12-08 - array to store rigid segment info
    int                        *atom_segment_ids;

    // bond information
    int            *bonds_origin_atom;
    int            *bonds_target_atom;
    std::string    *bond_types;

    // ring info
    bool           *atom_ring_flags;
    bool           *bond_ring_flags;

    // activation info
    int             num_active_atoms;
    int             num_active_bonds;
    bool           *atom_active_flags;
    bool           *bond_active_flags;

    // AMBER atom type info
    bool            amber_at_assigned;
    int            *amber_at_id;
    float          *amber_at_radius;
    float          *amber_at_well_depth;
    int            *amber_at_heavy_flag;
    int            *amber_at_valence;
    int            *amber_at_bump_id;

    // Molecule Descriptors for database filter
    unsigned int    rot_bonds;
    unsigned int    heavy_atoms;
    float           mol_wt;
    float           formal_charge;
    //unsigned int    hb_donors;
    //unsigned int    hb_acceptors;
    //float           xlogp;


    // Pharmacophore Type info
    bool            chem_types_assigned;
    std::string    *chem_types;

    // AMBER bond type info
    bool            amber_bt_assigned;
    int            *amber_bt_id;
    int            *amber_bt_minimize;
    int            *amber_bt_torsion_total;
    FLOATVec       *amber_bt_torsions;

    // GB radius & scale factor
    float          *gb_hawkins_radius;
    float          *gb_hawkins_scale;

    // Multi-Grid Ensemble
    int             grid_num;

    // Amber Score Identity
    std::string     amber_score_ligand_id;

    // Arrays obtained from within xlogp.cpp
    // Removed by sudipto and trent until xlogp works again
  /*  int             number_of_H_Donors;
    int             number_of_H_Acceptors;
    int             *H_bond_donor;
    int             *H_bond_acceptor;  */

    // General scoring information
    float           current_score;
    float           internal_energy;
    float           intral_energy; //add by wenyizng
    float           energy_tor_total;
    float           contact_energy;   
    float           hb_energy; 
	float           dis_energy;

    std::string     current_data;  //components of score from scoring function
    std::string     primary_data;
    std::string     hbond_text_data;


    // for footprints score
//    std::vector < double > footprint_vdw;
//    std::vector < double > footprint_es;
//    std::vector < int >    footprint_hb;
//    std::vector <FOOTPRINT_ELEMENT> footprints;
    // Neighbor lists
    INTVec         *neighbor_list;
    bool *ie_neighbor_list;
    INTVec         *atom_child_list;    // list of children atoms for each bond 
                                        // in mol (2x#bonds- for
                                        // directionality)

    // Internal Energy members
     bool use_internal_energy;
     int att_exp;
     int rep_exp;
     float dielectric;
    //add by wenyizhang
    //int centerflag;
    //int clusterflag;
    //bool removeflag;
    
    void allocate_arrays(int, int,int=1000); //jwu
    
    void id_ring_atoms_bonds();
    void find_rings(std::vector < int >, std::vector < int >,
                           std::vector < bool > &, std::vector < bool > &, int);
    std::vector < int > get_atom_neighbors(int);
    int get_bond(int, int);
    std::vector < int > get_atom_children(int, int);
    bool bond_is_rotor(int index) {return amber_bt_id[index] != -1;};
    bool atoms_are_one_three(int, int);
    bool atoms_are_one_four(int, int);

    void clear_molecule();
    void translate_mol(const DOCKVector & vec);
    void rotate_mol(double mat[9]);
    void rotate_mol(double mat[3][3]);
    void           prepare_molecule();
    
    bool operator <(DOCKMol anothermol) const{
        return (current_score<anothermol.current_score);
    };
    
};

class Rotatemol{
    public: 
        DOCKMol startmol; //will use the rotet vector and tanslate vector to change the mol
        DOCKMol endmol; // the mol after the change
        vector<float> axisA;
        vector<float> axisB;// rotate vector
        vector<float> tranlate;// translate vector
        float angle;
	int centerflag;
	int clusterflag;
	bool removeflag;
        bool operator<(Rotatemol mol) const {
            return (endmol.current_score<mol.endmol.current_score);
        };
    
};

/********************************************************************/
class           BREADTH_SEARCH {

  public:

    INTVec atoms;
    INTVec          nbrs;
    INTVec          nbrs_next;

    int             get_search_radius(DOCKMol &, int, int, vector<int> &);

};

typedef std::vector < Rotatemol > RotatemolVec;


void copy_molecule(DOCKMol &, DOCKMol &);
#endif 
