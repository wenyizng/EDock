#ifndef ENERGY_H
#define ENERGY_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#include "mol.h"
#include "base_score.h"
#include "flex_dock.h"
/***********************************/
class           ATOM_TYPE_NODE {

  public:
    char            type[6];
    int             include;
    int             next_total;
    int             vector_atom;
    int             multiplicity;
    float           weight;
    std::vector < ATOM_TYPE_NODE > next;
};

/***********************************/
class           ATOM_TYPE {

  public:
    char            name[100];
    char            atom_model;
    float           radius;
    float           well_depth;
    int             heavy_flag;
    int             valence;
    int             bump_id;
    float           gbradius;
    float           gbscale;
    std::vector < ATOM_TYPE_NODE > definitions;
    

                    //ATOM_TYPE();
                    //virtual ~ ATOM_TYPE();
};

/***********************************/

/***********************************/
class           ATOM_TYPER {

  public:
    std::vector < ATOM_TYPE > types;
    std::vector < int   >atom_types;

    void            get_vdw_labels(std::string fname);
    void print_vector(std::vector<ATOM_TYPE>& types);
    int             assign_vdw_labels(DOCKMol &, int);
    // float vdw_radius_from_type(int type){return types[type].radius;};
    // float well_depth_from_type(int type){return types[type].well_depth;};

};
/***********************************/
class           BOND_TYPE {

  public:
    char            name[100];
    int             drive_id;
    int             minimize;
    int             torsion_total;
    std::vector < float >torsions;
    ATOM_TYPE_NODE  definition[2];
    
};

/***********************************/
class           BOND_TYPER {

  public:
    std::vector < BOND_TYPE > types;
    std::vector < int   >flex_ids;
    int             total_torsions;

    void get_flex_search(std::string fname);
    void print_vector();
    void get_flex_labels(std::string fname);
    void apply_flex_labels(DOCKMol &);

    //bool            is_rotor(int);
};
/***********************************/
class           AMBER_TYPER {

  public:

    // parameter file locations
    std::string     vdw_defn_file;
    std::string     flex_defn_file;
    std::string     flex_drive_tbl;
    std::string     chem_defn_file;
    char            atom_model;
    int             verbose;

    // amber atom and bond typing classes
    ATOM_TYPER      atom_typer;
    BOND_TYPER      bond_typer;
    //CHEM_TYPER      chem_typer;

    //void            initialize(bool read_vdw, bool read_gb_parm, bool use_chem);
     void initialize(string, string, string);
    //void            input_parameters(Parameter_Reader & parm, bool read_vdw, bool use_chem);
    void            assign_hbond_labels(DOCKMol &);
    void            prepare_molecule(DOCKMol &);
    float           getMW (DOCKMol & mol);
    
};
/****************************************************************/
// +++++++++++++++++++++++++++++++++++++++++
class           XYZCRD {
  public:
    float           x,
                    y,
                    z;

    XYZCRD(float _x = 0.0, float _y = 0.0, float _z = 0.0)
        : x(_x), y(_y), z(_z) {
    };
    ~XYZCRD() {
        x = y = z = 0.0;
    };

    void assign_vals(float a, float b, float c) {
        x = a;
        y = b;
        z = c;
    };
    float distance_squared( const XYZCRD & p ) const {
        return (p.x - x)*(p.x - x) + (p.y - y)*(p.y - y) + (p.z - z)*(p.z - z);
    };
    float distance( const XYZCRD & p ) const {
        return sqrt( distance_squared( p ) );
    };

    XYZCRD & operator=(const XYZCRD & xyz) {
        x = xyz.x;
        y = xyz.y;
        z = xyz.z;
        return (*this);
    };

};

typedef         std::vector < XYZCRD > XYZVec;


class           Base_Grid {

  public:

    //Base_Grid();
    //virtual ~Base_Grid() ;

    int             span[3];          //number of bins for each dimension
    int             size;             //total number of items in array
    float           spacing,          //grid space size (in Ang)
                    origin[3];        //starting corner of grid
    float           x_min, x_max,     //min and max coordinates
                    y_min, y_max,     //in xyz dimensions
                    z_min, z_max;  
    XYZCRD          corners[8];       //xyz coordinates of corners of box
    int             neighbors[8];     //array index of all cubes surrounding
                                      //atom of interest
    float           cube_coords[3];   //xyz coordinates of cube of interest
                                      //translated to origin of 0,0,0
    int             nearest_neighbor; //closest corner of cube to atom
                                      //of interest
    unsigned char  *bump;             //array to hold bumps
                                      //needed as base variable because
                                      //header information about all grids
                                      //is contained in bump grid and thus
                                      //must be read in by all grids

    void            read_header(std::string filename);
    void            calc_corner_coords();
    bool            is_inside_grid_box(float x, float y, float z);
    void            find_grid_neighbors(float x, float y, float z);
    float           interpolate(float *grid);

  private:

    int             find_grid_index(int, int, int);

};
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class Energy_Grid is the energy grid from the grid program .nrg file.
// The energy grid is a community resource.
// This class uses the Singleton pattern:
// use Energy_Grid :: get_instance(std::string filename)
// to access the energy grid.
//
class       Energy_Grid:public Base_Grid {

public:

    float  *avdw, *bvdw, *es;
    int     atom_model, att_exp, rep_exp;
    // 2011-01-07 -- trent e balius
    // move from private to public inorder to
    // reset got_the_grid so that more than 
    // one grid can be read in.
    static bool got_the_grid;

    // functions
    ~Energy_Grid() ;
    Energy_Grid();
    void   clear_grid();
    void get_instance(std::string filename);

private:

    void   read_energy_grid(std::string file_prefix);

    // functions

};


class    Energy_Score:public Base_Score {
 
    public:

    Energy_Grid *   energy_grid;

    std::string     grid_file_name;
    float           vdw_scale;
    float           es_scale;

    float           vdw_component;
    float           es_component;

                    //Energy_Score();
                    //virtual ~ Energy_Score();
    //void            input_parameters(Parameter_Reader & parm,
                                     //bool & primary_score,
                                     //bool & secondary_score);
    void            initialize(AMBER_TYPER &);

    bool            compute_score(DOCKMol & mol,float cutoff);
    //std::string     output_score_summary(DOCKMol & mol);
};

/*****************************************************************/
class    Descriptor_Energy_Score:public Base_Score {

  public:

    std::string     receptor_filename;
    DOCKMol         receptor;
    float           vdw_scale;
    float           es_scale;
    float           hb_scale;
    float           fp_vdw_scale;
    float           fp_es_scale;
    float           fp_hb_scale;
    float           internal_scale;

/*
    float           rot_bonds_scale;
    float           heavyatoms_scale;
    float           hb_don_scale;
    float           hb_acc_scale;
    float           molecular_wt_scale;
    float           formal_charge_scale;
    float           XlogP_scale;
*/

    float           rep_exp, att_exp;
    float           diel_screen;
    bool            use_ddd;

    float           vdw_component;
    float           es_component;
    //hbond add
    int             hbond;
    //footprintadd 
    double          vdw_foot_dist;
    double          es_foot_dist;
    double          hbond_foot_dist;
    std::string     desc_foot_compare_type;

    // footprint reference
    bool            bool_footref_mol2;
    bool            bool_footref_txt;
    bool            desc_foot_comp_all_residue;
    // std::string     desc_range;

    bool            desc_normalize_foot;
    
    //RANGE           ref_ranges;
    //RANGE           pose_ranges;
    std::string     desc_fp_info;
    //std::string     desc_range_vdw;
    //std::string     desc_range_es;
    //std::string     desc_range_hb;

    // this are only used if desc_foot_comp_all_residue is false, i.e. not all residues are used.
    bool            desc_foot_specify_a_range; 
    bool            desc_foot_specify_a_threshold;
    bool            desc_use_remainder; 

    double          desc_vdw_threshold;
    double          desc_es_threshold;
    double          desc_hb_threshold;

    int             desc_vdw_num_resid;
    int             desc_es_num_resid;
    int             desc_hb_num_resid;


    DOCKMol         footprint_reference;
    bool            constant_footprint_ref;
    std::string     constant_footprint_ref_file;
    // stop footprint

                    Descriptor_Energy_Score();
                    virtual ~ Descriptor_Energy_Score();
    //                Descriptor_Energy_Score();

    // ADD IN footprint reader the reads in a file in the 
    // same formate as the output and stors the footprints
    // will read in h-bond, vdw, and coul.  
    // footprints are calculated outside dockrun.


   
    // one can pass a reference molecule to calculate the
    // reference footprints with in dockrun. 
    // footprint reference
    //void            submit_footprint_reference(AMBER_TYPER &);
    //bool            compute_footprint(DOCKMol &);
    //double          calc_fp_similarity(DOCKMol &,char,char);
    //double          Euclidean_distance(std::vector <double> *,std::vector <double> *); 
    //double          Pearson_correlation(std::vector <double> *,std::vector <double> *); 
    // stop footprint

    //void            input_parameters(Parameter_Reader & parm,
    //bool & primary_score,
    //bool & secondary_score);

    //void            initialize(AMBER_TYPER &);
    //void            close();
    //bool            compute_score(DOCKMol &);
    //std::string     output_score_summary(DOCKMol &);
    //bool            on_range(int,std::string);
//    void            range_satifying_threshold(std::string &, double, std::string);
    //RANGE           range_satifying_threshold_all(std::vector <FOOTPRINT_ELEMENT>);
    //std::string     union_of_ranges(std::string,std::string);
};


bool energy (AMBER_TYPER & amber, DOCKMol & mol, DOCKMol & receptor, Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,
            string grid_file, float cutoff,vector<vector<int>> biolip_matrix,vector<float>sphx,vector<float>sphy,vector<float> sphz,float bindsite_weight,
            vector<int> protein_atom_index,vector<vector <float> > ave, vector<vector <float> > std);
float tor_energy_total (AMBER_TYPER &,vector<vector<int>> ,vector<TORSION> &);
float tor_energy(AMBER_TYPER & amber,vector<vector<int>> biolip_matrix,TORSION torsion);
bool energy_internal(AMBER_TYPER & amber, DOCKMol & mol);
bool hbond_cal(DOCKMol &, DOCKMol &, int, int ,double &, double &);
float hbond_energy (DOCKMol &mol, DOCKMol &receptor);
float contact_energy (DOCKMol & mol, vector<float> sphx, vector<float> sphy, vector<float> sphz);
float distance_energy (DOCKMol & mol, DOCKMol & receptor, vector<int> protein_atom_index,vector<vector <float> > ave, vector<vector <float> > std);
#endif
