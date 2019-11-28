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

float contact_energy (DOCKMol & mol, vector<float> sphx, vector<float> sphy, vector<float> sphz);
bool energy (AMBER_TYPER & amber, DOCKMol & mol,Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file, float cutoff,vector<float> sphx, vector<float> sphy, vector<float> sphz);

#endif
