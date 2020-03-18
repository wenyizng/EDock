//
#ifndef MATCH_H
#define MATCH_H 

#include <string>
#include <vector>
#include "mol.h"
#include "energy.h"


class           Sphere {

  public:
    XYZCRD          crds;
    float           radius;
    int             surface_point_i;
    int             surface_point_j;
    int             critical_cluster;
    std::string     color;

    Sphere() {
        clear();
    };
    ~Sphere() {
        clear();
    };

    void            clear();
    float           distance(Sphere &);
    Sphere & operator=(const Sphere & s) {
        crds = s.crds;
        radius = s.radius;
        surface_point_i = s.surface_point_i;
        surface_point_j = s.surface_point_j;
        critical_cluster = s.critical_cluster;
        color = s.color;
        return (*this);
    };

};

typedef         std::vector< Sphere > SphereVec;




class Active_Site_Spheres {
public:
    typedef SphereVec :: const_iterator const_iterator ;
    static SphereVec & get_instance(const char *);
    //static void set_sphere_file_name( Parameter_Reader & parm );
private:
    Active_Site_Spheres();
    static std::string sphere_file_name;
    static SphereVec the_instance;
};

int read_spheres( std::string sphere_file_name, SphereVec & spheres );

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           CLIQUE {
  public:
    INTVec nodes;
    double           residual;

    bool            operator<(CLIQUE clique) const {
        return (residual < clique.residual);
    };
};

typedef         std::vector < CLIQUE > CLIQUEVec;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           CLIQUE_EDGE {

  public:
    int             spheres[2];
    int             centers[2];
    double           residual;

    bool            operator<(CLIQUE_EDGE clnode) const {
        return (residual < clnode.residual);
    };
    bool            operator>(CLIQUE_EDGE clnode) const {
        return (residual > clnode.residual);
    };
    void            operator=(CLIQUE_EDGE cn) {
        spheres[0] = cn.spheres[0];
        spheres[1] = cn.spheres[1];
        centers[0] = cn.centers[0];
        centers[1] = cn.centers[1];
        residual = cn.residual;
    };
};


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           CRITICAL_CLUSTER {

  public:
    int             index;
    INTVec          spheres;
};

typedef         std::vector < CRITICAL_CLUSTER > CLUSTERVec;


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           Orient {

  public:
    // General data members
    bool            orient_ligand;      // flag to use orienting routines
    DOCKMol         original;           // Original molecule
//    DOCKMol         cached_orient;      // Cached orientation
    bool            last_orient_flag;   // flag to control iterations
                                        // (false=more orients, true=last
                                        // orient)
    int             num_orients;        // number of orientations produced
    bool            critical_points;    // flag for use of critical point
                                        // matching
    bool            use_ligand_spheres;
    std::string     lig_sphere_filename;

    // Sphere & Center data
    SphereVec       centers;
    SphereVec       spheres;
    int             num_spheres;//binding sites number 
    int             num_centers;// non-h atoms number
    int             num_nodes;
    double         *sph_dist_mat; // binding sites distance matrix
    double         *lig_dist_mat;//atoms distance matrix
    double         *residual_mat;//
    std::vector < CLIQUE_EDGE > edges;

    // matching data structures

/*
    bool           *candset;
    bool           *notset;
    int            *state;
    double         *level_residuals;
*/

    int             level;
    XYZVec          clique_spheres;
    XYZVec          clique_centers;
    int             clique_size;
    DOCKVector      spheres_com;
    DOCKVector      centers_com;
    double          rotation_matrix[3][3];

    CLIQUEVec       cliques;
    int             current_clique;

    // File I/O members
    std::string     sphere_file_name;   // Sphere infile name

    // Matching parameters
    double          tolerance; // 
    double          orig_tolerance;
    double          dist_min;
    int             min_nodes;
    int             max_nodes;
    int             max_orients;
    bool            automated_matching;
    int             am_iteration_num;

   

    // critical sphere parameters
    CLUSTERVec      receptor_critical_clusters;

    bool            verbose;

    // Member functions
    void            prepare_receptor(const char *); // read in receptor params
    void            get_spheres(const char *);      // read in receptor spheres
    void            match(DOCKMol &, const char *);    // generate information needed
    void            calculate_sphere_distance_matrix(); // calculate matrix of
                                                        // sphere-sphere
                                                        // diststances
    void            get_centers(DOCKMol &);     // Generate ligand centers from
    void            calculate_ligand_distance_matrix(); // generate matrix of
                                                        // center-center
                                                       // distances
    void            id_all_cliques();   // test fxn to find all cliques
    bool            check_clique_critical_points(CLIQUE &);
    bool            new_next_orientation(DOCKMol &);
    void            new_extract_coords_from_clique(CLIQUE &);
    void            calculate_translations();
    void            translate_clique_to_origin();
    void            calculate_rotation();
    //bool            check_clique_chemical_match(CLIQUE &);
    /*
    Orient();
    virtual ~ Orient();
    void            input_parameters(Parameter_Reader & parm);  // Read in
                                                                // parameters
    void            initialize(int argc, char **argv);  // initialize
    
    
    
    
                                                // for matching
    void            get_centers(DOCKMol &);     // Generate ligand centers from 
                                                // heavy atom positions
    void            calculate_ligand_distance_matrix(); // generate matrix of
                                                        // center-center
                                                        // distances
    void            clean_up(); // frees allocated memory
    bool            next_orientation(DOCKMol &);        // calls for a new
                                                        // orientation, manages 
                                                        // control flow for
                                                        // last orient
    bool            generate_orientation(DOCKMol &);    // generates a new
                                                        // orientation
    void            extract_coords_from_clique();       // extracts the spheres 
                                                        // and centers from a
                                                        // clique
    void            calculate_translations();
    void            translate_clique_to_origin();
    void            calculate_rotation();
    bool            more_orientations();
    
    bool            new_next_orientation(DOCKMol &);
    void            new_extract_coords_from_clique(CLIQUE &);
    
    void            read_chem_match_tbl();
    */
};



//void match(DOCKMol &, const char *);
#endif  // MATCH_H
