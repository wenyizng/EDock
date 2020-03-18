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
#include <algorithm>
#include <stdio.h>


#include "mol.h"
#include "MathTools.h"
#include "match.h"
#include "output.h"
#include "energy.h"

/*** matfit stuff ***/
#define SMALL  1.0e-20
#define SMALSN 1.0e-10
#define ABS(x)   (((x)<0)   ? (-(x)) : (x))

// +++++++++++++++++++shpere class++++++++++++++++++++++
void
Sphere::clear()
{
    crds.assign_vals(0.0, 0.0, 0.0);
    radius = 0.0;
    surface_point_i = 0;
    surface_point_j = 0;
    critical_cluster = 0;
    color.clear();
}


// +++++++++++++++++++++++++++++++++++++++++
float
Sphere::distance(Sphere & sph)
{
    return crds.distance(sph.crds);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class Active_Site_Spheres

// static member initializers
string Active_Site_Spheres :: sphere_file_name = string();
SphereVec Active_Site_Spheres :: the_instance = vector< Sphere >();

// +++++++++++++++++++++++++++++++++++++++++
SphereVec &
Active_Site_Spheres :: get_instance(const char * recfilename)
{
    sphere_file_name = recfilename;
    if ( 0 == the_instance.size() ) {
        assert( 0 != sphere_file_name.size() );
        int count = read_spheres( sphere_file_name, the_instance );
	cout<<"sphere_file_name is :"<<sphere_file_name<<endl;
	//cout <<"the_instance is: " <<the_instance<<endl;
        assert( the_instance.size() == count );
    }
    return the_instance;
}

int
read_spheres( string sphere_file_name, SphereVec & spheres )
{
    char            junk[50],
                    line[500],
                    label[50];
    Sphere          tmp;
    FILE           *sphere_file;
    int             i;
    string          color_label;
    int             color_int;
    vector < string > site_color_labels;  // sphere file chem type labels
    INTVec          site_color_ints;      // sphere file chem type numbers

    int num_spheres = 0;

    // open the sphere file
    
    sphere_file = fopen(sphere_file_name.c_str(), "r");

    if (sphere_file == NULL) {
        cout << "\n\nCould not open " << sphere_file_name <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    // read in the color table (if any)
    while (fgets(line, 500, sphere_file)) {
        sscanf(line, "%s", junk);
	//cout<<line<<endl;

        // if line is a color table entry
        if (strcmp(junk, "color") == 0) {
            sscanf(line, "%s %s %d", junk, label, &color_int);
            color_label = label;

            site_color_labels.push_back(color_label);
            site_color_ints.push_back(color_int);
        }
        // if color table is finished
        if (strcmp(junk, "cluster") == 0)
            break;
    }

    // read in the spheres
    while (fgets(line, 500, sphere_file)) {

        sscanf(line, "%s", junk);
        if (strcmp(junk, "cluster") == 0)
            break;
        else
            sscanf(line, "%d %f %f %f %f %d %d %d", &tmp.surface_point_i,
                   &tmp.crds.x, &tmp.crds.y, &tmp.crds.z, &tmp.radius,
                   &tmp.surface_point_j, &tmp.critical_cluster, &color_int);
	    //cout<<line<<endl;//radius important !!!!!
        // assign the proper color label
        tmp.color = "null";

        if (color_int > 0) {
            for (i = 0; i < site_color_ints.size(); i++) {
                if (color_int == site_color_ints[i]) {
                    tmp.color = site_color_labels[i];
                    break;
                }
            }
        }

        spheres.push_back(tmp);
        num_spheres++;
    }

    fclose(sphere_file);

    return num_spheres;
}

/************************************************************************/
static void
minimized_fit(double umat[3][3], double rm[3][3])
{
    // Sudipto: could someone add comments to explain what each varible
    // here describes, and more details about the algorithm used?

    double          rot[3][3],
                    turmat[3][3],
                    c[3][3],
                    coup[3],
                    dir[3],
                    step[3],
                    v[3],
                    rtsum,
                    rtsump,
                    rsum,
                    stp,
                    stcoup,
                    ud,
                    tr,
                    ta,
                    cs,
                    sn,
                    ac,
                    delta,
                    deltap,
                    gfac,
                    cle,
                    clep;
    int             i,j,k,l,m,  //loop counters
                    jmax,
                    ncyc,
                    nsteep,
                    nrem;

    /*
     * Rotate repeatedly to reduce couple about initial direction to zero.
     * Clear the rotation matrix
     */
    for (l = 0; l < 3; l++) {
        for (m = 0; m < 3; m++)
            rot[l][m] = 0.0;
        rot[l][l] = 1.0;
    }

    /*
     * Copy vmat[][] (sp) into umat[][] (dp) 
     */
    jmax = 30;
    rtsum = umat[0][0] + umat[1][1] + umat[2][2];
    delta = 0.0;

    for (ncyc = 0; ncyc < jmax; ncyc++) {
        /*
         * Modified CG. For first and every NSTEEP cycles, set previous step as 
         * zero and do an SD step 
         */
        nsteep = 3;
        nrem = ncyc - nsteep * (int) (ncyc / nsteep);

        if (!nrem) {
            for (i = 0; i < 3; i++)
                step[i] = 0.0;
            clep = 1.0;
        }

        /*
         * Couple 
         */
        coup[0] = umat[1][2] - umat[2][1];
        coup[1] = umat[2][0] - umat[0][2];
        coup[2] = umat[0][1] - umat[1][0];
        cle = sqrt(coup[0] * coup[0] + coup[1] * coup[1] + coup[2] * coup[2]);

        /*
         * Gradient vector is now -coup 
         */
        gfac = (cle / clep) * (cle / clep);

        /*
         * Value of rtsum from previous step 
         */
        rtsump = rtsum;
        deltap = delta;
        clep = cle;
        if (cle < SMALL){
            //cout <<ncyc<< ": cle < SMALL:" << cle << "<" << SMALL << endl;
            break;
        }

        /*
         * Step vector conjugate to previous 
         */
        stp = 0.0;
        for (i = 0; i < 3; i++) {
            step[i] = coup[i] + step[i] * gfac;
            stp += (step[i] * step[i]);
        }
        stp = 1.0 / sqrt(stp);

        /*
         * Normalised step 
         */
        for (i = 0; i < 3; i++)
            dir[i] = stp * step[i];

        /*
         * Couple resolved along step direction 
         */
        stcoup = coup[0] * dir[0] + coup[1] * dir[1] + coup[2] * dir[2];

        /*
         * Component of UMAT along direction 
         */
        ud = 0.0;
        for (l = 0; l < 3; l++)
            for (m = 0; m < 3; m++)
                ud += umat[l][m] * dir[l] * dir[m];


        tr = umat[0][0] + umat[1][1] + umat[2][2] - ud;
        ta = sqrt(tr * tr + stcoup * stcoup);
        cs = tr / ta;
        sn = stcoup / ta;

        /*
         * If cs<0 then posiiton is unstable, so don't stop 
         */
        if ((cs > 0.0) && (ABS(sn) < SMALSN)){
            //cout <<ncyc <<": (cs > 0.0) && (ABS(sn) < SMALSN):"<< cs << ">"<< 0.0 << " && " << ABS(sn) << "<" << SMALSN << endl;
            break;
        }

        /*
         * Turn matrix for correcting rotation:
         * 
         * Symmetric part 
         */
        ac = 1.0 - cs;
        for (l = 0; l < 3; l++) {
            v[l] = ac * dir[l];
            for (m = 0; m < 3; m++)
                turmat[l][m] = v[l] * dir[m];
            turmat[l][l] += cs;
            v[l] = dir[l] * sn;
        }

        /*
         * Asymmetric part 
         */
        turmat[0][1] -= v[2];
        turmat[1][2] -= v[0];
        turmat[2][0] -= v[1];
        turmat[1][0] += v[2];
        turmat[2][1] += v[0];
        turmat[0][2] += v[1];

        /*
         * Update total rotation matrix 
         */
        for (l = 0; l < 3; l++) {
            for (m = 0; m < 3; m++) {
                c[l][m] = 0.0;
                for (k = 0; k < 3; k++)
                    c[l][m] += turmat[l][k] * rot[k][m];
            }
        }

        for (l = 0; l < 3; l++)
            for (m = 0; m < 3; m++)
                rot[l][m] = c[l][m];

        /*
         * Update umat tensor 
         */
        for (l = 0; l < 3; l++)
            for (m = 0; m < 3; m++) {
                c[l][m] = 0.0;
                for (k = 0; k < 3; k++)
                    c[l][m] += turmat[l][k] * umat[k][m];
            }

        for (l = 0; l < 3; l++)
            for (m = 0; m < 3; m++)
                umat[l][m] = c[l][m];

        rtsum = umat[0][0] + umat[1][1] + umat[2][2];
        delta = rtsum - rtsump;

        /*
         * If no improvement in this cycle then stop 
         */
        if (ABS(delta) < SMALL){
            //cout <<ncyc<< ": ABS(delta) < SMALL: " << ABS(delta) << "<" << SMALL << endl;
            break;
        }

        /*
         * Next cycle 
         */
    }

    rsum = rtsum;

    /*
     * Copy rotation matrix for output 
     */

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            rm[i][j] = rot[i][j];       // can be transposed
}

/*************************************************************************/
int
compute_rot_matrix(XYZVec & x1, XYZVec & x2, double rm[3][3], int n)
{
    int             i,
                    j;
    double          umat[3][3];


    if (n < 2) {
        return (0);
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            umat[i][j] = 0.0;
    }

    for (j = 0; j < n; j++) {
        umat[0][0] += x1[j].x * x2[j].x;
        umat[1][0] += x1[j].y * x2[j].x;
        umat[2][0] += x1[j].z * x2[j].x;

        umat[0][1] += x1[j].x * x2[j].y;
        umat[1][1] += x1[j].y * x2[j].y;
        umat[2][1] += x1[j].z * x2[j].y;

        umat[0][2] += x1[j].x * x2[j].z;
        umat[1][2] += x1[j].y * x2[j].z;
        umat[2][2] += x1[j].z * x2[j].z;
    }

    minimized_fit(umat, rm);

    return (1);
}

/************************************************/
void
Orient::get_spheres(const char * recfilename)
{
    bool            found_cluster;
    CRITICAL_CLUSTER tmp_cluster;//orient.h

    spheres = Active_Site_Spheres :: get_instance(recfilename); // save binding sites in sphere.cpp
    num_spheres = spheres.size();

    //if (verbose) cout << "Read in " << num_spheres 

    // process critical points
    //cout<<"critical_points is:"<<critical_points<<endl;
    //critical_points is 0
    /*
    if (critical_points) {
        receptor_critical_clusters.clear();

        // loop over spheres
        for (int i = 0; i < spheres.size(); i++) {

            // if sphere belongs to a critical cluster
            if (spheres[i].critical_cluster > 0) {

                found_cluster = false;

                for (int j = 0; j < receptor_critical_clusters.size(); j++) {
                    // if the cluster has been identified previously
                    if (spheres[i].critical_cluster ==
                        receptor_critical_clusters[j].index) {
                        // add the sphere to the existing cluster
                        receptor_critical_clusters[j].spheres.push_back(i);
                        found_cluster = true;
                        break;
                    }
                }

                // else create new cluster entry
                if (!found_cluster) {
                    tmp_cluster.index = spheres[i].critical_cluster;
                    tmp_cluster.spheres.clear();
                    tmp_cluster.spheres.push_back(i);
                    receptor_critical_clusters.push_back(tmp_cluster);
                }
            }
        }
    }
    */
}

void
Orient::prepare_receptor(const char* recfilename)
{
    get_spheres(recfilename);// get the binding site
    calculate_sphere_distance_matrix();
    
}

/************************************************/
void
Orient::calculate_sphere_distance_matrix()
{

    sph_dist_mat = new double[num_spheres * num_spheres];

    for (int i = 0; i < num_spheres; i++)
        for (int j = 0; j < num_spheres; j++) {
            sph_dist_mat[num_spheres * i + j] = spheres[i].distance(spheres[j]); 
	    // caclute the i j distance //utils.h
	    //cout<<sph_dist_mat[num_spheres * i + j]<<endl;
        }
}

/************************************************/
void
Orient::get_centers(DOCKMol & mol)
{
    centers.clear();
    num_centers = 0;
    //cout << use_ligand_spheres << endl; // use_ligand_spheres is 0
    

    Sphere tmp;
    //cout<<"mol.num_atoms is:" <<mol.num_atoms<<endl;
    
    for (int atom = 0; atom < mol.num_atoms; atom++) {
        //cout<<"mol.amber_at_heavy_flag[atom] is:"<<mol.amber_at_heavy_flag[atom]<<endl;
        if (mol.amber_at_heavy_flag[atom]) {
            if (mol.atom_active_flags[atom]) {
                tmp.crds.x = mol.x[atom];
                tmp.crds.y = mol.y[atom];
                tmp.crds.z = mol.z[atom];
                tmp.radius = 0.0;
                tmp.surface_point_i = 0;
                tmp.surface_point_j = 0;
                tmp.critical_cluster = 0;
                centers.push_back(tmp);
                num_centers++;
            }
        }
     }
   
}
    

/************************************************/
void
Orient::calculate_ligand_distance_matrix()
{

    lig_dist_mat = new double[num_centers * num_centers];

    for (int i = 0; i < num_centers; i++)
        for (int j = 0; j < num_centers; j++) {
            lig_dist_mat[num_centers * i + j] = centers[i].distance(centers[j]);
            //cout<<lig_dist_mat[num_centers * i + j]<<endl;
        }
}


/************************************************/
bool
Orient::check_clique_critical_points(CLIQUE & clique)
{
    Sphere          tmp;
    int             index,
                    sphere;
    int             i,
                    j,
                    k;
    INTVec          tmp_clique_spheres;
    INTVec          hits;
    int             hit_sum;

    if (critical_points) {

        tmp_clique_spheres.clear();

        for (i = 0; i < clique.nodes.size(); i++) {
            index = clique.nodes[i];
            sphere = index % num_spheres;    // correct orinting: sudipto & DTM
            tmp_clique_spheres.push_back(sphere);
        }

        hits.clear();
        hits.resize(receptor_critical_clusters.size(), 0);

        // loop over the clusters
        for (i = 0; i < hits.size(); i++) {

            // loop over the cluster spheres
            for (j = 0; j < receptor_critical_clusters[i].spheres.size(); j++) {

                // loop over the clique spheres
                for (k = 0; k < tmp_clique_spheres.size(); k++) {

                    if (receptor_critical_clusters[i].spheres[j] ==
                        tmp_clique_spheres[k]) {
                        hits[i] = 1;
                        break;
                    }
                }

                if (hits[i] == 1)
                    break;
            }
        }

        hit_sum = 1;

        for (i = 0; i < hits.size(); i++) {
            hit_sum *= hits[i];
        }

        if (hit_sum == 0)
            return false;
        else
            return true;

    } else {
        return true;
    }
}


/************************************************/
// called in match_ligand()

void
Orient::id_all_cliques()
{
    int             i,
                    j,
                    k,
                    size;
    int             next_cand,
                    index;
    // int notcount;
    CLIQUE          tmp_clique;
    double          *new_level_residuals;
    double           tmp_resid;

    int limit_cliques = 0;
    am_iteration_num = 1;
    cliques.clear();

    // init arrays
    size = num_centers * num_spheres;
    cout<<"size is:"<<size<<endl;
    min_nodes = 3;
    max_nodes = 10;
    //cout<<"max_nodes is :" <<max_nodes<<endl; //10
    bool *new_candset = new bool[size * max_nodes];
    bool *new_notset = new bool[size * max_nodes];
    int *new_state = new int[max_nodes];
    new_level_residuals = new double[max_nodes];

    
    // Loop over automated matching loop &&(am_iteration_num < 2)
    //cout<<"cliques.size is:"<<cliques.size()<<endl;
    //cout<<automated_matching<<endl;
    automated_matching = 1;
    //cout<<"am_iteration_num is:"<<am_iteration_num<<endl;
    //cout<<"cliques.size() is:"<<cliques.size()<<endl;
    //cout<<"am_iteration_num is:"<<am_iteration_num<<endl;
    max_orients = 1000;
    orig_tolerance = 0.25;
    //cout<<"min_nodes is:"<<min_nodes<<endl;
    //cout<<"max_nodes is:"<<max_nodes<<endl;
    
    while ((cliques.size() <= max_orients)
           && ((automated_matching) || (am_iteration_num == 1))
           && (am_iteration_num < 10)) {

        tolerance = am_iteration_num * orig_tolerance;
        //cout<<"tolerance is :"<< tolerance<<endl;
        
        for (i = 0; i < size * max_nodes; i++) {
            if (i < size)
                new_candset[i] = true;
            else
                new_candset[i] = false;

            new_notset[i] = false;

            if (i < max_nodes) {
                new_state[i] = -1;
                new_level_residuals[i] = 0.0;
            }
        }

//      
        int new_level = 0;

        // main loop over levels
        while (new_level > -1) {
//            cout << ".";
//            if (new_level <= 0 || new_level >= 10) 
//               cout << endl <<"new_level:" <<new_level << endl;
            // find next true in candset
            next_cand = -1;
            for (i = new_state[new_level] + 1; (i < size) && (next_cand == -1);
                 i++) {
                index = new_level * size + i;

//                if (index > size * max_nodes) cout <<"index:" << index << endl;

                if ((new_candset[index] == true)
                    && (new_notset[index] == false))
                    next_cand = i;

                // compute residuals
                if (next_cand != -1) {

                    // compute residuals
                    if (new_level > 0)
                        new_level_residuals[new_level] =
                            new_level_residuals[new_level - 1];
                    else if (new_level == 0)
                        new_level_residuals[new_level] = 0.0;
                    else 
                        cout << "ERROR new_level is negative." << endl;

                    for (j = new_level - 1; j > -1; j--) {
                        k = next_cand * size + new_state[j];
//                        if (k > num_nodes*num_nodes) cout << "k:" << k << ">" << num_nodes*num_nodes << endl;
                        // sum total resid method
                        new_level_residuals[new_level] += residual_mat[k];
                        //cout<<"new_level_residuals"<<"["<<new_level<<"] is "<<new_level_residuals[new_level]<<endl;

                        // single max resid method
                        // if(residual_mat[k] > new_level_residuals[new_level])
                        // new_level_residuals[new_level] = residual_mat[k];
                    }

                    if (new_level_residuals[new_level] > tolerance) {
                        next_cand = -1;
                        new_candset[index] = false;
                        new_notset[index] = true;
                    }

                }
                // end compute residuals
            }

            // if a candidate node is found at the current level
            if (next_cand != -1) {

                new_state[new_level] = next_cand;
                new_candset[new_level * size + new_state[new_level]] = false;

                new_level++;
                if (new_level < max_nodes)
                    new_state[new_level] = new_state[new_level - 1];
                if (new_level >= max_nodes) { 
                    // add state to cliques if tolerance is proper
                    if ((new_level_residuals[new_level - 1] >
                         orig_tolerance * (am_iteration_num - 1))
                        && (new_level_residuals[new_level - 1] <= tolerance)) {
                        limit_cliques++;

                        tmp_clique.nodes.clear();
                        tmp_clique.residual =
                            new_level_residuals[new_level - 1];
                        for (i = 0; i < new_level; i++)
                            tmp_clique.nodes.push_back(new_state[i]);

                        if (check_clique_critical_points(tmp_clique)){
                                cliques.push_back(tmp_clique);
                                //cout<< "I am in check_clique_chemical_match!"<<endl;
                                // not
			    }
                                

                    }
                    // end clique add code

                    // if you've reached the end of the tree and still have
                    // nodes to add
                    new_level--;
//                    if (new_level <= 0 || new_level >= 10) 
//                       cout << endl <<"new_level:" <<new_level << endl;

                    if (new_level >= 0)
                        new_notset[new_level * size + new_state[new_level]] =
                            true;

                } else {

                    //cout << "=";
                    // recompute candset and notset
                    for (i = 0; i < size; i++) {
                        index = new_state[new_level - 1] * size + i;
                        new_candset[new_level * size + i] =
                            ((new_candset[(new_level - 1) * size + i]));
                        new_notset[new_level * size + i] =
                            ((new_notset[(new_level - 1) * size + i]));

                        // compute residuals
                        if ((new_candset[new_level * size + i])
                            || (new_notset[new_level * size + i])) {

                            tmp_resid = new_level_residuals[new_level - 1];

                            for (j = new_level - 1; j > -1; j--) {
                                k = i * size + new_state[j];

                                // sum total resid method
                                tmp_resid += residual_mat[k];

                                // single max resid method
                                // if(residual_mat[k] > tmp_resid)
                                // tmp_resid = residual_mat[k];

                            }

                            if (tmp_resid > tolerance) {
                                new_candset[new_level * size + i] = false;
                                new_notset[new_level * size + i] = false;
                            }
                        }
                        // End residual comp
                    }

                }

            } else {            // if no candidates are found
                //cout << "a";
                // if notset.size == 0 && level > min_nodes
                //cout<<"min_nodes is:"<<min_nodes<<endl;
                //min_nodes = 3;
                if (new_level >= min_nodes) {
                    // 
                    // notcount = 0;
                    // for(i=0;i<size;i++) // check that notset is empty
                    // if(new_notset[new_level*size+i] == true)
                    // notcount++;
                    // if(notcount == 0) {

                    // add state to cliques
                    if ((new_level_residuals[new_level - 1] >
                         orig_tolerance * (am_iteration_num - 1))
                        && (new_level_residuals[new_level - 1] <= tolerance)) {

                        tmp_clique.nodes.clear();
                        tmp_clique.residual =
                            new_level_residuals[new_level - 1];
                        for (i = 0; i < new_level; i++)
                            tmp_clique.nodes.push_back(new_state[i]);

                        // perform critical point checking
                        if (check_clique_critical_points(tmp_clique))
                                cliques.push_back(tmp_clique);

                    }
                    // end clique add code

                    // }

                }

                new_state[new_level] = -1;
                new_level--;
//                if (new_level <= 0 || new_level >= 10) 
//                   cout << endl <<"new_level:" <<new_level << endl;
            }

            // end if 100X the # of orients are found.  This is mostly for the
            // dense tree cases.
            if (limit_cliques >= 100 * max_orients) {
                cout << "Warning:  Match Search Truncated due to too many " <<
                    max_nodes << " cliques." << endl;
                break;
            }

        }

        am_iteration_num++;

    }                           // End automated matching loop
    
    // sort by residuals
    //cout<<cliques.begin()<<endl;
    sort(cliques.begin(), cliques.end());
    
  
    current_clique = 0;

    // clean up arrays
    delete[]new_candset;
    delete[]new_notset;
    delete[]new_state;
    delete[]new_level_residuals;
    new_candset = NULL;
    new_notset = NULL;
    new_state = NULL;
    new_level_residuals = NULL;
   
}

void 
Orient::match(DOCKMol & mol, const char* recfilename){
    prepare_receptor(recfilename);
    //cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
    //cout<<spheres.size()<<endl;
    int             s1,
                    s2,
                    c1,
                    c2;
    int             i,
                    j,
                    k;
    INTVec          cmt_sphere_idx,
                    cmt_center_idx;
    num_orients = 0;
    last_orient_flag = false;
    //original.clear_molecule();
    centers.clear(); // centers is SphereVec
    clique_spheres.clear(); // XYZVec 
    clique_centers.clear(); // XYZVec
    copy_molecule(original, mol);
    //outlig(original);
    get_centers(original); // save the each atom's x, y, z into the centers
    calculate_ligand_distance_matrix();
    cout<<"num_centers"<<num_centers<<endl;
    cout<<"num_spheres"<<num_spheres<<endl;
    num_nodes = num_spheres * num_centers;
    residual_mat = new double[num_nodes * num_nodes];
    k = 0;
    // populate residual matrix
    for (i = 0; i < num_nodes; i++) {
        s1 = i % num_spheres;
        c1 = i / num_spheres;

        for (j = 0; j < num_nodes; j++) {
            s2 = j % num_spheres;
            c2 = j / num_spheres;

            residual_mat[k] =
                fabs(sph_dist_mat[s1 * num_spheres + s2] -
                lig_dist_mat[c1 * num_centers + c2]);
                //cout<<residual_mat[k]<<endl;

                k++;
            }
        }
    level = 0;
    am_iteration_num = 1;
    id_all_cliques();
    //cout<<"cliques.size()"<<cliques.size()<<endl;
    
}

/************************************************/
void
Orient::new_extract_coords_from_clique(CLIQUE & clique)
{
    // XYZCRD tmp;
    Sphere          tmp;
    int             index,
                    center,
                    sphere;
    int             i;

    clique_spheres.clear();
    clique_centers.clear();
    clique_size = clique.nodes.size();

    for (i = 0; i < clique_size; i++) {
        index = clique.nodes[i];
/*
        sphere = index / num_centers;
        center = index % num_centers;
*/
        // replace it with the proper sphere/center indexing
        sphere = index % num_spheres;
        center = index / num_spheres;
        // DTM - End removal of faulty code - 1/30/07

        tmp = spheres[sphere];
        clique_spheres.push_back(tmp.crds);

        tmp = centers[center];
        clique_centers.push_back(tmp.crds);
    }

}

/************************************************/
void
Orient::calculate_translations()
{

    double sph_com_x = 0.0;
    double sph_com_y = 0.0;
    double sph_com_z = 0.0;
    double cen_com_x = 0.0;
    double cen_com_y = 0.0;
    double cen_com_z = 0.0;

    for (int i = 0; i < clique_size; i++) {
        cen_com_x += clique_centers[i].x;
        cen_com_y += clique_centers[i].y;
        cen_com_z += clique_centers[i].z;
        sph_com_x += clique_spheres[i].x;
        sph_com_y += clique_spheres[i].y;
        sph_com_z += clique_spheres[i].z;
    }

    sph_com_x = sph_com_x / clique_size;
    sph_com_y = sph_com_y / clique_size;
    sph_com_z = sph_com_z / clique_size;
    cen_com_x = cen_com_x / clique_size;
    cen_com_y = cen_com_y / clique_size;
    cen_com_z = cen_com_z / clique_size;

    spheres_com.x = sph_com_x;
    spheres_com.y = sph_com_y;
    spheres_com.z = sph_com_z;
    centers_com.x = cen_com_x;
    centers_com.y = cen_com_y;
    centers_com.z = cen_com_z;
}

/************************************************/
void
Orient::translate_clique_to_origin()
{

    for (int i = 0; i < clique_size; i++) {

        clique_centers[i].x = clique_centers[i].x - centers_com.x;
        clique_centers[i].y = clique_centers[i].y - centers_com.y;
        clique_centers[i].z = clique_centers[i].z - centers_com.z;

        clique_spheres[i].x = clique_spheres[i].x - spheres_com.x;
        clique_spheres[i].y = clique_spheres[i].y - spheres_com.y;
        clique_spheres[i].z = clique_spheres[i].z - spheres_com.z;
    }
}


/************************************************/
void
Orient::calculate_rotation()
{
    int n = clique_size;
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            rotation_matrix[i][j] = 0.0;

    compute_rot_matrix(clique_centers, clique_spheres, rotation_matrix, n);
    
    //for (i = 0; i < 3; i++)
        //for (j = 0; j < 3; j++)
            //cout<<rotation_matrix[i][j]<<endl;
}

/************************************************/
// Called in main loop in dock.cpp.  Is a condition in while loop.

bool
Orient::new_next_orientation(DOCKMol & mol)
{
    //cout << "new_next_orientation" << endl;
    //cout<<"orient_ligand is:"<<orient_ligand<<endl; //0
    //cout<<"last_orient_flag is:"<<last_orient_flag<<endl; //1
   

        // in case no cliques could be found
    if (cliques.size() == 0) {
        cout << "No orients found for current anchor" << endl;
        return false;
    }

	//cout<<"current_clique is:" <<current_clique<<endl;
    new_extract_coords_from_clique(cliques[current_clique]);
    //cout<<"new_extract_coords_from_clique finished!"<<endl;
    calculate_translations();
    //cout<<"calculate_translations finished!"<<endl;
    translate_clique_to_origin();
    calculate_rotation();
    
    copy_molecule(mol, original); // original --> mol, mol is changed
    mol.translate_mol(-centers_com);
    mol.rotate_mol(rotation_matrix);
    mol.translate_mol(spheres_com);

    // DTM - change this line to ensure all cliques are examined as orients (stop skipping the last one) - 1/30/07
        // if ((current_clique == max_orients) || (current_clique == cliques.size() - 1))
    //if ((current_clique == max_orients) || (current_clique == cliques.size())){
	    //if (verbose) cout << "Current clique:" << current_clique << endl;
       //last_orient_flag = true;
    //} else{
       //last_orient_flag = false;
            
        //}

    return true;

   

}


