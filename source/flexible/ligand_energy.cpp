#include <sstream> 
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <random>

#include "readligand.cpp"
#include "output.h"
#include "match.h"
#include "energy.h"
#include "translate_rotate.h"
#include "MC_simulation.h"
#include "REMC.h"
#include "precon.h"
#include "initialmol.h"
#include "flex_dock.h"
#include "generaterandnum.h"
using namespace std;

/*example:
docking outfile tmin tmax grid sphere 
edock_MC re1 2 100
*/


int main(int argc, char ** argv){

    random_device rd;
    mt19937 gen(rd());
    //uniform_real_distribution<> randgeneration(-1.0, 1.0);
    //uniform_real_distribution<> rand0to1(0, 1.0);
    DOCKMol mol;
    Orient c_orient;
    AMBER_TYPER amber;
    Energy_Score c_nrg;
    
    float cutoff = 100.0;
    //char * smin=argv[2];
    //float Tmin=atof(smin);
    //char * smax=argv[3];
    //float Tmax=atof(smax);
    //string outfile=string(argv[1]);
    bool flexible_flag =1; //whethe do flexible docking
    


    //string grid_file = argv[4];
    //const char * ligand_in = "flex_123.mol2";
    //const char * ligand_in = "lig_charge.mol2";
    char * ligand_in = argv[1];
    //char * receptor_in = argv[5]; // is a .sph file
    string filevdw = "/nfs/amino-home/wenyizng/dock6/dock6/parameters/vdw_AMBER_parm99.defn";
    string fileflex = "/nfs/amino-home/wenyizng/dock6/dock6/parameters/flex.defn";
    string fileflex_drive_tbl = "/nfs/amino-home/wenyizng/dock6/dock6/parameters/flex_drive.tbl";
    //string filevdw = "../../parameters/vdw_AMBER_parm99.defn";
    //string fileflex = "../../parameters/flex.defn";
    //string fileflex_drive_tbl = "../../parameters/flex_drive.tbl";

    int initialnum=10; //REMC number
    int swapnumber = 100;
    int MC_steps = 101;

    //char * vdw_weight_c = argv[6]; //user input vdw weight
    //float vdw_weight = atof(vdw_weight_c);
    //char * bindsite_weight_c = argv[7];//user input binding sites weight
    //float  bindsite_weight = atof(bindsite_weight_c);
    //c_nrg.vdw_scale = vdw_weight;



    /*******************************************************************************************/
    //read the ligand informaiton about the title, mol_info_line, comment1, comment2, comment3,     
    //score_text_data,energy, simplex_text, mol_data, num_atoms, num_bonds, num_residues, and so on.
    //and compulate the ring information
    
    
    cout<<"********read ligand part*****************"<<endl;
    readlig(mol,ligand_in);
    /*******************************************************************************************/

    /*******************************************************************************************/
    //this part read some parameters, if input is native ligand, can output the energy, 
    //this part is used    
    //for train the energy 
    //add by wenyizhang 07052017
    
    
    cout<<"********amber initialize****************"<<endl;
    amber.initialize(filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str());
    cout<<"*******amber.prepare_molecule***********"<<endl;
    // read the vdw_AMBER_parm99.defn, save in latom
    
    amber.prepare_molecule(mol);// prepare the ligand mol for the vdw
    mol.prepare_molecule();//neighbors informations add by wenyi

   


    //flexible part: bond checking
    vector<TORSION> torsions;
    FLOATVec  vertex;//this vector I don't use     
    id_torsions(mol, vertex, torsions); 
    //cout<<"torsions.size is"<<torsions.size()<<endl;
    if (torsions.size() == 0) flexible_flag = 0;
    
    //for (int tt=0;tt<torsions.size();tt++)
    //{
    //    cout<<tt<<"tpye:"<<amber.bond_typer.types[amber.bond_typer.flex_ids[torsions[tt].bond_num]].drive_id<<" ihedral"<<torsions[tt].ihedral<<endl;
    //}
    
    /******idea 13*******/
    //cout<<"ligand number "<<mol.num_atoms<<endl;
    
    if(mol.num_atoms>40 && torsions.size()>10)
    {
        initialnum=20; //REMC number
        swapnumber = 200;
        MC_steps = 201; 
    }
    
    //c_nrg.grid_file_name = grid_file.c_str();
    //c_nrg.initialize(amber);
   
    //read the BioLip_ligand_torsion_lib
    int torsiontypes[8]={1,2,3,32,6,62,63,122};
    ifstream biolip_file("/nfs/amino-home/wenyizng/dock6/dock6/code/flexible_docking/flexible_ligand_initial/BioLip_ligand_torsion_lib");
    //ifstream biolip_file("../../torsion/BioLip_ligand_torsion_lib");
    int row=0;
    int tempnumber=0;
    //int biolip_matrix[8][361]={0};
    vector<vector<int>> biolip_matrix(8,vector<int>(361));
    string pattern = "\t";
    for (string s; getline(biolip_file, s);)
    {
      vector< string> tortimes = split(s, pattern);
      for(int i=0;i<tortimes.size();i++)
      {
        tempnumber = atoi(tortimes[i].c_str());
        biolip_matrix[row][i]=tempnumber;
      }
	    row=row+1;  
		
    }

   //##################################################//
   Base_Score ligscore;
   ligscore.initialize_internal_energy(mol);
   ligscore.compute_ligand_internal_energy (mol);
   cout<<"internal energy "<<mol.internal_energy<<endl;
   mol.energy_tor_total = tor_energy_total (amber,biolip_matrix,torsions);
   cout<<"torsion angles energy "<<mol.energy_tor_total<<endl;
   cout<<"ligand energy:"<<mol.internal_energy + mol.energy_tor_total<<endl;
  return 1;
}



