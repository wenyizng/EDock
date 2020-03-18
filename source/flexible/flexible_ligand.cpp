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
//#include "MC_simulation.h"
//#include "REMC.h"
//#include "precon.h"
//#include "initialmol.h"
#include "flex_dock.h"
using namespace std;

/*example:
flexible_ligand lig outputname
*/

int main(int argc, char ** argv){

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> randgeneration(-1.0, 1.0);
    DOCKMol mol;
    Orient c_orient;
    AMBER_TYPER amber;
    Energy_Score c_nrg;
    char * ligand_in = argv[1];
    string outfile=string(argv[2]);
    string filevdw = "/nfs/amino-home/wenyizng/dock6/dock6/parameters/vdw_AMBER_parm99.defn";
    string fileflex = "/nfs/amino-home/wenyizng/dock6/dock6/parameters/flex.defn";
    string fileflex_drive_tbl = "/nfs/amino-home/wenyizng/dock6/dock6/parameters/flex_drive.tbl";
    float energy_inter_score= 0.0;
    float internal_energy = 0.0 ;


    bool flexible_flag =1; //whethe do flexible docking

    
    /*******************************************************************************************/
    //read the ligand informaiton about the title, mol_info_line, comment1, comment2, comment3,     
    //score_text_data,energy, simplex_text, mol_data, num_atoms, num_bonds, num_residues, and so on.
    //and compulate the ring information
    readlig(mol,ligand_in);
 

    
    //cout<<"********amber initialize****************"<<endl;
    amber.initialize(filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str());
    //cout<<"*******amber.prepare_molecule***********"<<endl;
    // read the vdw_AMBER_parm99.defn, save in latom
    
    amber.prepare_molecule(mol);// prepare the ligand mol for the vdw
    mol.prepare_molecule();//neighbors informations add by wenyi
   
   
    /*************************flexible docking*******************************/
    vector<TORSION> torsions;
    Base_Score score;
    FLOATVec  vertex;//this vector I don't use     
    id_torsions(mol, vertex, torsions); 
    score.initialize_internal_energy(mol);
    score.compute_ligand_internal_energy(mol);//
    //cout<<ligandname<<"\t"<<internal_energy<<"\t"<<heavy_atoms<<"\t"<<bond_num<<endl;
    //if(mol.internal_energy>1.0e-4 && mol.internal_energy<1000.0 && mol.heavy_atoms>1)
    //    cout<<argv[1]<<"\t"<<setprecision(3)<<mol.internal_energy<<"\t"<<mol.heavy_atoms<<"\t"<<torsions.size()<<endl;
    
    
    if(flexible_flag==1)
    {
    
        /******************************************************************************************/
        //torsionkinds = ['1.tor','2.tor','3.tor','32.tor','6.tor','62.tor','63.tor','122.tor'] //
        int torsiontypes[8]={1,2,3,32,6,62,63,122};
        ifstream biolip_file("/nfs/amino-home/wenyizng/dock6/dock6/code/flexible_docking/flexible_ligand_initial/BioLip_ligand_torsion_lib");
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
        cout<<"row number "<<row<<endl;
        cout<<biolip_matrix[0][1]*1.0<<"^^^"<<biolip_matrix[0][0]*1.0<<"&&&&"<<-(log((biolip_matrix[0][1]*1.0)/(biolip_matrix[0][0]*1.0)))<<endl;
  
        
        /******************torsion angle*****************************/
	cout<<"testing testing testing"<<endl;
	//calculate the torsion angle energy
	
	tor_energy_total(amber,biolip_matrix,torsions);
	
        for(int i=0;i<torsions.size();i++)
        {
        //cout<<mol.amber_at_id[i]<<endl;
        //cout<<amber.bond_typer.types[i].drive_id<<"\t"<<torsions[i].ihedral<<endl;
	//cout<<amber.bond_typer.types[amber.bond_typer.flex_ids[torsions[i].bond_num]].drive_id<<'\t'<<torsions[i].ihedral<<endl;
            cout<<"each torsion angle energy "<<tor_energy (amber,biolip_matrix,torsions[i]);
        }
        /*******************************************************************************************/
        
        
        
        score.compute_ligand_internal_energy(mol);
        cout<<"flexible docking part"<<endl;
        cout<<"native flexible energy is:"<<mol.internal_energy<<endl;
        //uniform_int_distribution <> randint(0,torsions.size()-1);
        //for (int i=0;i<100;i++)
        //    cout<<randint(gen)<<" ";
        //cout<<"\n";
        //flexible docking test
         
         /******************************************************************************************/
        string flexiblefile = "flex_" +  outfile + ".mol2";
        ofstream flexsible(flexiblefile.c_str());
        
        for(int i=0;i<torsions.size();i++)
        {
           
            vector<bool> flexibleatoms(mol.num_atoms,0);
            //rotate the Left atoms in bond
            //cout<<"torsions[i].Latoms.size() is:"<<torsions[i].Latoms.size()<<endl;
            //if(torsions[i].Latoms.size()<=1 && torsions[i].Latoms.size() > mol.num_atoms)
            //    continue;
        
           for(int j=0;j<torsions[i].Latoms.size();j++)
        
           {   
               flexibleatoms[torsions[i].Latoms[j]]=1;
           }   
           vector<float> axisA(3,0);
           vector<float> axisB(3,0);
           float angle=180.0*randgeneration(gen); //randnum from-1 to 1
           cout<<"random angle ******* "<<angle<<endl;
           cout<<amber.bond_typer.types[i].drive_id<<endl;
           //cout<<"random number is:"<<randgeneration(gen)<<endl;
           axisA[0]=mol.x[torsions[i].atom2];axisA[1]=mol.y[torsions[i].atom2];axisA[2]=mol.z[torsions[i].atom2];
           axisB[0]=mol.x[torsions[i].atom3];axisB[1]=mol.y[torsions[i].atom3];axisB[2]=mol.z[torsions[i].atom3];
           GroupRotation(axisA, axisB, angle, mol, flexibleatoms);
           float tmp_internal_energy;
           tmp_internal_energy=score.compute_ligand_internal_energy(mol);
           cout<<i<<" flexible internal_energy is:"<<mol.internal_energy<<","<<tmp_internal_energy<<endl;
           flexibleatoms.clear();
           axisA.clear(); axisB.clear();
           //break;
           //rotate the right atoms in bond
          
        
        }
        
        flexsible<<"########## Name:"<<"\t"<<mol.internal_energy<<endl;
        score.compute_ligand_internal_energy(mol);
        cout<<"after flexible ligand, the ligand internal energy is:"<<mol.internal_energy<<endl;
        stringstream ss;
        ss << "_flex_lig";
        string str=string(argv[2])+ss.str();
        Write_Mol2(mol,flexsible,str);
	      flexsible.close();
        	
	
    }
    /////////////////////flexible docking end/////////////////////////////
    
    
   return 1;
}



