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
    
    cout << "edock start"<<endl;
    //cout<<"testing "<<largeangle(gen)<<endl;
    //cout << argv[1]<<"\n"<<argv[2]<<endl;
    string grid_file = argv[4];
    //string grid_file = "native";

    //string grid_file = "native";
    //char * ligand_in = argv[2];
    //const char * ligand_in = "crystal_ligand.mol2";
    const char * ligand_in = "lig_charge.mol2";
    char * receptor_in = argv[5]; // is a .sph file
    //const char * receptor_in = "selected_spheres.sph";
    //const char * receptor_in = "select_spheres.sph";

    string filevdw = "../../parameters/vdw_AMBER_parm99.defn";
    string fileflex = "../../parameters/flex.defn";
    string fileflex_drive_tbl = "../../parameters/flex_drive.tbl";
    
    float energy_inter_score= 0.0;
    float internal_energy = 0.0 ;
    //char * cutoffs = argv[4];
    //char * cutoffs =1;
    float cutoff = 100.0;
    char * smin=argv[2];
    float Tmin=atof(smin);
    char * smax=argv[3];
    float Tmax=atof(smax);
    string outfile=string(argv[1]);
    bool flexible_flag =0; //whethe do flexible docking
    int initialnum=40;
    //finding docking center
    vector<float> sphx;
    vector<float> sphy;
    vector<float> sphz; 
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

    //native energy 

    c_nrg.grid_file_name = grid_file.c_str();
    c_nrg.initialize(amber);
    if((energy(amber, mol, c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,sphx,sphy,sphz)))
        
         cout<<"native ligand energy is:"<<mol.current_score<<"\t"<<"internal_energy is:" <<mol.internal_energy<<"\t"<<"intral_energy is:"<<mol.intral_energy<<endl;   

 
    /*************************flexible docking*******************************/
    vector<TORSION> torsions;
    FLOATVec  vertex;//this vector I don't use     
    id_torsions(mol, vertex, torsions); 
    cout<<"torsions.size is"<<torsions.size()<<endl;

    //cout<<"random number is:"<<randgeneration(gen)<<endl;
    /*******************************************************************************************/
    vector<float> centers(3, 0);
    double t[3];
    for(int i=0;i<mol.num_atoms;i++)
        {
            centers[0]+=mol.x[i];
            centers[1]+=mol.y[i];
            centers[2]+=mol.z[i];
        }
    centers[0]=centers[0]/mol.num_atoms;
    centers[1]=centers[1]/mol.num_atoms;
    centers[2]=centers[2]/mol.num_atoms;
  
    /*******************************************************************************************/
    //randomly translate and rotate the ligand
    float tspacing, rspacing;
    tspacing=30.0;rspacing=80.0;
    //translation operator
    double a,b,c;
    vector<float> axis1(3, 0);
    vector<float> axis2(3, 0);
    for(int i=0;i<3;i++)
    {
        //t[i] = ((rand()*1.0/(RAND_MAX*1.0))*2.0-1) * spacing;
        t[i] = randgeneration(gen) * tspacing;
        //cout<<t[i]<<endl;
    }     
    for(int atom=0;atom<mol.num_atoms;atom++)
    {
         mol.x[atom] +=t[0];
         mol.y[atom] +=t[1];
         mol.z[atom] +=t[2];
    }
    
    //rotation operator
    axis1[0]=centers[0]; axis1[1]=centers[1]; axis1[2]=centers[2];
    a = rand0to1(gen);
    b = rand0to1(gen);
    double theta = 2 * 3.1415926 * a;
    double phi = acos(2*b-1.0);
    
    axis2[0] = sin(phi) * cos(theta) + centers[0];
    axis2[1] = sin(phi) * sin(theta) + centers[1];
    axis2[2]= cos(phi) + centers[2];
         
    c=randgeneration(gen);
    float angle = c* rspacing;
    //cout<<"angle is:"<<angle<<endl;
    //cout<<temprotatemols[0].startmol.current_score<<endl;
    GroupRotation(axis1, axis2, angle, mol); // random rotate the mol    
    vector <float>(axis1).swap(axis1);
    vector <float>(axis2).swap(axis2);
    vector <float>(centers).swap(centers);
  
    /*******************************************************************************************/

    cout<<"*******match part***********************"<<endl;
    c_orient.match(mol,receptor_in);    
    cout<<"orient fininshed!!!!!!!"<<endl;
    //cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
    cout<<"spheres.size()"<<c_orient.spheres.size()<<endl;

    
    float spherecenterx,spherecentery,spherecenterz;
    spherecenterx=spherecentery=spherecenterz=0;
    for(int i=0;i<c_orient.spheres.size();i++)
    {
        //cout<<c_orient.spheres[i].crds.x<<"\t"<<c_orient.spheres[i].crds.y<<"\t"<<c_orient.spheres[i].crds.z<<endl;
        spherecenterx += c_orient.spheres[i].crds.x;
        spherecentery += c_orient.spheres[i].crds.y;
        spherecenterz += c_orient.spheres[i].crds.z;
	sphx.push_back(c_orient.spheres[i].crds.x);
        sphy.push_back(c_orient.spheres[i].crds.y);
	sphz.push_back(c_orient.spheres[i].crds.z);
                
	
    }
    spherecenterx=spherecenterx/c_orient.spheres.size();
    spherecentery=spherecentery/c_orient.spheres.size();
    spherecenterz=spherecenterz/c_orient.spheres.size();
    //cout<<spherecenterx<<"\t"<<spherecentery<<"\t"<<spherecenterz<<endl;
    
    //if((energy(amber, mol, c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,sphx,sphy,sphz)))
        
    //     cout<<"native ligand energy is:"<<mol.current_score<<"\t"<<"internal_energy is:" <<mol.internal_energy<<"\t"<<"intral_energy is:"<<mol.intral_energy<<endl;  
    
    vector <DOCKMol> molvec;
    vector <DOCKMol> matchvec;
    /*c_orient.cliques saves the match mols */
    //cout<<c_orient.cliques.size()<<endl;
    int matching_count=0;
    int allmatchcount=0;
    cout<<"c_orient.cliques.size() is:"<<c_orient.cliques.size()<<endl;
    //
    int orient_iterations=0;
    if(c_orient.cliques.size() < 5000 && c_orient.cliques.size()>0) //it means that the mol can be oriented
	orient_iterations=c_orient.cliques.size();
    else if(c_orient.cliques.size() >= 5000)
	orient_iterations=5000;
    for(c_orient.current_clique = 0; c_orient.current_clique<orient_iterations;c_orient.current_clique++) 
    {
	//cout<<"c_orient.current_clique is:"<<c_orient.current_clique<<endl;
        if(c_orient.new_next_orientation(mol))
	    
            if((energy(amber, mol, c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,sphx,sphy,sphz)))
            {
		//cout<<"score is:"<<mol.current_score<<endl;
		if (mol.current_score < 1e+6)
		{
		    molvec.push_back(mol); 
		    if(molvec.size()>1000)
			break;
		    // molvec save the matching conformations
		    //cout<<molvec[molvec.size()-1].current_score<<endl;
		}    
            }
	/*    	    
        {
	    allmatchcount++;
	    //cout<<(energy(amber, mol, c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff))<<endl;
	    cout<<"allmatchcout "<<allmatchcount<<endl;
	    outpdb(mol,allmatchmol,allmatchcount); 
        }
	*/        
    }
    cout<<"the size is:"<<molvec.size()<<endl;



    //allmatchmol.close();   
   bool REMC_flag;
    if(molvec.size()>=initialnum)
        REMC_flag=1;
    else if(molvec.size()==0)
    {
        int num_molvec=molvec.size();
        cout<<"num_molvec "<<num_molvec<<endl;
        vector<vector<float> > ini_x(initialnum-num_molvec, vector<float> (mol.num_atoms,0.0));
        vector<vector<float> > ini_y(initialnum-num_molvec, vector<float> (mol.num_atoms,0.0));
        vector<vector<float> > ini_z(initialnum-num_molvec, vector<float> (mol.num_atoms,0.0));
        REMC_flag = intialmol2(ini_x,ini_y,ini_z,num_molvec, mol, spherecenterx, spherecentery, spherecenterz, initialnum, gen);

        //cout<<"ini_x len"<<ini_x.size()<<endl;
        for(int i=0;i<initialnum-num_molvec;i++)
        {
            DOCKMol moltemp;
            copy_molecule(moltemp,mol);
            for(int atom=0;atom<mol.num_atoms;atom++)
            {
                moltemp.x[atom]=ini_x[i][atom];
                moltemp.y[atom]=ini_y[i][atom];
                moltemp.z[atom]=ini_z[i][atom];
                
            }
            //cout<<ini_x[i][0]<<"\t"<<moltemp.x[0]<<"\t";
            //cout<<molvec[molvec.size()-1].x[0]<<endl;

            if (!(energy(amber, moltemp, c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,sphx,sphy,sphy)))
            {
                moltemp.current_score=1e+6;
                moltemp.internal_energy=1e+6;
                moltemp.intral_energy=1e+6;
            }
            molvec.push_back(moltemp);
            //cout<<molvec[molvec.size()-1].x[0]<<endl;
       }
        
   }
   else if(molvec.size()<initialnum && molvec.size()>0)
   {
       int num_molvec=molvec.size();
       cout<<"num_molvec "<<num_molvec<<endl;
       vector<float> axis1(3, 0);
       vector<float> axis2(3, 0);
       double a,b,c;
       for (int i=0;i<initialnum-num_molvec;i++)
       {   
	   int index=i%num_molvec;
           
           cout<<i<<"\t"<<index<<endl;
           for(int m=0;m<mol.num_atoms;m++)
           {
               mol.x[m]=molvec[index].x[m];
               mol.y[m]=molvec[index].y[m];
               mol.z[m]=molvec[index].z[m];
           }
           DOCKMol moltemp;
           copy_molecule(moltemp,mol);
           molvec.push_back(moltemp);
           cout<<"now molvec.size()"<<molvec.size()<<endl;
           cout<<"before mol"<<molvec[0].x[0]<<"\t"<<molvec[0].y[0]<<"\t"<<molvec[0].z[0]<<endl;

           {
               float ligcenterx,ligcentery,ligcenterz;
               ligcenterx=ligcentery=ligcenterz=0;
               for(int m=0;m<mol.num_atoms;m++)
               {
                   ligcenterx +=mol.x[m];
                   ligcentery +=mol.y[m];
                   ligcenterz +=mol.z[m];
               }
               ligcenterx=ligcenterx/mol.num_atoms;
               ligcentery=ligcentery/mol.num_atoms;
               ligcenterz=ligcenterz/mol.num_atoms;
               cout<<ligcenterx<<"\t"<<ligcentery<<"\t"<<ligcenterz<<endl;
               
               axis1[0]=ligcenterx; axis1[1]=ligcentery; axis1[2]=ligcenterz;
               
               
               a = rand0to1(gen);
               b = rand0to1(gen);
               cout<<"a"<<a<<" "<<"b"<<b<<" "<<"c"<<c<<endl;
               double theta = 2 * 3.1415926 * a;
               double phi = acos(2*b-1.0);
               axis2[0] = sin(phi) * cos(theta) + ligcenterx;
               axis2[1] = sin(phi) * sin(theta) + ligcentery;
               axis2[2]= cos(phi) + ligcenterz;
               c=randgeneration(gen);
               float angle = c*180.0;
               
               GroupRotation(axis2,axis1,angle,mol); // random rotate the mol 
               
                if (!(energy(amber, mol, c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,sphx,sphy,sphz)))
               {
                   mol.current_score=1e+5;
                   mol.internal_energy=1e+5;
                   mol.intral_energy=1e+5;
               }
               for(int m=0;m<mol.num_atoms;m++)
               {
                  molvec[molvec.size()-1].x[m] = mol.x[m];
                  molvec[molvec.size()-1].y[m] = mol.y[m];
                  molvec[molvec.size()-1].z[m] = mol.z[m];

               }
               molvec[molvec.size()-1].current_score=mol.current_score;
               molvec[molvec.size()-1].internal_energy=mol.internal_energy;
               molvec[molvec.size()-1].intral_energy = mol.intral_energy;
               cout<<"after mol"<<mol.x[0]<<"\t"<<mol.y[0]<<"\t"<<mol.z[0]<<endl;
               cout<<"after molvec"<<molvec[molvec.size()-1].x[0]<<"\t"<<molvec[molvec.size()-1].y[0]<<"\t"<<molvec[molvec.size()-1].z[0]<<endl;
               //vector <float>(axis1).swap(axis1);  
               //vector <float>(axis2).swap(axis2);
           }
          
       }
       REMC_flag=1;
   }

    //cout<<"REMC_flag is:"<<REMC_flag<<endl; 
    if(initialnum > molvec.size()) 
        initialnum=molvec.size();
    cout<<"final checking initialnum:"<<initialnum<<endl;
   vector<int> molindex(initialnum);
   vector<DOCKMol> tempmachmols;
   for(int i=0;i<initialnum;i++)
   {
        tempmachmols.push_back(molvec[i]);
        molindex[i]=i;
   }
   sortmol(tempmachmols,molindex);
   /*
   for(int i=0;i<molindex.size();i++)
   {
       cout<<molindex[i]<<"\t";
       //cout<<tempmachmols[molindex[i]].current_score<<",";
   }
   */    
   //cout<<"\n"; 
   for(int i=0;i<initialnum;i++)
   {
       matchvec.push_back(tempmachmols[molindex[i]]);    
   }
   //cout<<"\n";
   ////////////////////////////////////////
   //string matchfile = "match_" + outfile+".pdb";
   //ofstream matchmol(matchfile.c_str());
   //for(int i=0;i<matchvec.size();i++)
   //{
   //   outpdb(matchvec[i],matchmol,i); 
   //}
   //matchmol.close();
   string mol2file="match_" + outfile + ".mol2";
   ofstream mol2(mol2file.c_str());
   for(int i=0;i<matchvec.size();i++)
   {
       Write_Mol2(matchvec[i], mol2, "match");
   }
   mol2.close();
   /*******************************************************************************************/
   //REMC_flag = false;
   if(REMC_flag==true)
   {
       vector<DOCKMol> clustermols;
       vector<vector<float> > x,y,z;
       vector<float> accenergy;
       if(REMC(amber, matchvec, c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,1,x,y,z,accenergy,Tmin,Tmax,outfile,gen,torsions,flexible_flag,sphx,sphy,sphz))
   {
	//string pdbfile = string(argv[1])+".pdb";
	string molfile = string(argv[1])+".mol2";
	//ofstream clusterpdb(pdbfile.c_str());
	ofstream clustermol(molfile.c_str());

        for(int i=0;i<x.size();i++)
        {
	   //cout<<"i="<<i<<endl;
	
	   clustermols.push_back(mol);
	   for(int j=0;j<mol.num_atoms;j++)
	   {
	      //cout<<i<<"\t"<<j<<"\t"<<clustermols[i].atom_names[j]<<endl;
	      //cout<< "ATOM  " << setw(5) << j+1 << " " << setw(4)<<mol.atom_names[j] << " " << setw(3)<<"ALA" <<"  "<<setw(4) << j+1 <<"   " ;
	      //cout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<x[i][j]<<setw(8)<<y[i][j]<<setw(8)<<z[i][j]<<"\n";
	     clustermols[clustermols.size()-1].x[j]=x[i][j];
	     clustermols[clustermols.size()-1].y[j]=y[i][j];
	     clustermols[clustermols.size()-1].z[j]=z[i][j];
	   }
	  //cout<<"\n";
	  clustermols[clustermols.size()-1].current_score=accenergy[i];
	  
	  //outpdb(clustermols[i],clusterpdb,i);
          stringstream ss;
          ss << i;
          string str=string(argv[1])+ss.str();
	  clustermol<<"########## Name:"<<"\t"<<i<<"\t"<<clustermols[i].current_score<<endl;
	  Write_Mol2(clustermols[i],clustermol,str);  
       }
       cout<<"total conformations "<<x.size()<<endl;
       //no conformations can be generated by REMC,using the initial conformations as the final results
       if(x.size()==0)
       {
	   
           for(int i=0;i<matchvec.size();i++)
	   {
	     //outpdb(matchvec[i],clusterpdb,i);
	     stringstream ss;
             ss << i;
             string str=string(argv[1])+ss.str();
	     clustermol<<"########## Name:"<<"\t"<<i<<"\t"<<matchvec[i].current_score<<endl;
	     Write_Mol2(matchvec[i],clustermol,str);  
	   } 
       }
       cout<<"REMC simulation is successful!"<<endl;
       //clusterpdb.close();
       //clustermol.close();	
   }
   else
        {
            cout<<"MC simulation is not successful! Don't cry! fighting!"<<endl;
            REMC_flag=false;
        }
        
  
   }

   
   if(REMC_flag==false)
   {
        //string pdbfile = string(argv[1])+".pdb";
	string molfile = string(argv[1])+".mol2";
	//ofstream clusterpdb(pdbfile.c_str());
	ofstream clustermol(molfile.c_str());
        //outpdb(mol,clusterpdb,1);
        for(int i=0;i<matchvec.size();i++)
        {
            stringstream ss;
            ss << i;
            string str=string(argv[1])+ss.str();
	    clustermol<<"########## Name:"<<"\t"<<i<<"\t"<<matchvec[i].current_score<<endl;
	    Write_Mol2(matchvec[i],clustermol,str);
        }
        
   }
      
   
   return 1;
}



