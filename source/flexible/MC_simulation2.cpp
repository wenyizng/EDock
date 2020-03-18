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

#include "energy.h"
#include "MC_simulation.h"
#include "mol.h"
#include "translate_rotate.h"
#include "output.h"
#include "generaterandnum.h"

//bool MC(AMBER_TYPER & amber, DOCKMol &mol,Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,int cutoff,vector<vector <float> > &x,vector<vector <float> > &y,vector<vector <float> > &z,vector <float> &tmpenergy, float temprature,float angel,mt19937 &gen,vector<int> &MCsteps)
bool MC(AMBER_TYPER & amber, DOCKMol &mol,DOCKMol & receptor, Energy_Score & c_nrg,string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,float cutoff,vector<vector <float> > &x,vector<vector <float> > &y,vector<vector <float> > &z,vector <float> &tmpenergy,float temprature,mt19937 &gen,vector<float> &replica_MC_energy, bool flexible_flag,vector<vector<int>> biolip_matrix,int MC_steps,int & acceprate,vector<float> sphx,vector<float> sphy,vector<float> sphz,float bindsite_weight,vector<int> protein_atom_index,vector<vector <float> > ave, vector<vector <float> > std)
{

    DOCKMol original,endmol;
    copy_molecule(original,mol); // input is mol, copy it to the original
    //int maxsimulations=501; // max simulations times
    int maxsimulations = MC_steps;
    //cout<<"start:temprature is:"<<temprature<<"\t"<<"method is:"<<method<<"\t"<<mol.current_score<<endl;
    cout<<"start:temprature is:"<<temprature<<"\t"<<mol.current_score<<endl;
    int tcout=0;
    copy_molecule(endmol,mol);
    int totaltimes=1000;
    bool flexible_option;
    int tmpacc=0;
    if(flexible_flag==0) 
        flexible_option=0;
    
    DOCKMol tmpmol;
    copy_molecule(tmpmol, mol);
    vector<TORSION> rotatorsions;
    FLOATVec  vertex;//this vector I don't use     
    id_torsions(endmol, vertex, rotatorsions); //
    //int flexible_steps = rotatorsions.size();
    int flexible_steps = 1;
    
    bool movementway = 0; //0 means rigid body, 1 means flexible body
    while(maxsimulations > 1)
    {
    	//cout<<maxsimulations<<"**************************"<<endl;
        //cout<<"totaltimes"<<totaltimes<<"**************************"<<endl;
    	//cout<<"temprotatemols.size()\t"<<temprotatemols.size()<<endl;
    	//cout<<temprotatemols[temprotatemols.size()-1].current_score<<endl;
    	//DOCKMol endmol;
    	//cout<<"MC"<<maxsimulations<<endl;

    	/**************************flexible movement**********************************/
    	//if(flexible_flag == 1 && maxsimulations%2==0)//
        if(flexible_flag ==1) 
        {
            //cout<<"in flexible"<<endl;
            if(totaltimes <0) break;
            movementway = 1;
            for(int mm = 0;mm <flexible_steps; mm++)
            {
                //cout<<"flexible movment testing **********************"<<endl;
                id_torsions(endmol, vertex, rotatorsions); 
                //swapmol(tmpmol,endmol);//copy the current mol to the tmpmol, and flexible tmpmol
                int bondindex=randint(gen,rotatorsions.size());//randomly generate the rotatable bond id
    		    int id_bond = amber.bond_typer.types[amber.bond_typer.flex_ids[rotatorsions[bondindex].bond_num]].drive_id;
        		float flexible_angle=180.0; 
        		float angle;
        		bool type1flag = 0;
        		bool type2flag = 0;
        		bool type3flag = 0;
        		float type2_angle[4] = {-60.0,0,60.0,179.5};
        		float type3_angle[3] = {-60.0,60.0,179.5};
                int tpye2rand;
                int tpye3rand;
                if (maxsimulations < MC_steps / 3.0)
                {
                    
                    switch(id_bond)
                    {
                        case 1:
                            //cout<<"1tpye:"<<amber.bond_typer.types[amber.bond_typer.flex_ids[rotatorsions[bondindex].bond_num]].drive_id<<" ihedral"<<rotatorsions[bondindex].ihedral<<endl;
                            if (abs(rotatorsions[bondindex].ihedral - 180.0)<1.0 && abs(rotatorsions[bondindex].ihedral + 180.0)<1.0)
                                type1flag = 1; 
        			        else
        				        angle = 179.5;
    			        break;
                        case 2:
                            //cout<<"2tpye:"<<amber.bond_typer.types[amber.bond_typer.flex_ids[rotatorsions[bondindex].bond_num]].drive_id<<" ihedral"<<rotatorsions[bondindex].ihedral<<endl;
                            tpye2rand = randint(gen,sizeof(type2_angle));
                            if (abs(rotatorsions[bondindex].ihedral - type2_angle[tpye2rand])>1e-4)
                                angle = type2_angle[tpye2rand]-rotatorsions[bondindex].ihedral;
                            else
                                type2flag = 1;
                            //cout<<"new angle "<<angle<<endl;
                            break;
                        case 3:
                            //cout<<"3tpye:"<<amber.bond_typer.types[amber.bond_typer.flex_ids[rotatorsions[bondindex].bond_num]].drive_id<<" ihedral"<<rotatorsions[bondindex].ihedral<<endl;
    			            tpye3rand = randint(gen,sizeof(type3_angle));
                            if (abs(rotatorsions[bondindex].ihedral - type3_angle[tpye3rand])>1e-4)
                                angle = type3_angle[tpye3rand]-rotatorsions[bondindex].ihedral;
                            else
                                type3flag = 1;
                            //cout<<"new angle "<<angle<<endl;
                        break;
                        default:
                            //cout<<"othertpye:"<<amber.bond_typer.types[amber.bond_typer.flex_ids[rotatorsions[bondindex].bond_num]].drive_id<<" ihedral"<<rotatorsions[bondindex].ihedral<<endl;
                            angle=flexible_angle*randgeneration(gen); //randnum from-1 to 1
    		                //cout<<"random generate angle :"<<angle<<endl;  
                         	
                    }

    		        if ((type1flag ==1 || type2flag ==1 || type3flag == 1 ) && original.current_score < 1e+6) // it means the same angle, not flexible movement
    		        {
                        //cout<<"the same angle in type1 or type2 or type3"<<endl;
                        //cout<<amber.bond_typer.types[amber.bond_typer.flex_ids[rotatorsions[bondindex].bond_num]].drive_id<<" ihedral"<<rotatorsions[bondindex].ihedral<<endl;
                        //no acceptable, record the last one*****.
                        vector <float> tmp(original.num_atoms,0.0);
                        x.push_back(tmp);
                        y.push_back(tmp);
                        z.push_back(tmp);
                        for(int i=0;i<original.num_atoms;i++)
                        {
                          x[x.size()-1][i]=original.x[i];
                          y[y.size()-1][i]=original.y[i];
                          z[z.size()-1][i]=original.z[i];
                        }
                        tmpenergy.push_back(original.current_score );
                        replica_MC_energy.push_back(original.current_score);
                        //swapmol(endmol,original);
                        break;
                    }      
                }
                //else
                {
                    angle=flexible_angle*randgeneration(gen); //randnum from-1 to 1
    		        //cout<< "in random flexible movement"<<endl;   
                }
                float delta_tor_energy = 0.0;
    	        //flexible ligand movement
                //cout<<"before "<<rotatorsions[bondindex].ihedral<<"rotate angle "<<angle<<endl;

                delta_tor_energy = flexible_mol(endmol,rotatorsions[bondindex], angle, amber,biolip_matrix);

                //cout<<"delta_tor_energy"<<delta_tor_energy<<endl;
                //cout<<"after sampling "<<rotatorsions[bondindex].ihedral<<" delta_tor_energy "<<delta_tor_energy<<endl;
                //vector<TORSION> tmptorsions;
                //FLOATVec  vertex1;
                //id_torsions(endmol, vertex1, tmptorsions);
                //cout<<"double check"<<tmptorsions[bondindex].ihedral<<endl;
                
                if (energy(amber, endmol,receptor, c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,biolip_matrix,sphx,sphy,sphz,bindsite_weight,protein_atom_index,ave,std))
                {
                    //swapmol(endmol,tmpmol);
                    //cout<<"mc flexible part:"<<endmol.internal_energy<<" "<<endmol.current_score<<endl;
                    //go to the MC to justy whether the flexible can be acceptable
                    int id = 1;		
                    //cout<<"endmol.internal_energy "<<endmol.internal_energy<<"original.internal_energy "<<original.internal_energy<<endl;

                    //in the docking part, the detaenergy is all of energy. **************
                    
    	            double detaenergy = endmol.current_score - original.current_score; 
    	            //cout<<"dateenergy is:"<<detaenergy<<endl;
                    if(detaenergy>=0)
    	            {
    	                double r = rand0to1(gen);
                        double condation = exp((-detaenergy)/(temprature));
    		            //cout<< "r and condation "<<r<<"\t"<<condation<<endl;
    		            if(r>condation)
		                {
    		              id = 3;
    		              //cout<<"before:"<<endmol.current_score<<endl;

                          //*****no acceptable, record the last one*****/.
                          if(original.current_score < 1e+6)
                          {
                           vector <float> tmp(original.num_atoms,0.0);
                           x.push_back(tmp);
                           y.push_back(tmp);
                           z.push_back(tmp);
                           for(int i=0;i<original.num_atoms;i++)
                           {
                             x[x.size()-1][i]=original.x[i];
                             y[y.size()-1][i]=original.y[i];
                             z[z.size()-1][i]=original.z[i];
                           }
                           tmpenergy.push_back(original.current_score );
                           replica_MC_energy.push_back(original.current_score);  
                          }
                          
    		              swapmol(endmol,original);
		                }
    		            //else
                        //{
                        //    cout<<r<<"\t"<<condation<<"\t"<<endmol.current_score<<endl;
                        //} 
	                }
                    //accept
                    if(id == 1 )
	                {

                        
                        //cout<< "last mol " <<"original.current_score:" <<original.current_score<<" original.intral_energy:" <<original.intral_energy<<" original.internal_energy:"<<original.internal_energy<<" energy_tor_total:"<<original.energy_tor_total<<endl;
                        //cout<<"detaenergy "<<detaenergy<<endl; 
                        /**************add the endmol x y z coordate*************************/
                        if(endmol.current_score<1e+6)
                        {

                            cout<<maxsimulations<< " flexible movment" <<"orginal.current_score:" <<original.current_score<<" original.intral_energy:" <<original.intral_energy<<" original.internal_energy:"<<original.internal_energy<<" energy_tor_total:"<<original.energy_tor_total<<" original.hb_energy:"<<original.hb_energy<< " original.dis_energy:"<<original.dis_energy<<endl;
                            //cout<<"original.current_score"<<original.current_score<<" original.internal_energy:"<<original.internal_energy<<"rotatorsions[bondindex], angle"<<angle<<endl;
                            cout<<maxsimulations<< " flexible movment" <<"endmol.current_score:" <<endmol.current_score<<" endmol.intral_energy:" <<endmol.intral_energy<<" endmol.internal_energy:"<<endmol.internal_energy<<" energy_tor_total:"<<endmol.energy_tor_total<<" endmol.hb_energy:"<<endmol.hb_energy<<"endmol.dis_energy:"<<endmol.dis_energy<<endl;
                            vector <float> tmp(endmol.num_atoms,0.0);
                            x.push_back(tmp);
                            y.push_back(tmp);
                            z.push_back(tmp);
                            for(int i=0;i<endmol.num_atoms;i++)
                            {
                                x[x.size()-1][i]=endmol.x[i];
                                y[y.size()-1][i]=endmol.y[i];
                                z[z.size()-1][i]=endmol.z[i];
                            }
                        
                            tmpenergy.push_back(endmol.current_score );
                            //cout <<"I am in flexible movement accept" <<"flexible accept energy"<<endmol.current_score<<endl;
                            replica_MC_energy.push_back(endmol.current_score);   
                        }
        		        
    	                swapmol(original, endmol);//put endmol in original
                        acceprate ++;	
                        tmpacc++;

	                }
                    //movementway=0; //continue rigid body movement
                    //maxsimulations--; 
                    totaltimes--;

                } //if energy end
                else //can not calculate the energy
                {
                    if(original.current_score < 1e+6)
                    {
                        vector <float> tmp(original.num_atoms,0.0);
                        x.push_back(tmp);
                        y.push_back(tmp);
                        z.push_back(tmp);
                        for(int i=0;i<original.num_atoms;i++)
                        {
                            x[x.size()-1][i]=original.x[i];
                            y[y.size()-1][i]=original.y[i];
                            z[z.size()-1][i]=original.z[i];
                        }
                        tmpenergy.push_back(original.current_score );
                        replica_MC_energy.push_back(original.current_score);
                    }
                    
                    swapmol(endmol,original);
                    totaltimes--;
                }
            } //for end
           
             
        }
    
        /**************************rigid body movement**********************************/ 
        
        int method = 1;     
		if(generatationmol(original,endmol,method,gen))//
	    { 

            if(totaltimes<0) // if 100*1000 times also can finish, break;
            break;
            //cout<<"in rigid"<<endl;
	        if(energy(amber, endmol,receptor,c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,biolip_matrix,sphx,sphy,sphz,bindsite_weight,protein_atom_index,ave,std))
	        {        
    		    int id = 1;
    	        double detaenergy = endmol.current_score - original.current_score;
		        //cout<<"dateenergy is:"<<detaenergy<<endl;
		
	            if(detaenergy>=0)
	            {
	               double r = rand0to1(gen);
                   double condation = exp((-detaenergy)/(temprature));
                   //cout<<r<<"\t"<<condation<<endl;
		           //cout<< "r and condation    "<<r<<"\t"<<condation<<endl;
		           if(r>condation)
		           {
			         id = 3;
			         //cout<<"before:"<<endmol.current_score<<endl;

                     //no acceptable, record the last one
                     if (original.current_score < 1e+6)
                     {
                        vector <float> tmp(original.num_atoms,0.0);
                        x.push_back(tmp);
                        y.push_back(tmp);
                        z.push_back(tmp);
                        for(int i=0;i<original.num_atoms;i++)
                         {
                            x[x.size()-1][i]=original.x[i];
                            y[y.size()-1][i]=original.y[i];
                            z[z.size()-1][i]=original.z[i];
                        }
                        tmpenergy.push_back(original.current_score );
                        replica_MC_energy.push_back(original.current_score);
                     }
                     
			         swapmol(endmol,original);
			         //cout<<"after:"<<endmol.current_score<<endl;
		           }
		          //else
                  //{
                  //    cout<<r<<"\t"<<condation<<"\t"<<endmol.current_score<<endl;
                  //} 
		    
	            }
		        //else if(detaenergy<0)
                //{
                //    cout<<"0.00000"<<"\t"<<"1.0000"<<"\t"<<endmol.current_score<<endl;
                //}
               
	            //accept
		        if(id == 1)
		        {
                    //if(method==2)
                    //{
                    //    cout<<"turn upside down "<<endmol.current_score<<"     accept"<<endl;
                    //}
	                
		            //cout<< "last mol " <<"original.current_score:" <<original.current_score<<" original.intral_energy:" <<original.intral_energy<<" original.internal_energy:"<<original.internal_energy<<" energy_tor_total:"<<original.energy_tor_total<<endl;
                    //add the endmol x y z coordate
            	    //if (endmol.current_score < 1e+6)
                    {
                        //cout<<"endmol.internal_energy "<<endmol.internal_energy<<"original.internal_energy "<<original.internal_energy<<endl;
                        cout<< "rigid movement" <<"original.current_score:"<<original.current_score<<" original.intral_energy:" <<original.intral_energy<<" original.dis_energy:"<<original.dis_energy<<endl;
                        cout<< "rigid movement" <<"endmol.current_score:"<<endmol.current_score<<" endmol.intral_energy:" <<endmol.intral_energy<<" endmol.dis_energy:"<<endmol.dis_energy<<endl;
                        vector <float> tmp(endmol.num_atoms,0.0);
                        x.push_back(tmp);
                        y.push_back(tmp);
                        z.push_back(tmp);
                        for(int i=0;i<endmol.num_atoms;i++)
                        {
                            x[x.size()-1][i]=endmol.x[i];
                            y[y.size()-1][i]=endmol.y[i];
                            z[z.size()-1][i]=endmol.z[i];
                        }
                        tmpenergy.push_back(endmol.current_score);
                        replica_MC_energy.push_back(endmol.current_score);   
                    }
                    
                    swapmol(original, endmol);//put endmol in original
        		    acceprate ++;
                    tmpacc++;
        		    
        		    //for(int j=0;j<temprotatemols[temprotatemols.size()-1].num_atoms;j++)
        	        //    {
        	        //        cout<< "ATOM  " << setw(5) << j+1 << " " << setw(4)<<temprotatemols[temprotatemols.size()-1].atom_names[j] << " " << setw(3)<<"ALA" <<"  "<<setw(4) << j+1 <<"   " ;
        	        //        cout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<temprotatemols[temprotatemols.size()-1].x[j]<<setw(8)<<temprotatemols[temprotatemols.size()-1].y[j]<<setw(8)<<temprotatemols[temprotatemols.size()-1].z[j]<<"\n";
        	        //    }
        		    //cout<<"\n"<<endl;
        		    	
		        }
                maxsimulations--;
                totaltimes--;
	        }
            else // can not calculate the energy 
            {
                //if (original.current_score < 1e+6)
                {
                    vector <float> tmp(original.num_atoms,0.0);
                    x.push_back(tmp);
                    y.push_back(tmp);
                    z.push_back(tmp);
                    for(int i=0;i<original.num_atoms;i++)
                    {
                        x[x.size()-1][i]=original.x[i];
                        y[y.size()-1][i]=original.y[i];
                        z[z.size()-1][i]=original.z[i];
                    }
                    tmpenergy.push_back(original.current_score );
                    replica_MC_energy.push_back(original.current_score);  
                }
               
                swapmol(endmol,original);
                maxsimulations--;
                totaltimes--;
            }
	        continue;    
	    }
        	
    }

    //accepted rate
    //cout<<tmpenergy.size()<<endl;
    //if(tmpenergy.size()>0)
    //    cout<<"end:acceptance "<<tmpenergy.size()<<"\t"<<tmpenergy.size()*1.0/500.0<<" final energy:"<<tmpenergy[tmpenergy.size()-1]<<endl;
    //if(tmpenergy.size()==0)
    //{
    //    cout<<"end:acceptance 0"<<" final energy:"<<mol.current_score<<"\n I want check the program, xsize is "<<x.size()<<endl;
    //    return false;
    //}
    if(x.size() == 0)
	return false;
    //cout<<"temprotatemols.size():"<<temprotatemols.size()<<endl;
    //cout<<"tcout"<<"\t"<<tcout<<endl;

    else
    {
        cout<<"current mc acc times "<<tmpacc<<endl;
        return true;
    }
   
}

bool MC_flexible_ligand(AMBER_TYPER & amber, DOCKMol &mol,float temprature,mt19937 &gen,vector<vector<int>> biolip_matrix)
{

    float current_energy = energy_internal(amber,mol);
    if(current_energy <100.0)
    {
        return true;
    }
    DOCKMol original,endmol;
    energy_internal(amber,mol);//
    copy_molecule(original,mol); // input is mol, copy it to the original
    
    //cout<<"start:temprature is:"<<temprature<<"\t"<<"method is:"<<method<<"\t"<<mol.current_score<<endl;
    cout<<"start:temprature is:"<<temprature<<"\t"<<mol.internal_energy<<endl;
    int tcout=0;
    copy_molecule(endmol,mol);
    DOCKMol tmpmol;
    copy_molecule(tmpmol, mol);
    vector<TORSION> rotatorsions;
    FLOATVec  vertex;//this vector I don't use     
    id_torsions(endmol, vertex, rotatorsions); 
    if (rotatorsions.size() == 0) return false;
    int maxsimulations; // max simulations times
    if (rotatorsions.size()>10)
        maxsimulations=1001;
    else
        maxsimulations=501; 
    while(maxsimulations > 1)
    {
	
	//cout<<"temprotatemols.size()\t"<<temprotatemols.size()<<endl;
	//cout<<temprotatemols[temprotatemols.size()-1].current_score<<endl;
	//DOCKMol endmol;
	//cout<<"MC"<<maxsimulations<<endl;
	/**************************generatationmol**********************************/
        swapmol(tmpmol,endmol);//copy the current mol to the tmpmol, and flexible tmpmol
       
        int bondindex=randint(gen,rotatorsions.size());//randomly generate the rotatable bond id
        float flexible_angle=180.0;
        float angle=flexible_angle*randgeneration(gen); //randnum from-1 to 1
       
        float delta_tor_energy = 0.0;
        //flexible ligand movement
        //cout<<"before "<<rotatorsions[bondindex].ihedral<<"rotate angle "<<angle<<endl;
        
        delta_tor_energy = flexible_mol(tmpmol,rotatorsions[bondindex], angle, amber,biolip_matrix);
        
        //cout<<"after sampling "<<rotatorsions[bondindex].ihedral<<" delta_tor_energy "<<delta_tor_energy<<endl;
        //vector<TORSION> tmptorsions;
        //FLOATVec  vertex1;
        //id_torsions(tmpmol, vertex1, tmptorsions);
        //cout<<"double check"<<tmptorsions[bondindex].ihedral<<endl;
        
        //calculate the torenergy and internal energy
        if (energy_internal(amber,tmpmol))
        {
            swapmol(endmol,tmpmol);
            //cout<<"mc flexible part:"<<endmol.internal_energy<<" "<<endmol.current_score<<endl;
            //go to the MC to justy whether the flexible can be acceptable
            int id = 1;		
            //cout<<"endmol.internal_energy "<<endmol.internal_energy<<"original.internal_energy "<<original.internal_energy<<endl;
	        //double detaenergy = endmol.internal_energy - original.internal_energy + delta_tor_energy;
            double detaenergy = endmol.internal_energy - original.internal_energy; 
            //double detaenergy =delta_tor_energy;
	        //cout<<"dateenergy is:"<<detaenergy<<endl;
            if(detaenergy>0)
	        {
	            double r = rand0to1(gen);
                double condation = exp((-detaenergy)/(temprature));
		        //cout<< "r and condation    "<<r<<"\t"<<condation<<endl;
		        if(r>condation)
		        {
		          id = 3;
		          //cout<<"before:"<<endmol.current_score<<endl;
		          swapmol(endmol,original);
                    
                  //replica_MC_energy.push_back(endmol.current_score);
		          //cout<<"after:"<<endmol.current_score<<endl;
		        }
		        //else
                //{
                //    cout<<r<<"\t"<<condation<<"\t"<<endmol.current_score<<endl;
                //} 
	        }
            //accept
            if(id == 1)
	        {
		      cout<<"acceptor "<<endmol.internal_energy<<endl;
	          swapmol(original, endmol);//put endmol in original	
	        }

            //***********************************************************//  
        }
        maxsimulations--;    
   
    }
    //save the lastes endmol;
    //copy_molecule(mol,endmol);
    //cout<<"before checking"<<mol.x[0]<<" "<<endmol.x[0]<<endl;
    for(int nn=0;nn<endmol.num_atoms;nn++)
    {
        mol.x[nn] = endmol.x[nn];
        //cout<<matchvec[i].x[m]<<endl;
        mol.y[nn] = endmol.y[nn];
        mol.z[nn] = endmol.z[nn];
    }
    mol.current_score=endmol.current_score;
    mol.internal_energy=endmol.internal_energy;
    mol.intral_energy = endmol.intral_energy;
    mol.energy_tor_total = endmol.energy_tor_total;
    //cout<<"after checking"<<mol.x[0]<<" "<<endmol.x[0]<<endl;
    return true;
}

bool MC_rigid_ligand(AMBER_TYPER & amber, DOCKMol &mol, DOCKMol &receptor, Energy_Score & c_nrg, string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file, float cutoff, float temprature,mt19937 &gen,vector<vector<int>> biolip_matrix,vector<float> sphx,vector<float> sphy,vector<float> sphz,float bindsite_weight, vector<int> protein_atom_index,vector<vector <float> > ave, vector<vector <float> > std)
{
    DOCKMol original,endmol;
    copy_molecule(original,mol); // input is mol, copy it to the original
    int maxsimulations=501; // max simulations times
    //cout<<"start:temprature is:"<<temprature<<"\t"<<"method is:"<<method<<"\t"<<mol.current_score<<endl;
    cout<<"start:temprature is:"<<temprature<<"\t"<<mol.current_score<<endl;
    copy_molecule(endmol,mol);
    int totaltimes=1000;
    int tmpacc=0;    
    while(maxsimulations > 1)
    {
        if(totaltimes <0) break;
        int method = 2;     
        if(generatationmol(original,endmol,method,gen))//
        { 

            if(totaltimes<0) // if 100*1000 times also can finish, break;
            break;
            //cout<<"in rigid"<<endl;
            if(energy(amber, endmol,receptor,c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,cutoff,biolip_matrix,sphx,sphy,sphz,bindsite_weight,protein_atom_index,ave,std))
            {        
                int id = 1;
                double detaenergy = endmol.current_score - original.current_score;
                //cout<<"dateenergy is:"<<detaenergy<<endl;
        
                if(detaenergy>=0)
                {
                   double r = rand0to1(gen);
                   double condation = exp((-detaenergy)/(temprature));
                   //cout<<r<<"\t"<<condation<<endl;
                   //cout<< "r and condation    "<<r<<"\t"<<condation<<endl;
                   if(r>condation)
                   {
                     id = 3;
                     //cout<<"before:"<<endmol.current_score<<endl;

                     //no acceptable, record the last one
                     swapmol(endmol,original);
                     //cout<<"after:"<<endmol.current_score<<endl;
                   }
                  //else
                  //{
                  //    cout<<r<<"\t"<<condation<<"\t"<<endmol.current_score<<endl;
                  //} 
            
                }
                //else if(detaenergy<0)
                //{
                //    cout<<"0.00000"<<"\t"<<"1.0000"<<"\t"<<endmol.current_score<<endl;
                //}
               
                //accept
                if(id == 1)
                {
                    cout<< maxsimulations<<" rigid movement" <<"endmol.current_score:"<<endmol.current_score<<" endmol.intral_energy:" <<endmol.intral_energy<<" endmol.internal_energy:"<<endmol.internal_energy<<" endmol.energy_tor_total:"<<endmol.energy_tor_total<<" endmol.hb_energy:"<<endmol.hb_energy<<endl;
                    //cout<< "last mol " <<"original.current_score:" <<original.current_score<<" original.intral_energy:" <<original.intral_energy<<" original.internal_energy:"<<original.internal_energy<<" energy_tor_total:"<<original.energy_tor_total<<endl;
                    swapmol(original, endmol);//put endmol in original
                    tmpacc++;
                    //for(int j=0;j<temprotatemols[temprotatemols.size()-1].num_atoms;j++)
                    //    {
                    //        cout<< "ATOM  " << setw(5) << j+1 << " " << setw(4)<<temprotatemols[temprotatemols.size()-1].atom_names[j] << " " << setw(3)<<"ALA" <<"  "<<setw(4) << j+1 <<"   " ;
                    //        cout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<temprotatemols[temprotatemols.size()-1].x[j]<<setw(8)<<temprotatemols[temprotatemols.size()-1].y[j]<<setw(8)<<temprotatemols[temprotatemols.size()-1].z[j]<<"\n";
                    //    }
                    //cout<<"\n"<<endl;
                        
                }
                maxsimulations--;
                totaltimes--;
            }
            else // can not calculate the energy 
            {
                swapmol(endmol,original);
                maxsimulations--;
                totaltimes--;
            }
            continue;    
        }
          
    }
    for(int nn=0;nn<endmol.num_atoms;nn++)
    {
        mol.x[nn] = endmol.x[nn];
        //cout<<matchvec[i].x[m]<<endl;
        mol.y[nn] = endmol.y[nn];
        mol.z[nn] = endmol.z[nn];
    }
    mol.current_score=endmol.current_score;
    mol.internal_energy=endmol.internal_energy;
    mol.intral_energy = endmol.intral_energy;
    mol.energy_tor_total = endmol.energy_tor_total;
    
}

vector<int> sort_index(RotatemolVec temprotatemols)
{
    vector<int> original_index;
    vector<float> molscores;
    for(int i =0;i<temprotatemols.size()-1;i++)
    {
        original_index.push_back(i);
        molscores.push_back(temprotatemols[i].endmol.current_score);
    }
	


    for(int i=0;i<molscores.size();i++)
    {
	   for(int j=i+1;j<molscores.size();j++)
	   {
	     int temp;
	     float tempscore;
	     if(molscores[j] <molscores[i])
	     {
		 tempscore = molscores[j];
		 molscores[j] = molscores[i];
		 molscores[i] = tempscore;
         temp=original_index[j];
		 original_index[j]=original_index[i];
		 original_index[i]=temp; 
	     }
	     
	   }        
		
    }
    //for(int i=0;i<original_index.size();i++)
    //{
    //    cout<<i<<":"<<original_index[i]<<temprotatemols[original_index[i]].endmol.current_score<<endl;
    //}
    vector <float>(molscores).swap(molscores);
    return original_index;
}


bool generatationmol(DOCKMol & original,DOCKMol & endmol,int method,mt19937 &gen)
{
    
    swapmol(endmol,original); //put original in endmol;
    double a,b,c,t[3],angle;
    vector<float> centers(3, 0);
    for(int i=0;i<endmol.num_atoms;i++)
    {
        centers[0]+=endmol.x[i];
        centers[1]+=endmol.y[i];
        centers[2]+=endmol.z[i];
    }
    centers[0]=centers[0]/endmol.num_atoms;
    centers[1]=centers[1]/endmol.num_atoms;
    centers[2]=centers[2]/endmol.num_atoms;

    float tspacing, rspacing;
    if (method==1)
    {
        tspacing = 0.2; // translate spacing
        rspacing = 3.0;
        c=randgeneration(gen);
        angle = c* rspacing;
        
	
    }
    if (method==2)
    {
        tspacing = 1.0;
        //rspacing = largeangle(gen);
        //angle=rspacing;
        angle = 0.0;
	//cout<<"turn angle " <<angle<<endl;
    
    }
    
    //translation operator
    
    vector<float> axis1(3, 0);
    vector<float> axis2(3, 0);
    for(int i=0;i<3;i++)
    {
        //t[i] = ((rand()*1.0/(RAND_MAX*1.0))*2.0-1) * spacing;
        t[i] = randgeneration(gen) * tspacing;
        //cout<<t[i]<<endl;
    } 
    //cout<<t[0]<<" "<<t[1]<<" "<<t[2]<<endl;      
    for(int atom=0;atom<endmol.num_atoms;atom++)
    {
         endmol.x[atom] +=t[0];
         endmol.y[atom] +=t[1];
         endmol.z[atom] +=t[2];
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
    
    //cout<<"angle is:"<<angle<<endl;
    //cout<<temprotatemols[0].startmol.current_score<<endl;
    GroupRotation(axis1, axis2, angle, endmol); // random rotate the mol    
    vector <float>(axis1).swap(axis1);
    vector <float>(axis2).swap(axis2);
    vector <float>(centers).swap(centers);
    
    //old 
    /*	
    for(int i=0;i<3;i++)
    {
        //a = (rand()*1.0/(RAND_MAX*1.0))*2.0-1;
        //b = (rand()*1.0/(RAND_MAX*1.0))*2.0-1;
        a=randgeneration(gen);
        b=randgeneration(gen);
        axis1[i] = a+centers[i];
        axis2[i] = b+centers[i];  
    }
    */
    //new
    return true;
}

float flexible_mol(DOCKMol & mol, TORSION &torsion, float angle, AMBER_TYPER & amber,vector<vector<int>> biolip_matrix)
{
    float delta_tor_energy=0.0;
    float old_tor_energy=0.0;
    float new_tor_energy=0.0;
    //before rotating the bond, the torsion energy
    old_tor_energy = tor_energy (amber,biolip_matrix,torsion);
    //cout<<"old_tor_energy "<<old_tor_energy<<endl;
    //cout<<"old torsion angle "<<torsion.ihedral<<endl;
    vector<bool> flexibleatoms(mol.num_atoms,0);
    
    //int bondindex=randint(gen,torsions.size());
    //cout<<"bondindex is"<<bondindex<<endl;  
    
    if (torsion.Latoms.size()>1)
    {
         for(int j=0;j<torsion.Latoms.size();j++)   
        {   
             flexibleatoms[torsion.Latoms[j]]=1;
             //cout<<" atom name "<<mol.atom_names[j];

        }
    }  
    else
    {
         for(int j=0;j<torsion.Ratoms.size();j++)   
         {   
             flexibleatoms[torsion.Ratoms[j]]=1;
             //cout<<" atom name "<<mol.atom_names[j];

         }
    }  
    

   
    vector<float> axisA(3,0);
    vector<float> axisB(3,0);
    axisA[0]=mol.x[torsion.atom2];axisA[1]=mol.y[torsion.atom2];
    axisA[2]=mol.z[torsion.atom2];
    axisB[0]=mol.x[torsion.atom3];axisB[1]=mol.y[torsion.atom3];
    axisB[2]=mol.z[torsion.atom3];
    GroupRotation(axisA, axisB, angle, mol, flexibleatoms);//rotate angle 
    
    //compulate the 2Dihedral for a1 a2 a3 a4
    vector<float> c1(3,0);
    vector<float> c2(3,0);
    vector<float> c3(3,0);
    vector<float> c4(3,0);
    c1[0]=mol.x[torsion.atom1]; c1[1]=mol.y[torsion.atom1]; c1[2]=mol.z[torsion.atom1];
    c2[0]=mol.x[torsion.atom2]; c2[1]=mol.y[torsion.atom2]; c2[2]=mol.z[torsion.atom2];
    c3[0]=mol.x[torsion.atom3]; c3[1]=mol.y[torsion.atom3]; c3[2]=mol.z[torsion.atom3];
    c4[0]=mol.x[torsion.atom4]; c4[1]=mol.y[torsion.atom4]; c4[2]=mol.z[torsion.atom4];
    //cout<<"2dihedralis:"<<Points2Dihedral(c1,c2,c3,c4)<<endl;
    torsion.ihedral=Points2Dihedral(c1,c2,c3,c4);
    if((torsion.ihedral-(-180.0))<=1e-4)
	torsion.ihedral=180.0;
    //cout<<"flexible_mol "<<torsion.ihedral<<" angle "<<angle<<endl;
    vector <float> (c1).swap(c1);vector <float> (c2).swap(c2);vector <float> (c3).swap(c3);vector <float> (c4).swap(c4);
    vector <float>(axisA).swap(axisA);
    vector <float>(axisB).swap(axisB);
    vector <bool>(flexibleatoms).swap(flexibleatoms); 
    //after rotating the bond, the torsion energy
    new_tor_energy = tor_energy (amber,biolip_matrix,torsion);
    //cout<<"new_tor_energy "<<new_tor_energy<<" old_tor_energy "<<old_tor_energy<<endl;
    delta_tor_energy = new_tor_energy - old_tor_energy;
    return delta_tor_energy;
}

void swapmol(DOCKMol & mol1, DOCKMol & mol2)//put mol2 in mol1
{
    for(int i=0;i<mol1.num_atoms;i++)
    {
	mol1.x[i]=mol2.x[i];
	mol1.y[i]=mol2.y[i];
	mol1.z[i]=mol2.z[i];	
    }
    mol1.current_score=mol2.current_score;
    mol1.internal_energy=mol2.internal_energy;
    mol1.intral_energy=mol2.intral_energy;
    mol1.energy_tor_total=mol2.energy_tor_total;
    mol1.hb_energy = mol2.hb_energy;
    mol1.contact_energy = mol2.contact_energy;
    mol1.dis_energy = mol2.dis_energy;
}












