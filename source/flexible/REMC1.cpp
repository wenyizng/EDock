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
#include <time.h>

#include "energy.h"
#include "MC_simulation.h"
#include "mol.h"
#include "translate_rotate.h"
#include "output.h"
#include "REMC.h"
#include "generaterandnum.h"
//########################################
//amber: VDW forces field parameters,  REmolvec: the replicas, c_nrg/filevdw/fileflexfileflex_drive_tbl/grid_file: for caculation of MC simulation energy, cutoff: for the energy soft(not using), x/y/z/accenergy: for the acceptable conformations coordinates informations and energy, Tmin/Tmax: min and max temparature, outfile: the name of all output files, gen: is the random seed
//########################################
bool REMC(AMBER_TYPER & amber, vector <DOCKMol> & REmolvec, DOCKMol & receptor, Energy_Score & c_nrg,string filevdw, string fileflex,string fileflex_drive_tbl,string grid_file,int cutoff,vector<vector <float> > &x,vector<vector <float> > &y, vector<vector <float> > &z,vector<float> &accenergy,float Tmin, float Tmax,string outfile,mt19937 &gen, bool flexible_flag,vector<vector<int>> biolip_matrix,int swapnumber, int MC_steps,vector<float>sphx,vector<float>sphy,vector<float> sphz,float bindsite_weight,vector<int> protein_atom_index,vector<vector <float> > ave, vector<vector <float> > std)
{
    vector <float> temperature(REmolvec.size(),0); //REmolvec.size() is the replica number
    vector <float> Rotationangle(REmolvec.size(),0);
    vector <float> Rotationangle_initial(REmolvec.size(),0);
    vector<DOCKMol> Tvecmol;
    //Tvecmol saves the conformations which are the MC initial conformation in each replica
    vector<vector<float> > Tvecenergy(REmolvec.size());
    vector<vector<long int>> Tvecsteps(REmolvec.size());
    //vector<vector<float>> plotenergy(REmolvec.size());
    //cout<<"Tvecmol size:"<<Tvecmol.size()<<endl;
    
    vector<float> conformationsacc(REmolvec.size(),0);//total accept probality for each temperature;

    //float Rmax=30.0;
    //float Rmin=2.0;
    int N=REmolvec.size();
    //N is the replica number 
    //cout<<"N size is:"<<N<<"\t"<<"Tvecmol size:"<<Tvecmol.size()<<endl;
    //int swapnum=200;
    int swapnum = swapnumber;
    int ligand_atoms_number=REmolvec[0].num_atoms;
    //vector <float> p(REmolvec.size(),0);
    vector<float> swapacc(N,0);//swap probality for each temperature
    vector<int> swaptimes(N,0);
    vector<int> randnum;
    vector<int> count(N,0);//for the count the replica conformations, if the count[i]==100, the i's replica acceptable conformations into the cluster
    vector <vector<float> >replica_MC_energy(N);
    vector <int> acceprate(REmolvec.size(),0);//record the accepted times

    //###################generate the initial temperature and replica initial conformations############//
    
    for(int i=1;i<=N;i++)
    {
	
	   temperature[i-1]=(Tmin*pow((Tmax/Tmin),(i-1)*1.0/(N-1)));
	   Tvecmol.push_back(REmolvec[i-1]);
	   cout<<i<<"\t"<<temperature[i-1]<<"\t"<<REmolvec[i-1].current_score<<"\t"<<REmolvec[i-1].x[0]<<endl;
	
    }   
    //###################start REMC operation####################################//
    clock_t t1, t2;
    float diff;
    t1 = clock();
    for(int iteration = 0;iteration<swapnum;iteration++)
    {
	
	   cout<<"**************"<<iteration<<"**************"<<endl;	
	   for(int i=0;i<Tvecmol.size();i++)
	   {
    	    cout<<i<<"\t"<<Tvecmol[i].current_score<<endl;
    	    //vector <DOCKMol> MCrotatemols; // the acceptable conformations after the MC
    	    vector <float> tmpenergy; // acceptable conformations energy
            vector<vector<float> > tmpx,tmpy,tmpz; //acceptable conformations x y z coordinates
            vector <int> MCsteps; // the MC step accepting the conformation
            //int method;
            //if(iteration<=swapnum/2)
            //	method=1; // at the begining of the swap, movement big
            //if(iteration>swapnum/2)
            //	method=2; // at the midium of the swap, movement small
            //##########Tvecmol[i][Tvecmol[i].size()-1] means that at each MC, the last one as the conformations
            //if(MC(amber, Tvecmol[i], c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,1,tmpx,tmpy,tmpz, tmpenergy,temperature[i],Rotationangle[i],gen,MCsteps))   
            if(MC(amber, Tvecmol[i], receptor,c_nrg,filevdw.c_str(), fileflex.c_str(),fileflex_drive_tbl.c_str(),grid_file,1,tmpx,tmpy,tmpz, tmpenergy,temperature[i],gen,replica_MC_energy[i],flexible_flag,biolip_matrix,MC_steps,acceprate[i],sphx,sphy,sphz,bindsite_weight,protein_atom_index,ave,std))
            {
                //selection 1
                cout<<i<<"replica acceprate"<<acceprate[i]<<endl;
                cout<<"tmpx size:"<<tmpx.size()<<"\t"<<"tmpenergy size:"<<tmpenergy.size()<<endl;
                cout<<"################################"<<replica_MC_energy[i].size()<<"###############"<<endl;
                for(int k=0;k<tmpenergy.size();k++)
                {
                
                    Tvecenergy[i].push_back(tmpenergy[k]);//for the output 
                    //cout<<"iteration is:"<<iteration<<endl;
                    //cout<<"MCsteps[k] is:"<<MCsteps[k]<<endl;
                    //long int tempsteps=MCsteps[k]+iteration*20000;
                    //cout<<tempsteps<<"\t"<<tmpenergy[k]<<endl;
                    //Tvecsteps[i].push_back(tempsteps);//for the output
                    //plotenergy[i].push_back(tmpenergy[k]);
                    	    
                }
                
        		cout<<"the "<<iteration<<"swap, "<<i<<" acceptor conformations number"<<Tvecenergy[i].size()<<endl;
        		//for the clusters
                //cout<<"x.size()"<<x.size()<<endl;
                //cout<<"writing mol ing......"<<endl;
        		for(int k=0;k<tmpx.size();k++)
	            {
                    /*remcmols.push_back(REmolvec[0]);
                    for(int j=0;j<ligand_atoms_number;j++)
                    {
                    remcmols[remcmols.size()-1].x[j]=tmpx[k][j];
                    remcmols[remcmols.size()-1].y[j]=tmpy[k][j];
                    remcmols[remcmols.size()-1].z[j]=tmpz[k][j];
                    }
                    stringstream ss;
                    ss << k;
                    string str=ss.str();
                    remc<<"########## Name:"<<"\t"<<iteration<<" "<<i<<" "<<k<<"\t"<<tmpenergy[k]<<endl;
                    Write_Mol2(remcmols[remcmols.size()-1],remc,str);
                    */
		            count[i]++;
                    //if(count[i]==50 && i<4 && x.size()<2500)
		            if(count[i]==50 && i<4 && x.size()<2500)
                    {
                	    x.push_back(tmpx[k]);
                	    y.push_back(tmpy[k]);
                	    z.push_back(tmpz[k]);
                	    accenergy.push_back(tmpenergy[k]);
                	    count[i]=0;//acceptable for cluster and re-count the number of i's replica conformations
                    }
			
		        }   
	            
		    
        		//save the last one model for change, code changed by wenyizng 27/11/2017
        		//rm ################Tvecmol[i].push_back(REmolvec[0]);
        		//cout<<"the last energy is:"<<tmpenergy[tmpenergy.size()-1]<<endl;
        		
                //cout<<"I am here"<<endl;
		        for(int k=0;k<ligand_atoms_number;k++)
		        {
        		    Tvecmol[i].x[k]=tmpx[tmpx.size()-1][k];
        		    Tvecmol[i].y[k]=tmpy[tmpy.size()-1][k];
        		    Tvecmol[i].z[k]=tmpz[tmpz.size()-1][k];
		        }
        		Tvecmol[i].current_score=tmpenergy[tmpenergy.size()-1];
        		cout<<i<<" current_score:"<<"\t"<<Tvecmol[i].current_score<<endl;
	        }
            //vector <DOCKMol>().swap(MCrotatemols);
            vector <float>().swap(tmpenergy); // acceptable conformations energy
	        vector<vector<float> >().swap(tmpx);vector<vector<float> >().swap(tmpy); vector<vector<float> >().swap(tmpz);
	        vector <int>().swap(MCsteps); // the MC step accepting the conformation
	        
	    }

        //swap the conformation, calculate the swap rate 
	    int startindex;
	    if(iteration%2 ==0)
	       startindex =0;
	    if(iteration%2==1)
	       startindex = 1;
	
	    //judge exchange the decoys based on temperature
	    cout<<"startindex is:"<<startindex<<endl;
        cout<<"Tvecmol size"<<Tvecmol.size()<<endl;
        float p,rand1,deta;
	    for(int i=startindex;i<Tvecmol.size()-1;i=i+2)
	    {
            
    	    cout<<i<<endl;
    	    float detaenergy = Tvecmol[i+1].current_score - Tvecmol[i].current_score;
    	    cout<<i+1<<" energy:"<<Tvecmol[i+1].current_score<<","<<i <<" energy:"<<Tvecmol[i].current_score<<endl;
    	    cout<<i+1 <<"temperature"<<temperature[i+1]<<","<<i<<"temperature"<<temperature[i]<<endl;
    	    deta = exp((1.0/temperature[i+1]-1.0/temperature[i])*detaenergy);
    	    cout<<"deta is:"<<deta<<endl;
    	    if(deta >= 1.0)
    		  p =1.0;
    	    else
		      p = deta;
            rand1 = rand0to1(gen);
            
            //cout<<"rand1 is "<<rand1<<endl;
            if(rand1<p)
	        {
                cout<<"swap******"<<endl;
                cout<<"p "<<p<<endl;
		        cout<<"rand1 "<<rand1<<endl;
		        swaptimes[i]++;
		        swaptimes[i+1]++;
		        //cout<<"before exchange:"<<endl;
		        //cout<<i+1<<"\t"<<Tvecmol[i+1][Tvecmol[i+1].size()-1].current_score<<"\t"
                //    <<i<<"\t"<<Tvecmol[i][Tvecmol[i].size()-1].current_score<<endl;

		
                float changex[ligand_atoms_number];
                float changey[ligand_atoms_number];
                float changez[ligand_atoms_number];
                
        		float change_energy = Tvecmol[i].current_score;
        		Tvecmol[i].current_score=Tvecmol[i+1].current_score;
        		Tvecmol[i+1].current_score=change_energy;
        		for(int m=0;m<ligand_atoms_number;m++)
		        {
        		    changex[m]=Tvecmol[i].x[m];
        		    Tvecmol[i].x[m] = Tvecmol[i+1].x[m];
        		    Tvecmol[i+1].x[m] = changex[m];

        		    changey[m]=Tvecmol[i].y[m];
        		    Tvecmol[i].y[m] = Tvecmol[i+1].y[m];
        		    Tvecmol[i+1].y[m] = changey[m];

        		    changez[m]=Tvecmol[i].z[m];
        		    Tvecmol[i].z[m] = Tvecmol[i+1].z[m];
        		    Tvecmol[i+1].z[m] = changez[m];
		        }
		        //for output
                // it means there is no conformation acceptable in the last MC simulation
                cout<<"I am here"<<endl;
                if(replica_MC_energy[i].size()==0 && replica_MC_energy[i+1].size()==0)
                {
                    cout<<"checking "<<replica_MC_energy[i].size()<<" "<<replica_MC_energy[i+1].size()<<endl;
                    continue;
       
                } 
                
                if(replica_MC_energy[i].size()==0 && replica_MC_energy[i+1].size()!=0)
                {
                    
                    replica_MC_energy[i+1][replica_MC_energy[i+1].size()-1]=1e+6;
                    continue;
                }
                if(replica_MC_energy[i].size()!=0 && replica_MC_energy[i+1].size()==0)
                {
                    replica_MC_energy[i][replica_MC_energy[i].size()-1]=1e+6;
                    continue;
                }
                
                float replica_change_energy=replica_MC_energy[i][replica_MC_energy[i].size()-1];
                replica_MC_energy[i][replica_MC_energy[i].size()-1]=replica_MC_energy[i+1][replica_MC_energy[i+1].size()-1];
		        replica_MC_energy[i+1][replica_MC_energy[i+1].size()-1]=replica_change_energy;
        		//cout<<"after exchange:"<<endl;
        		//cout<<i+1<<"\t"<<Tvecmol[i+1][Tvecmol[i+1].size()-1].current_score<<"\t"
                //    <<i<<"\t"<<Tvecmol[i][Tvecmol[i].size()-1].current_score<<endl;
		        cout<<"after swap:"<<i+1<<" energy:"<<Tvecmol[i+1].current_score<<","<<i <<" energy:"<<Tvecmol[i].current_score<<endl;
		
	        }
            else
            {
                cout<<"not swap &&&&";
                cout<<rand1<<"\t"<<p<<endl;
            }
	        
	    }
        
        t2 = clock();
        diff = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
        cout<<"time*********"<<diff<<endl;
        if(diff > 3*3600.0)
        //if(diff > 36.0)
        {
        	cout<<"stop!"<<diff<<endl;
        	break;
        }
           
       
    }
    
    //cout<<"Total running time is "<<diff<<" seconds"<<endl;
    printf("\nTotal running time is %5.2f seconds\n", diff); 
    //############## output for checking REMC ###########################//
    string accout=outfile+".acc"; //acceptance file
    ofstream fpacc(accout.c_str());
    fpacc<<"index\ttemperature\tswapnum\tswapacc\tconformations_num\tconformationsacc\tfinal_energy\tave_energy\n";
    //should be changed later
    vector <float> ave_energy(Tvecenergy.size(),0);
    for(int i = 0;i<REmolvec.size();i++)
    {	
    	double tempenergy=0;
    	for(int j=0;j<Tvecenergy[i].size();j++)
	   {
	       tempenergy += Tvecenergy[i][j];
	   } 
	   cout<<i<<" replica size:"<<acceprate[i]<<" total energy:"<<tempenergy<<endl;
	   ave_energy[i] = tempenergy*1.0/Tvecenergy[i].size();
    }
    cout<<"REmolvec size"<<REmolvec.size()<<endl;
    for(int i = 0;i<REmolvec.size();i++)
    {
	
    	//if((i+1)%2==0)
    	//	swapacc[i]=swaptimes[i]*1.0/(swapnum/2.0+1.0);
    	//if((i+1)%2==1)
    	///	swapacc[i]=swaptimes[i]*1.0/(swapnum/2.0);
	    swapacc[i]=swaptimes[i]*1.0/(swapnum);
        conformationsacc[i] = acceprate[i]*1.0/Tvecenergy[i].size();
        //cout<<i<<" re acceprate"<<acceprate[i]<<endl;
        if(acceprate[i] !=0)
        {
            cout<<i<<"\t"<<temperature[i]<<"\t"<<swaptimes[i]<<"\t"<<swapnum<<"\t"<<swapacc[i]<<"\t"<<acceprate[i]<<"\t"<<conformationsacc[i]<<"\t"<<Tvecenergy[i][Tvecenergy[i].size()-1]<<"\t"<<ave_energy[i]<<endl;
            fpacc<<i<<"\t"<<temperature[i]<<"\t"<<swaptimes[i]<<"\t"<<swapacc[i]<<"\t"<<acceprate[i]<<"\t"<<conformationsacc[i]<<"\t"<<Tvecenergy[i][Tvecenergy[i].size()-1]<<"\t"<<ave_energy[i]<<endl;
        }
        
        if(acceprate[i]==0)
        {
            cout<<i<<"\t"<<temperature[i]<<"\t"<<swaptimes[i]<<"\t"<<swapnum<<"\t"<<swapacc[i]<<"\t"<<acceprate[i]<<"\t"<<conformationsacc[i]<<"\t"<<"nan"<<"\t"<<"nan"<<endl;
            fpacc<<i<<"\t"<<temperature[i]<<"\t"<<swaptimes[i]<<"\t"<<swapnum<<"\t"<<swapacc[i]<<"\t"<<acceprate[i]<<"\t"<<conformationsacc[i]<<"\t"<<"nan"<<"\t"<<"nan"<<endl;
        }
	
	}
    fpacc.close();
    string energyout=outfile+".energy";//acceptable energy file
    ofstream fpenergy(energyout.c_str());
    for(int i = 0;i<REmolvec.size();i++)
    {
	   for(int j=0;j<replica_MC_energy[i].size();j++)
	   {
            //fpenergy<<Tvecenergy[i][j]<<" "<<Tvecsteps[i][j]<<","; 
            //fpenergy<<plotenergy[i][j]<<" "<<Tvecsteps[i][j]<<",";
            fpenergy<<replica_MC_energy[i][j]<<"\t";
	   }
	   fpenergy<<"\n";
    }
    fpenergy.close();
    //remc.close();
    return true;
}

void swap_conformations(DOCKMol & mol1, DOCKMol & mol2)
{
    for(int i=0;i<mol1.num_atoms;i++)
    {
	   cout<<mol1.atom_names[i]<<"\t"<<mol2.atom_names[i]<<"\t";
    }
    cout<<"\n"<<"after"<<endl;
    DOCKMol temp;
    copy_molecule(temp,mol1);
    mol1.clear_molecule();
    copy_molecule(mol1,mol2);
    mol2.clear_molecule();
    copy_molecule(mol2,temp);
    temp.clear_molecule();
    for(int i=0;i<mol1.num_atoms;i++)
    {
	   cout<<mol1.atom_names[i]<<"\t"<<mol2.atom_names[i]<<"\t";
    }
    cout<<"\n"<<endl;
}

//num is the random number, maxvalue is the max index;
vector<int> randomgeneratenum(int num,int maxvalue)
{
    vector<int> randnum;
    srand((int)time(0));
    while(randnum.size()<num)
    {
	   int index=int(rand()/(RAND_MAX+1.0)*maxvalue);
	   bool flag=1;
	   for(int i=0;i<randnum.size();i++)
       {
	       if(index==randnum[i])
	       {
		      flag=0;
		      break;
	       }
	   }
	   if(flag==1)
	       randnum.push_back(index);
    }
    
    return randnum;
}
bool sortmol(vector<DOCKMol> &vectormols, vector<int> & vectorindex)
{
    //the ranked results are saved in vector index from low to high
    //cout<<"vectormols.size():"<<vectormols.size()<<endl;
    vector<float> energy(vectormols.size());
    for(int i=0;i<vectormols.size();i++)
    {
        energy[i]=vectormols[i].current_score;
    }
    for(int i=0;i<vectormols.size()-1;i++)
        for(int j=i+1;j<vectormols.size();j++)
        {
            int temp;
            float tempenergy;
            if(energy[j]<energy[i])
            {
                temp=vectorindex[i];
                vectorindex[i]=vectorindex[j];
                vectorindex[j]=temp;
                tempenergy=energy[i];
                energy[i]=energy[j];
                energy[j]=tempenergy;
         
            }
         }

    cout<<"########checking#############"<<endl;
    cout<<"vectorindex.size():"<<vectorindex.size()<<endl;
    for(int i=0;i<vectormols.size();i++)
        cout<<vectorindex[i]<<"\t"<<vectormols[vectorindex[i]].current_score<<",\t";
    cout<<"\n";
    cout<<"########checking#############"<<endl;

    return true;
}
