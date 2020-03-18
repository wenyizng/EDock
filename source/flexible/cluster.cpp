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

#include "mol.h"
#include "output.h"
#include "cluster.h"
#include "precon.h"

//cluster inputfile cutoff output
//using namespace std;
bool sort_neighboursnum(vector<int> &neighboursnum,vector<int> &arrayindex)
{
    for(int i=0;i<neighboursnum.size()-1;i++)
        for(int j=i+1;j<neighboursnum.size();j++)
        {
            int temp;
            if(neighboursnum[i]<neighboursnum[j])
            {
                temp=neighboursnum[i];
                neighboursnum[i]=neighboursnum[j];
                neighboursnum[j]=temp;
            }
        }
    for(int i=0;i<neighboursnum.size()-1;i++)
        cout<<neighboursnum[i]<<"\t";
    cout<<"\n";
}
int main(int argc, char ** argv)
{
    string input = argv[1];
    const char *filein = input.c_str();

    vector<vector<float> > x,y,z;
    vector<float> energy;
    vector<DOCKMol> clustermols;
    DOCKMol mol;
    float cutoff = atof(argv[2]);
    string output = argv[3];
    readconformations(filein, x,y,z,energy,mol);
    
/*    
    for(int i=0;i<energy.size();i++)
    {
        
        //clustermols.push_back(mol);
        //clustermols[i].current_score=energy[i];
        //cout<<i<<"\t"<<clusters[i].current_score<<endl;
        for(int j=0;j<mol.num_atoms;j++)
        {
            //clustermols[i].x[j]=x[i][j];
            //clustermols[i].y[j]=y[i][j];
            //clustermols[i].z[j]=z[i][j];
	    cout<< "ATOM  " << setw(5) << j+1 << " " << setw(4)<<mol.atom_names[j] << " " << setw(3)<<"ALA" <<"  "<<setw(4) << j+1 <<"   " ;
	    cout<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<x[i][j]<<setw(8)<<y[i][j]<<setw(8)<<z[i][j]<<"\n";
        }
        cout<<"\n";
    }
*/
    //for(int i=0;i<energy.size();i++)
    //   cout<<energy[i]<<"\t";
    //cout<<"\n";
    cluster(x,y,z,energy,mol,cutoff,output);
    
    return true;

}
bool cluster(vector<vector<float> > &x,vector<vector<float> > &y,vector<vector<float> > &z,vector<float> &energy,DOCKMol &mol,float cutoff, string output)
{
     //***************************initial*******************************//
    //
    
    int clusterNum[5] = {1,2,3,4,5};//5 cluster
    int number_cluster=5;
    int clusterindex = 0; 
    int clustermols[5]= {-1,-1,-1,-1,-1};
    //float RMSDcut_original = 2.5; 
    float RMSDcut_original = cutoff;
     
    float RMSDcut = RMSDcut_original;
    float RMSDspacing=0.05;
    float RMSDmax = 5.0;
    float RMSDmin = 1.0;
    int addnum=0;
    int molnum;
    string logname = output+".log";
    if(x.size()<10000)
        molnum=x.size();
    else
        molnum=10000;
    cout<<"molnum is:"<<molnum<<endl; // total mol number for cluster

    vector<bool> removeflag(molnum,false);
    vector<int> clusterflag(molnum,-1);
    vector<vector<float> > rmsdmatrix(molnum,vector<float>(molnum));
    
    
    //caculate the rmsd matrix
    float tempmax=-1;
    float tempmin=100;
    vector<int> neighboursnum(molnum,0);
    int maxneighbour = -1;
    int maxindex=-1;
    //calculate the rmsd matrix and neighboursnum
    for(int i=0;i<molnum;i++)
    {
        //cout<<i<<"\t";
	for(int j=0;j<molnum;j++)
	{
	    //float rmsdtemp = RMSD(x[i],y[i],z[i],x[j],y[j],z[j]);
	    rmsdmatrix[i][j] = RMSD(x[i],y[i],z[i],x[j],y[j],z[j]);
             //if(i==39)
            //cout<<rmsdmatrix[i][j]<<"\t";
            
	    if(rmsdmatrix[i][j]<RMSDcut_original)
	        neighboursnum[i]++;
	    //cout<<rmsdmatrix[i][j]<<endl; 
	    if(rmsdmatrix[i][j]>tempmax)
	    	tempmax=rmsdmatrix[i][j];
	    else if(rmsdmatrix[i][j]<tempmin && rmsdmatrix[i][j]>0)
	    	tempmin=rmsdmatrix[i][j];
            
	}
        //cout<<"\n";
    }
    cout<<"\n";
    cout<<"tempmax is:"<<tempmax<<"\t"<<"tempmin is:"<<tempmin<<endl;
    
    for(int i=0;i<molnum;i++)
    {
        if(neighboursnum[i]>maxneighbour)
	{
	    maxneighbour=neighboursnum[i];
	    maxindex=i;
	}
    }
   
    cout<<"maxneighbour is:"<<maxneighbour<<endl;
    cout<<"maxindex is:"<<maxindex<<endl;
    //for(int i=0;i<x.size();i++)
    //    cout<<rmsdmatrix[maxindex][i]<<endl;
    //calculate the new rmsdcut based on current maxindex
    for(int iteration=0;iteration<10000;iteration++)
    {
        int tempmaxneighbour=0;
        
        for(int j=0;j<molnum;j++)
	{
	    if(rmsdmatrix[maxindex][j]<RMSDcut)
            {
                tempmaxneighbour++;
	    }
	}
        if((tempmaxneighbour*1.0)/molnum > 0.7 && RMSDcut>RMSDmin)
	{          
	    RMSDcut = RMSDcut-RMSDspacing;
	    //cout<<"condation1\t"<<RMSDcut<<endl;
	    continue;
	}
        else if((tempmaxneighbour*1.0)/molnum<0.2 && RMSDcut<RMSDmax)
	{
	     RMSDcut = RMSDcut+RMSDspacing;
             //cout<<"condation2\t"<<RMSDcut<<endl;
	     continue;
	}
        else
        {
             cout<<"tempmaxneighbour:"<<tempmaxneighbour<<"\t"<<tempmaxneighbour*1.0/molnum<<endl;
             break;
        }  
    }
    cout<<"current RMSDcut:"<<RMSDcut<<endl;
    //calculate the new maxindex with the new rmsdcut
    neighboursnum.assign(molnum,0.0);
    
    maxneighbour=-1;
    maxindex=-1;
    for(int i=0;i<molnum;i++)
    {
        for(int j=0;j<molnum;j++)
        {
            if(rmsdmatrix[i][j]<RMSDcut)
	        neighboursnum[i]++;
        }
        if(neighboursnum[i]>maxneighbour)
        {
	    maxneighbour=neighboursnum[i];
	    maxindex=i;
	}
    }
    
    cout<<"maxneighbour is:"<<maxneighbour<<endl;
    cout<<"maxindex is:"<<maxindex<<endl;
    vector<int> clusternum; // record each cluster number
    //cluster 
    for(int i=0;i<number_cluster;i++)
    {
        neighboursnum.assign(molnum,0);
        maxneighbour=-1;
        maxindex=-1;   
        if(addnum>=molnum)
        {
            break;
        }
        
        for(int m=0;m<molnum;m++)
        {
            if(removeflag[m] == true)
                continue;
            for(int k=0;k<molnum;k++)
            {
                if(rmsdmatrix[m][k]<RMSDcut && removeflag[k] == false)
	        neighboursnum[m]++;
            }
            //find the max index in current cluster
            if(neighboursnum[m]>maxneighbour)
            {
	        maxneighbour=neighboursnum[m];
	        maxindex=m;
	    }
        }
        
        //remove the elements in current cluster
        
        clustermols[i]=maxindex;
        cout<<"maxindex is:"<<maxindex<<"\t"<<"maxneighbour is\t"<<maxneighbour<<endl;
        clusternum.push_back(maxneighbour);
        addnum=addnum+maxneighbour;
	cout<<"addnum is:"<<addnum<<endl;
        for(int j=0;j<molnum;j++)
        {
            if(rmsdmatrix[maxindex][j]<RMSDcut && removeflag[j] == false )
            {           
                removeflag[j]=1;
                clusterflag[j]=i+1;
            }
        }
        clusterindex++;
        //continue the next cluster 
    }
    //cout<<clusterindex<<endl;
/*    
    for(int i=0;i<number_cluster;i++)
    {
        int temp=0;
        for(int j=0;j<molnum;j++)
        {
            if(clusterflag[j] == (i+1))
                temp++;
        }
        cout<<clustermols[i]<<":\t"<<temp<<endl;
    }
*/

    //generate the center combo and closc
    //cout<<"79 is the cluster:"<<clusterflag[73]<<" energy "<<energy[73]<<endl;
    /*
    string fp1name=output+"_center.pdb";
    ofstream fp1(fp1name.c_str());
    string fp2name=output+"_combo.pdb";
    ofstream fp2(fp2name.c_str());
    string fp3name=output+"_closc.pdb";
    ofstream fp3(fp3name.c_str());
    */
    string m1name= output+"_center.mol2";
    ofstream m1(m1name.c_str());
    string m2name = output+"_combo.mol2";
    ofstream m2(m2name);
    string m3name= output+"_closc.mol2";
    ofstream m3(m3name.c_str());
    ofstream log(logname.c_str());
    log<<"centerindex"<<"\t"<<"centerenergy"<<"\t"<<"closcindex"<<"\t"<<"closcenergy"<<"\t"<<"clusternum";
    log<<"\n";
    int centermol[clusterindex],closcmol[clusterindex];
    vector<vector<float> > combomolx(clusterindex,vector<float>(mol.num_atoms));
    vector<vector<float> > combomoly(clusterindex,vector<float>(mol.num_atoms));
    vector<vector<float> > combomolz(clusterindex,vector<float>(mol.num_atoms));
    for(int i=0;i<combomolx.size();i++)
        for(int j=0;j<mol.num_atoms;j++)
        {
            combomolx[i][j]=0.0;
            combomoly[i][j]=0.0;
            combomolz[i][j]=0.0;
        }
    //
    for(int i=0;i<clusterindex;i++)
    {
         //center
         centermol[i]=center(clusterflag,rmsdmatrix,i+1);
         
         cout<<"center"<<i<<" index "<<centermol[i]<<"\t"<<"energy "<<energy[centermol[i]]<<endl;
         log<<centermol[i]<<"\t"<<energy[centermol[i]]<<"\t";
         //write center.pdb file and center.mol2
         /*
	 fp1 << "MODEL     "<<setw(4)<<i+1<<endl;
         for(int m=0;m<mol.num_atoms;m++)
	    {
	        fp1 << "ATOM  " << setw(5) << m+1 << " " << setw(4)<<mol.atom_names[m] << " " << setw(3)<<"LIG" <<"  "<<setw(4) << m+1 <<"   " ;
                fp1<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<x[centermol[i]][m]<<setw(8)<<y[centermol[i]][m]<<setw(8)<<z[centermol[i]][m]<<"\n";
            }
        fp1<<"ENDMDL\n";
        */
         
        m1<<"########## Name:"<<"\t"<<i+1<<"\t"<<energy[centermol[i]]<<endl;
        Write_Mol2_2(mol,x[centermol[i]],y[centermol[i]],z[centermol[i]],m1);
        
	
         //combo
         
         combo(i+1,clusterflag,x,y,z, combomolx[i],combomoly[i],combomolz[i]);
         //write center.pdb file and center.mol2
         /*
         fp2 << "MODEL     "<<setw(4)<<i+1<<endl;
         for(int m=0;m<mol.num_atoms;m++)
	 {
	    fp2 << "ATOM  " << setw(5) << m+1 << " " << setw(4)<<mol.atom_names[m] << " " << setw(3)<<"LIG" <<"  "<<setw(4) << m+1 <<"   " ;
            fp2<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<combomolx[i][m]<<setw(8)<<combomoly[i][m]<<setw(8)<<combomolz[i][m]<<"\n";
         }
        fp2<<"ENDMDL\n";
        */
	m2<<"########## Name:"<<"\t"<<i+1<<"\t"<<"0.0"<<endl;
        Write_Mol2_2(mol,combomolx[i],combomoly[i],combomolz[i],m2);
       
        //closc
        closcmol[i] = closec(clusterflag,i+1,combomolx[i],combomoly[i],combomolz[i],x,y,z);
        cout<<"closc"<<i<<" index "<<closcmol[i]<<"\t"<<"energy "<<energy[closcmol[i]]<<"\t"<<"cluster num is "<<clusternum[i]<<endl;
        log<<closcmol[i]<<"\t"<<energy[closcmol[i]]<<"\t"<<clusternum[i]<<endl;
         //write center.pdb file and center.mol2
         /*
         fp3 << "MODEL     "<<setw(4)<<i+1<<endl;
         for(int m=0;m<mol.num_atoms;m++)
	    {
	        fp3 << "ATOM  " << setw(5) << m+1 << " " << setw(4)<<mol.atom_names[m] << " " << setw(3)<<"LIG" <<"  "<<setw(4) << m+1 <<"   " ;
                fp3<<setiosflags(ios::fixed)<<setprecision(3)<<setw(8)<<x[closcmol[i]][m]<<setw(8)<<y[closcmol[i]][m]<<setw(8)<<z[closcmol[i]][m]<<"\n";
            }
        fp3<<"ENDMDL\n";
        */
        m3<<"########## Name:"<<"\t"<<i+1<<"\t"<<energy[closcmol[i]]<<endl;
        Write_Mol2_2(mol,x[closcmol[i]],y[closcmol[i]],z[closcmol[i]],m3);
    }
    //fp1.close();fp2.close();fp3.close();
    m1.close();m2.close();m3.close();log.close();
   
    return true;

}

float RMSD(vector<float> &x1, vector<float> &y1,vector<float> &z1, vector<float> &x2,vector<float> &y2, vector<float> &z2)
{
   float distance2 = 0;
   float rmsd = 0;
   
   for(int i=0;i<x1.size();i++)
   {
	distance2 += (x1[i]-x2[i])*(x1[i]-x2[i])+(y1[i]-y2[i])*(y1[i]-y2[i])+(z1[i]-z2[i])*(z1[i]-z2[i]);
   }
   rmsd = sqrt(distance2/x1.size());
   
   return rmsd;
}

int center(vector<int> &clusterflag,vector<vector<float> > &rmsdmatrix,int cluster)
{
    int centerindex=-1;
    float sumrmsd=1000000;
    double temprmsd = 0;

    for(int i=0;i<clusterflag.size();i++)
    {
        
        if(clusterflag[i]!=cluster)
        {
            continue;
        }
        double temprmsd = 0;
	for(int j=0;j<clusterflag.size();j++)
	{
                if(clusterflag[j]==cluster)
		    temprmsd=temprmsd+rmsdmatrix[i][j];
	}
	if(temprmsd<sumrmsd)
	{
		centerindex=i;
	}
    }
    
    return centerindex;
}

bool combo(int cluster,vector<int> &clusterflag,vector<vector<float> > &x,vector<vector<float> > &y,vector<vector<float> > &z, vector<float> &combomolx,vector<float> &combomoly,vector<float> &combomolz)
{
    int count=0;
    for(int i=0;i<clusterflag.size();i++)
    {
         
        if(clusterflag[i]!=cluster)
            continue;
        count++;
	for(int j=0;j<combomolx.size();j++)
	{
	    combomolx[j]=combomolx[j]+x[i][j];
	    combomoly[j]=combomoly[j]+y[i][j];
	    combomolz[j]=combomolz[j]+z[i][j];
	}
	
    }
    for(int i=0;i<combomolx.size();i++)
    {
	combomolx[i]=combomolx[i]/count;
	combomoly[i]=combomoly[i]/count;
	combomolz[i]=combomolz[i]/count;
    }
    
    return true;
}


int closec(vector<int> & clusterflag,int cluster,vector<float> &combomolx,vector<float> &combomoly,vector<float> &combomolz,vector<vector<float> > &x,vector<vector<float> > &y,vector<vector<float> > &z)
{
    int closecmol;
    float rmsdmin=1000;
    for(int i=0;i<clusterflag.size();i++)
    {
        if(clusterflag[i]!=cluster)
            continue;

	float temprmsd=0;
	temprmsd=combo_rmsd(combomolx,combomoly,combomolz,x[i],y[i],z[i]);
        //cout<<"temprmsd:"<<temprmsd<<endl;
	if(temprmsd<rmsdmin)
	{
	    rmsdmin=temprmsd;
	    closecmol=i;
	}
	
    }
    return closecmol;
}

float combo_rmsd(vector<float> &combomolx,vector<float> &combomoly,vector<float> &combomolz,vector<float> &x,vector<float> &y,vector<float> &z)
{
    float distance2 = 0;
    float rmsd = 0;
    for(int i=0;i<combomolz.size();i++)
    {
	distance2 += (combomolx[i]-x[i])*(combomolx[i]-x[i])+(combomoly[i]-y[i])*(combomoly[i]-y[i])+(combomolz[i]-z[i])*(combomolz[i]-z[i]);
    }
    rmsd = sqrt(distance2/combomolz.size());
    return rmsd;
}


