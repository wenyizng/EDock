/* functions for numpy style matrix manipulation */
#ifndef MATHTOOLS_CPP
#define MATHTOOLS_CPP 1
#include <iostream> 
#include <iomanip>
#include <vector>

using namespace std;

/* print the content of vector */
void print_vector(vector<int>& vec)
{
    int vec_len=vec.size();
    cout<<"[";
    for(int i=0;i<vec_len;i++)
    {
        if (i==(vec_len-1)) cout<<vec[i]<<"]"<<endl;
        else if (i<10 || vec_len-i <10) cout<<vec[i]<<", ";
        else if (i==10 && vec_len>20) cout<<" ... ";
    }
    if (vec_len==0) cout<<"]"<<endl;
}


/* print the content of vector */
void print_vector(vector<float>& vec)
{
    int vec_len=vec.size();
    cout<<"[";
    for(int i=0;i<vec_len;i++)
    {
        if (i==(vec_len-1)) cout<<setprecision(4)<<vec[i]<<"]"<<endl;
        else if (i<10 || vec_len-i <10) cout<<setprecision(4)<<vec[i]<<", ";
        else if (i==10 && vec_len>20) cout<<" ... ";
    }
    if (vec_len==0) cout<<"]"<<endl;
}

/* print the content of matrix */
void print_matrix(vector<vector<float> >& mat)
{
    int m=mat.size();
    cout<<"["<<endl;
    for(int i=0;i<m;i++)
    {
        if (i<10 || m-i <10) print_vector(mat[i]);
        else if (i==10 && m>20) cout<<" ... "<<endl;
    }
    cout<<"]"<<endl;
}

/* print the content of matrix */
void print_matrix(vector<vector<int> >& mat)
{
    int m=mat.size();
    cout<<"["<<endl;
    for(int i=0;i<m;i++)
    {
        if (i<10 || m-i <10) print_vector(mat[i]);
        else if (i==10 && m>20) cout<<" ... "<<endl;
    }
    cout<<"]"<<endl;
}

/* append or generate a list of interger number from start_num till end_num-1 */
void range(vector<int>& vec,int start_num, int end_num)
{
    for(int i=start_num;i<=end_num;i++) vec.push_back(i);
}
#endif
