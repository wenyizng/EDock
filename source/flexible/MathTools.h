/* functions for numpy style matrix manipulation */
#ifndef MATHTOOLS_H
#define MATHTOOLS_H 1
#include <iostream> 
#include <iomanip>
#include <vector>

using namespace std;

/* print the content of vector */
void print_vector(vector<int>& vec);
void print_vector(vector<float>& vec);

/* print the content of matrix */
void print_matrix(vector<vector<float> >& mat);
void print_matrix(vector<vector<int> >& mat);

/* append or generate a list of interger number from start_num till end_num-1 */
void range(vector<int>& vec,int start_num, int end_num);
#endif
