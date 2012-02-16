/*
 * This program creates an silo file for unstructured meshes. The program
 * needs to two arguments: the name of the input file and the name of the
 * output file. The input file must have the following structure: number of
 * cells, numbers of nodes, numbers of values of flux moments (multiple of the
 * number of nodes), the offset vector (position of the first node of cell 1,
 * position of the first node of cell 2, ..., number of nodes), coordinates (x,y) 
 * of the nodes, value of the flux moments at each nodes.
 */

#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include "POST_PROCESSING.hh"

using namespace std;

typedef vector<unsigned int> ui_vector;
typedef vector<double> d_vector;

int main(int argc,char **argv)
{
  assert(argc==3);
  unsigned int n_cells,n_nodes,n_values;
  ui_vector offset;
  d_vector x,y,flux_moments;
  string output_file(argv[2]);
  output_file.append(".silo");

  cout<<"Start reading the input file."<<endl;

  // Open the file to read it
  ifstream file(argv[1],ios::in);

  // Read the number of cells, the number of nodes and the number of values
  file>>n_cells>>n_nodes>>n_values;

  // Resize the vectors
  offset.resize(n_cells+1,0);
  x.resize(n_nodes,0.);
  y.resize(n_nodes,0.);
  flux_moments.resize(n_values,0.);

  // Read the offset
  for (unsigned int i=0; i<n_cells+1; ++i)
    file>>offset[i];
  // Read the coordinates of the nodes
  for (unsigned int i=0; i<n_nodes; ++i)
    file>>x[i]>>y[i];
  // Read the value of the flux moments
  for (unsigned int i=0; i<n_values; ++i)
    file>>flux_moments[i];

  // Close the file
  file.close();
  
  cout<<"Done reading the input file."<<endl;

  POST_PROCESSING post_process(&offset,&x,&y,&flux_moments,&output_file);

  cout<<"Start creating the silo file."<<endl;
  post_process.Create_silo_file();
  cout<<"Done creating the silo file."<<endl;

  return 0;
}
