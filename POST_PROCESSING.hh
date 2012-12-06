#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "silo.h"

#ifndef _POST_PROCESSING_HH_
#define _POST_PROCESSING_HH_

using namespace std;

typedef vector<unsigned int> ui_vector;
typedef vector<int> i_vector;
typedef vector<double> d_vector;

class POST_PROCESSING
{
  public :

    POST_PROCESSING(bool output_dose,unsigned int n_groups,ui_vector* offset,
        d_vector* x,d_vector* y,d_vector* flux_moments,d_vector* scalar_flux,
        d_vector* dose,string* output_file);
    ~POST_PROCESSING();
    /// Create the silo file.
    void Create_silo_file();
    /// Reorder the cells: triangular cells, quadrilateral cells, pentagonal cells, 
    /// etc.
    void Reorder_cells();

  private :

    bool output_dose;
    unsigned int n_groups;
    ui_vector* offset;
    d_vector* x;
    d_vector* y;
    d_vector* flux_moments;
    d_vector* scalar_flux;
    d_vector* dose;
    string* output_file;

    /// Coordinates of the nodes.
    double** coords;
    /// Connectivity.
    i_vector nodelist;
    /// Contains the number of edges for the different types of cells.
    i_vector shapesize;
    /// Contains the number of each types of cells.
    i_vector shapecounts;
    /// Contains the different types of cells.
    i_vector shapetypes;
};

#endif
