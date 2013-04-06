#include "POST_PROCESSING.hh"

POST_PROCESSING::POST_PROCESSING(bool output_dose_,unsigned int n_groups_,
    ui_vector* offset_,d_vector* x_,d_vector* y_,d_vector* flux_moments_,
    d_vector* scalar_flux_,d_vector* dose_,string* output_file_) :
  output_dose(output_dose_),
  n_groups(n_groups_),
  offset(offset_),
  c_offset(NULL),
  x(x_),
  y(y_),
  c_x(NULL),
  c_y(NULL),
  flux_moments(flux_moments_),
  c_flux_moments(NULL),
  scalar_flux(scalar_flux_),
  dose(dose_),
  output_file(output_file_)
{
  nodelist.resize(x->size());
}     

POST_PROCESSING::POST_PROCESSING(unsigned int n_groups_,ui_vector* offset_,
    ui_vector* c_offset_,d_vector* x_,d_vector* y_,d_vector* c_x_,d_vector* c_y_,
    d_vector* flux_moments_,d_vector* c_flux_moments_,string* output_file_) :
  output_dose(false),
  n_groups(n_groups_),
  offset(offset_),
  c_offset(c_offset_),
  x(x_),
  y(y_),
  c_x(c_x_),
  c_y(c_y),
  flux_moments(flux_moments_),
  scalar_flux(NULL),
  dose(NULL),
  output_file(output_file_)
{
  nodelist.resize(x->size());
}     

void POST_PROCESSING::Create_transport_silo_file()
{
  DBfile *dbfile = NULL;
  const unsigned int n_moments(flux_moments->size()/x->size());
  const unsigned int n_nodes(x->size());
  const unsigned int n_zones(offset->size()-1);
  const unsigned int n_dims(2);
  /// Coordinates of the nodes.
  double** coords;
  coords = new double*[2];

  // Create the silo file
  dbfile = DBCreate(output_file->c_str(),DB_CLOBBER,DB_LOCAL,
      "Flux moments on unstructured mesh",DB_HDF5);

  Reorder_cells();

  coords[0] = &((*x)[0]);
  coords[1] = &((*y)[0]);

  // Write out connectivity information
  DBPutZonelist2(dbfile,"zonelist",n_zones,n_dims,&nodelist[0],nodelist.size(),0,
      0,0,&shapetypes[0],&shapesize[0],&shapecounts[0],shapetypes.size(),NULL);
  
  // Write an unstructured mesh
  DBPutUcdmesh(dbfile,"mesh",n_dims,NULL,coords,n_nodes,n_zones,"zonelist",NULL,
      DB_DOUBLE,NULL);
  
  delete[] coords;

  // Write the scalar flux
  for (unsigned int g=0; g<n_groups; ++g)
  {
    stringstream group("group_");
    group.seekp(0,ios_base::end);
    group<<g;
    string scalar_flux_str(group.str()+= "_scalar_flux");
    DBPutUcdvar1(dbfile,scalar_flux_str.c_str(),"mesh",&(*scalar_flux)[0],
        n_nodes,NULL,0,DB_DOUBLE,DB_NODECENT,NULL);
  }

  // Write the dose if necessary
  if (output_dose==true)
    DBPutUcdvar1(dbfile,"dose","mesh",&(*dose)[0],n_nodes,NULL,0,DB_DOUBLE,
        DB_NODECENT,NULL);

  // Write the angular flux moments
  for (unsigned int g=0; g<n_groups; ++g)
  {
    stringstream group("group_");
    group.seekp(0,ios_base::end);
    group<<g;
    for (unsigned int i=0; i<n_moments; ++i)
    {
      const unsigned int flux_offset(i*n_nodes+g*n_moments*n_nodes);
      d_vector values(n_nodes);
      stringstream flux("_angular_flux_moment_");
      // Go at the end of the stringstream and append i
      flux.seekp(0,ios_base::end);
      flux<<i;
      string ang_flux_str(group.str()+flux.str());
      for (unsigned int j=0; j<n_nodes; ++j)
        values[j] = (*flux_moments)[flux_offset+j];
      DBPutUcdvar1(dbfile,ang_flux_str.c_str(),"mesh",&values[0],
          n_nodes,NULL,0,DB_DOUBLE,DB_NODECENT,NULL);
    }
  }

  // Close the silo file
  DBClose(dbfile);
}

void POST_PROCESSING::Create_diffusion_silo_file()
{
  DBfile *dbfile = NULL;
  const unsigned int n_nodes(x->size());
  const unsigned int n_zones(offset->size()-1);
  const unsigned int n_dims(2);
  /// Coordinates of the nodes.
  double** coords;
  coords = new double*[2];

  // Create the silo file
  dbfile = DBCreate(output_file->c_str(),DB_CLOBBER,DB_LOCAL,
      "Flux moments on unstructured mesh",DB_HDF5);

  Reorder_cells();

  coords[0] = &((*x)[0]);
  coords[1] = &((*y)[0]);

  // Write out connectivity information
  DBPutZonelist2(dbfile,"zonelist",n_zones,n_dims,&nodelist[0],nodelist.size(),0,
      0,0,&shapetypes[0],&shapesize[0],&shapecounts[0],shapetypes.size(),NULL);
  
  // Write an unstructured mesh
  DBPutUcdmesh(dbfile,"mesh",n_dims,NULL,coords,n_nodes,n_zones,"zonelist",NULL,
      DB_DOUBLE,NULL);

  // Write the scalar flux
  for (unsigned int g=0; g<n_groups; ++g)
  {
    stringstream group("group_");
    group.seekp(0,ios_base::end);
    group<<g;
    string scalar_flux_str(group.str()+= "_scalar_flux");
    DBPutUcdvar1(dbfile,scalar_flux_str.c_str(),"mesh",&(*flux_moments)[0],
        n_nodes,NULL,0,DB_DOUBLE,DB_NODECENT,NULL);
  }

  // Take care of the refined mesh
  Reorder_refined_cells();
  const unsigned int refined_n_nodes(nodelist.size());
  const unsigned int refined_n_zones(n_nodes/3);

  // Concatenate the coordinates and the flux moments
  d_vector* refined_flux_moments = new d_vector();
  d_vector* refined_x = new d_vector();
  d_vector* refined_y = new d_vector();

  refined_flux_moments->reserve(flux_moments->size()+c_flux_moments->size());
  refined_x->reserve(x->size()+c_x->size());
  refined_y->reserve(y->size()+c_y->size());

  refined_flux_moments->insert(refined_flux_moments->end(),flux_moments->begin(),flux_moments->end());
  refined_flux_moments->insert(refined_flux_moments->end(),c_flux_moments->begin(),c_flux_moments->end());
  refined_x->insert(refined_x->end(),x->begin(),x->end());
  refined_x->insert(refined_x->end(),c_x->begin(),c_x->end());
  refined_y->insert(refined_y->end(),y->begin(),y->end());
  refined_y->insert(refined_y->end(),c_y->begin(),c_y->end());

  coords[0] = &((*refined_x)[0]);
  coords[1] = &((*refined_y)[0]);

  // Write out connectivity information
  DBPutZonelist2(dbfile,"refined_zonelist",refined_n_zones,n_dims,&nodelist[0],nodelist.size(),0,
      0,0,&shapetypes[0],&shapesize[0],&shapecounts[0],shapetypes.size(),NULL);
  
  // Write an unstructured mesh
  DBPutUcdmesh(dbfile,"refined_mesh",n_dims,NULL,coords,n_nodes,refined_n_zones,"refined_zonelist",NULL,
      DB_DOUBLE,NULL);

  // Write the scalar flux
  for (unsigned int g=0; g<n_groups; ++g)
  {
    stringstream group("group_");
    group.seekp(0,ios_base::end);
    group<<g;
    string scalar_flux_str(group.str()+= "_refined_scalar_flux");
    DBPutUcdvar1(dbfile,scalar_flux_str.c_str(),"mesh",&(*refined_flux_moments)[0],
        refined_n_nodes,NULL,0,DB_DOUBLE,DB_NODECENT,NULL);
  }


  delete[] coords;

  // Close the silo file
  DBClose(dbfile);
}

void POST_PROCESSING::Reorder_cells()
{
  const unsigned int n_cells(offset->size()-1);
  unsigned int shapecounts_size(0),pos(0);
  vector<i_vector> reorder_offset;
  for (unsigned int i=0; i<n_cells; ++i)
  {
    bool existing_cell_type(false);
    int n_edges((*offset)[i+1]-(*offset)[i]);
    shapecounts_size = shapecounts.size();
    for (unsigned int j=0; j<shapecounts_size; ++j)
      if (shapesize[j]==n_edges)
      {
        reorder_offset[j].push_back((*offset)[i]);
        reorder_offset[j].push_back((*offset)[i+1]);
        shapecounts[j] += 1;
        existing_cell_type = true;
        break;
      }
    if (existing_cell_type==false)
    {
      i_vector cell_offset(2,0);
      cell_offset[0] = (*offset)[i];
      cell_offset[1] = (*offset)[i+1];
      reorder_offset.push_back(cell_offset);
      shapesize.push_back(n_edges);
      shapecounts.push_back(1);
    }
  }
  shapecounts_size = shapecounts.size();
  shapetypes.resize(shapecounts_size);
  for (unsigned int i=0; i<shapecounts_size; ++i)
  {
    if (shapesize[i]==3)
      shapetypes[i] = DB_ZONETYPE_TRIANGLE;
    if (shapesize[i]==4)
      shapetypes[i] = DB_ZONETYPE_QUAD;
    if (shapesize[i]>4)
      shapetypes[i] = DB_ZONETYPE_POLYGON;

    for (int j=0; j<shapecounts[i]; ++j)
    {
      int k_end(reorder_offset[i][2*j+1]);
      for (int k=reorder_offset[i][2*j]; k<k_end; ++k) 
      {
        nodelist[pos] = k;
        ++pos;
      }
    }
  }
}

void POST_PROCESSING::Reorder_refined_cells()
{
  const unsigned int n_cells(offset->size()-1);
  unsigned int n_refined_cells(0);
  i_vector reorder_offset;
  ui_vector n_edges(n_cells,0);
  shapetypes.resize(1, DB_ZONETYPE_TRIANGLE);
  nodelist.clear();

  for (unsigned int i=0; i<n_cells; ++i)
  {
    n_edges[i] = (*offset)[i+1]-(*offset)[i];
    n_refined_cells += n_edges[i];
  }
  
  shapecounts.resize(1,n_refined_cells);

  unsigned int pos(0),k(0),m(n_refined_cells);
  for (unsigned int i=0; i<n_cells; ++i)
  {
    for (unsigned int j=0; j<n_edges[i]; ++j)
    {
      nodelist[pos] = k;
      ++pos;
      nodelist[pos] = k+1;
      ++pos;
      nodelist[pos] = m;
      ++pos;
      ++k;
    }
    ++m;
  }
}
