#include "POST_PROCESSING.hh"

POST_PROCESSING::POST_PROCESSING(ui_vector* offset_,d_vector* x_,d_vector* y_,
    d_vector* flux_moments_,string* output_file_) :
  offset(offset_),x(x_),y(y_),flux_moments(flux_moments_),output_file(output_file_)
{
  nodelist.resize(x->size());
  coords = new double*[2];
}     

POST_PROCESSING::~POST_PROCESSING()
{
  delete[] coords;
}

void POST_PROCESSING::Create_silo_file()
{
  DBfile *dbfile = NULL;
  const unsigned int n_moments(flux_moments->size()/x->size());
  const unsigned int n_nodes(x->size());
  const unsigned int n_zones(offset->size()-1);
  const unsigned int n_dims(2);

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

  for (unsigned int i=0; i<n_moments; ++i)
  {
    const unsigned int flux_offset(i*n_nodes);
    d_vector values(n_nodes);
    stringstream flux("flux_moment_");
    // Go at the end of the stringstream and append i
    flux.seekp(0,ios_base::end);
    flux<<i;

    for (unsigned int j=0; j<n_nodes; ++j)
      values[j] = (*flux_moments)[flux_offset+j];
    DBPutUcdvar1(dbfile,flux.str().c_str(),"mesh",&values[0],n_nodes,NULL,0,
        DB_DOUBLE,DB_NODECENT,NULL);
  }

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
