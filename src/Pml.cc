#include "Pml.hh"

Pml::Pml()
{}

Pml::~Pml()
{}

void Pml::alloc_pml_fields(region_t r)
{

}

void Pml::free_pml_fields()
{

}

void Pml::apply(Face face, Grid &grid)
{
  region_t r = find_face(face, grid);

  switch (face)
  {
  case FRONT:
  case BACK:
    //condition<YZPlane>(r, grid);
    break;

  case LEFT:
  case RIGHT:
    //condition<XZPlane>(r, grid);
    break;

  case TOP:
  case BOTTOM:
    //condition<XYPlane>(r, grid);
    break;
  }
}

void Pml::pml_update_ex()
{

}

void Pml::pml_update_ey()
{

}

void Pml::pml_update_ez()
{

}

void Pml::pml_update_hx()
{
  // depends on h_z_coef1, h_z_coef2
  unsigned int idx; 

  for(i=0,it = bndyptr->i0;i<=bndyptr->dimx;i++,it++)
    for(j=0,jt = bndyptr->j0;j<=bndyptr->dimy-1;j++,jt++)
      for(k=0,kt = bndyptr->k0;k<=bndyptr->dimz-1;k++,kt++)
      {
        idx = pi(i,j,k);
        bndyptr->hxz[i][j][k]= 
          dataptr->h_z_coef1[kt]*
          gridptr->h_coef1[gridptr->index_h_x[it][jt][kt]]*
          bndyptr->hxz[i][j][k]
          +dataptr->h_z_coef2[kt]*
          gridptr->h_z_coef2[gridptr->index_h_x[it][jt][kt]]*
          (gridptr->e_y[it][jt][kt+1]-gridptr->e_y[it][jt][kt]);
        
        bndyptr->hxy[i][j][k] =
          dataptr->h_y_coef1[jt]*
          gridptr->h_coef1[gridptr->index_h_x[it][jt][kt]]*
          bndyptr->hxy[i][j][k]
          -dataptr->h_y_coef2[jt]*
          gridptr->h_y_coef2[gridptr->index_h_x[it][jt][kt]]*
          (gridptr->e_z[it][jt+1][kt]-gridptr->e_z[it][jt][kt]);
        
        gridptr->h_x[it][jt][kt] = 
          (bndyptr->hxz[i][j][k]+bndyptr->hxy[i][j][k]);
      }
}

void Pml::pml_update_hy()
{

}

void Pml::pml_update_hz()
{

}
