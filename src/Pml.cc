#include "Pml.hh"

Pml::Pml()
{}

Pml::~Pml()
{}

void Pml::alloc_pml_fields()
{
  if (thickness_ == 0)
    throw exception();
  
  region_t r = find_face(face, grid);

  unsigned int sz = (r.xmax - r.xmin) * (r.ymax - r.ymin) 
    * (r.zmax - r.zmin);
  
  exy_ = new field_t[sz];
  exz_ = new field_t[sz];

  eyx_ = new field_t[sz];
  eyz_ = new field_t[sz];

  ezx_ = new field_t[sz];
  ezy_ = new field_t[sz];

  hxy_ = new field_t[sz];
  hxz_ = new field_t[sz];

  hyx_ = new field_t[sz];
  hyz_ = new field_t[sz];

  hzx_ = new field_t[sz];
  hzy_ = new field_t[sz];
}

void Pml::set_thickness(unsigned int thickness)
{
  thickness_ = thickness;
}

void Pml::free_pml_fields()
{
  if (exy_) 
    delete[] exy_;

  if (exz_) 
    delete[] exz_;

  if (eyx_) 
    delete[] eyx_;

  if (eyz_) 
    delete[] eyz_;

  if (ezx_) 
    delete[] ezx_;

  if (ezy_) 
    delete[] ezy_;

  if (hxy_) 
    delete[] hxy_;

  if (hxz_) 
    delete[] hxz_;

  if (hyx_) 
    delete[] hyx_;

  if (hyz_) 
    delete[] hyz_;

  if (hzx_) 
    delete[] hzx_;

  if (hzy_) 
    delete[] hzy_;
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
