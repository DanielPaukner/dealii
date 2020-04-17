// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


// read a file in the MSH format used by the GMSH program

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_quad9.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>

#include <istream>
#include <string>

#include "../tests.h"

template <int dim, int spacedim = dim>
void
check_file(const std::string file_name)
{
  Triangulation<dim, spacedim> tria;
  GridIn<dim, spacedim>        gi;
  gi.attach_triangulation(tria);
  
  std::map<unsigned int, std::vector<Point<spacedim>>> map;
  
  std::filebuf fb;
  if (fb.open(file_name, std::ios::in))
    {
      std::istream is(&fb);
      std::string file_ending = file_name.substr(file_name.find_last_of(".") + 1);
      if(file_ending == "msh")
        gi.read_msh(is, map);
      else if(file_ending == "inp")
        gi.read_ucd(is, map, false);
      else
        AssertThrow(false, ExcNotImplemented ());
    }
  
  MappingQGeneric<dim, spacedim> mapping_1(1);
  MappingQuad9<dim, spacedim>    mapping_2(map);
  
  QGaussLobatto<dim> quad_1(2);
  QGaussLobatto<dim> quad_2(3);
  
  FE_Q<dim, spacedim> fe_1(1);
  FE_Q<dim, spacedim> fe_2(2);
  
  UpdateFlags flags = update_quadrature_points;
  
  FEValues<dim, spacedim> fe_values_1(mapping_1, fe_1, quad_1, flags);
  FEValues<dim, spacedim> fe_values_2(mapping_1, fe_2, quad_2, flags);
  FEValues<dim, spacedim> fe_values_3(mapping_2, fe_2, quad_2, flags);
  
  
  for(auto cell : tria.active_cell_iterators())
  {
    // print vertices
    for(unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v )
      deallog << cell->vertex(v) << std::endl;
    deallog << std::endl;
    
    deallog << "MappingQ(1)" << std::endl;
    fe_values_1.reinit(cell);
    for(auto q : fe_values_1.get_quadrature_points ())
        deallog << q << std::endl;
    deallog << std::endl;
    
    deallog << "MappingQ(2)" << std::endl;
    fe_values_2.reinit(cell);
    for(auto q : fe_values_2.get_quadrature_points ())
        deallog << q << std::endl;
    deallog << std::endl;
    
    deallog << "MappingQuad9()" << std::endl;
    fe_values_3.reinit(cell);
    for(auto q : fe_values_3.get_quadrature_points ())
        deallog << q << std::endl;
    deallog << std::endl;
    
  }
  
  
}

int
main()
{
  initlog();

  check_file<2, 3>(std::string(SOURCE_DIR "/quad9/quad9_1ele.msh"));
}
