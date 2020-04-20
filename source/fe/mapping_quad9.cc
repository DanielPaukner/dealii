// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2018 by the deal.II authors
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


#include <deal.II/base/derivative_form.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/std_cxx14/memory.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_quad9.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/full_matrix.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <iostream>

DEAL_II_NAMESPACE_OPEN



template <int dim, int spacedim>
MappingQuad9<dim, spacedim>::MappingQuad9(std::map<unsigned int, std::vector<Point<spacedim>>> &map_in)
  :	support_points_map(map_in),
    MappingQGeneric<dim, spacedim>(2)
{
  std::cout << "Constructed MappingQuad9!" << std::endl;


}



template <int dim, int spacedim>
std::unique_ptr<Mapping<dim, spacedim>>
MappingQuad9<dim, spacedim>::clone() const
{
  return std_cxx14::make_unique<MappingQuad9<dim, spacedim>>(*this);
}

template <int dim, int spacedim>
void
MappingQuad9<dim, spacedim>::print()
{
	std::cout << "MappingQuad9 exists!" << std::endl;
}

//template <int dim, int spacedim>
//std::vector<Point<spacedim>>
//MappingQuad9<dim, spacedim>::get_support_points(unsigned int cell_id) const
//{
//	return support_points_map[cell_id];
//}

template <int dim, int spacedim>
std::map<unsigned int, std::vector<Point<spacedim>>>
MappingQuad9<dim, spacedim>::get_support_points_map() const
{
	return support_points_map;
}

//template <int dim, int spacedim>
//std::vector<Point<spacedim>>
//MappingQuad9<dim, spacedim>::get_support_points_for_cell(const unsigned int cell_id) const
//{
//	return support_points_map[cell_id];
//}

template <int dim, int spacedim>
std::vector<Point<spacedim>>
MappingQuad9<dim, spacedim>::compute_mapping_support_points(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell) const
{
  const std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
    vertices = this->get_vertices(cell);

  // number of total points is 9 (3^2) for quad9 element
  // and 27 (3^3) for hex27 element
  unsigned int num_points = std::pow(3, dim);
  std::vector<Point<spacedim>> a(num_points);
  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    a[i] = vertices[i];

  const unsigned int cell_id = cell->active_cell_index();

  // get vector with points for current cell
  std::map<unsigned int, std::vector<Point<spacedim>>> temp = get_support_points_map();

  std::vector<Point<spacedim>> sup_points = temp[cell_id];
  // just copy the points you need, i.e. index 4 to 8
  // this is done to get the correct ordering of the first
  // four vertice which make up the actual cell in the triangulation
  for (unsigned int i = GeometryInfo<dim>::vertices_per_cell; i < num_points; ++i)
	a[i] = sup_points[i];

  // test output
  std::cout << "Computed mapping support points in MappingQuad9!" << std::endl;
  for(const auto &point : a)
  {
	  std::cout << "x: " << point(0) << " y: " << point(1) << " z: " << point(2) << std::endl;
  }

  return a;
}

//---------------------------------------------------------------------------

//
//template <int dim, int spacedim>
//MappingQGeneric<dim, spacedim>
//  StaticMappingQ1<dim, spacedim>::mapping = MappingQGeneric<dim, spacedim>(1);
//


//--------------------------- Explicit instantiations -----------------------
#include "mapping_quad9.inst"


DEAL_II_NAMESPACE_CLOSE
