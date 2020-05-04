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

#ifndef dealii_mapping_quad9_h
#define dealii_mapping_quad9_h


#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim = dim>
class MappingQuad9 : public MappingQGeneric<dim, spacedim>
{
public:
  MappingQuad9(std::vector<std::vector<Point<spacedim>>> &support_points);

  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;


private:
  virtual std::vector<Point<spacedim>>
  compute_mapping_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell)
    const override;

  std::vector<std::vector<Point<spacedim>>> support_points;
};

DEAL_II_NAMESPACE_CLOSE

#endif
