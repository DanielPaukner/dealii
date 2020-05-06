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

#ifndef dealii_mapping_q2_h
#define dealii_mapping_q2_h

#include <deal.II/base/config.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup mapping */
/*@{*/


/**
 * Implementation of a $d$-linear mapping from the reference cell to a general
 * quadrilateral/hexahedron.
 *
 * The mapping implemented by this class maps the reference (unit) cell to a
 * general grid cell with straight lines in $d$ dimensions. (Note, however,
 * that in 3D the <i>faces</i> of a general, trilinearly mapped cell may be
 * curved, even if the edges are not). This is the standard mapping used for
 * polyhedral domains. It is also the mapping used throughout deal.II for many
 * functions that come in two variants, one that allows to pass a mapping
 * argument explicitly and one that simply falls back to the MappingQ1 class
 * declared here. (Or, in fact, to an object of kind MappingQGeneric(1), which
 * implements exactly the functionality of this class.)
 *
 * The shape functions for this mapping are the same as for the finite element
 * FE_Q of polynomial degree 1. Therefore, coupling these two yields an
 * isoparametric element.
 *
 * @note This class is, in all reality, nothing more than a different name for
 * calling MappingQGeneric with a polynomial degree of one as argument.
 *
 * @author Guido Kanschat, 2000, 2001; Ralf Hartmann, 2000, 2001, 2005,
 * Wolfgang Bangerth, 2015
 */

template <int dim, int spacedim = dim>
class MappingQ2 : public MappingQGeneric<dim, spacedim>
{
public:
  /**
   * Constructor.
   *
   * @param[in]
   */
  MappingQ2(std::vector<std::vector<Point<spacedim>>> &support_points);

  /**
   * Return a pointer to a copy of the present object. The caller of this copy
   * then assumes ownership of it.
   */
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;


private:
  virtual std::vector<Point<spacedim>>
  compute_mapping_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell)
    const override;

  std::vector<std::vector<Point<spacedim>>> support_points;
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif
