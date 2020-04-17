/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, University of Heidelberg, 1999
 */


// @sect3{Include files}

// Again, the first few include files are already known, so we won't comment
// on them:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

// This one is new. We want to read a triangulation from disk, and the class
// which does this is declared in the following file:
#include <deal.II/grid/grid_in.h>

// We will use a circular domain, and the object describing the boundary of it
// comes from this file:
#include <deal.II/grid/manifold_lib.h>

// This is C++ ...
#include <fstream>
#include <iostream>
#include <map>
#include <vector>


// Finally, this has been discussed in previous tutorial programs before:
using namespace dealii;

template <int dim>
class Step5
{
public:
  Step5();
  void run();

private:
  void output_results() const;


};


template <int dim>
Step5<dim>::Step5()
{}


template <int dim>
void Step5<dim>::output_results() const
{
//  DataOut<dim> data_out;
//
//  data_out.attach_dof_handler(dof_handler);
//  data_out.add_data_vector(solution, "solution");
//
//  data_out.build_patches();
//
//  std::ofstream output("solution-" + std::to_string(cycle) + ".vtu");
//  data_out.write_vtu(output);
}




template <int dim>
void Step5<dim>::run()
{

	const int spacedim = 3;
	// normal triangulation
	Triangulation<dim, spacedim> triangulation;

	// create map<cell_id, vector with all the points>
	std::map<unsigned int, std::vector<Point<spacedim>>> map_in;

	GridIn<dim, spacedim> grid_in;
    grid_in.attach_triangulation(triangulation);
    std::ifstream input_file("quad9.inp");

    Assert(dim == 2, ExcInternalError());

  	grid_in.read_ucd(input_file, map_in);

  	// sanity check to see if all points are in map_in
  	std::vector<Point<spacedim>> test = map_in[0];

  	// loop
  	for(auto &point : test)
  	{
  		std::cout << point[0] << ", " << point[1] << ", " << point[2] << std::endl;
  	}

    // Create finite element
    FE_Q<dim,spacedim> fe(1);

    // Create quadrature rule
    QGauss<dim> quad(2);
    QGaussLobatto<dim> quad_GL(3);

    // Create mapping
    MappingQ<dim, spacedim> mapping(1);

    MappingQGeneric<dim, spacedim> mapping_generic(2);

    // Create FEValues (for a single set of FiniteElement, Quadrature, Mapping)
    FEValues<dim, spacedim> fe_values(mapping_generic,
                            fe,
                            quad_GL,
                            update_quadrature_points | update_JxW_values);


    // loop over cells in triangulation
    for (auto &cell : triangulation.cell_iterators())
    {
    	fe_values.reinit(cell);
    	for (auto &q_point : fe_values.get_quadrature_points())
        	std::cout << "qp_x: " << q_point[0] << " qp_y: " << q_point[1] << " qp_z: " << q_point[2] << std::endl;

    }







}


// @sect3{The <code>main</code> function}

// The main function looks mostly like the one in the previous example, so we
// won't comment on it further:
int main()
{
  Step5<2> laplace_problem_2d;
  laplace_problem_2d.run();
  return 0;
}
