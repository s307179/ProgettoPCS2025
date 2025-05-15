#pragma once

#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include "Utils.hpp"
#include "Eigen/Eigen"
using namespace Eigen;
namespace PolyhedronLibrary {
	// cubo già normalizzato
TEST(TestPolyhedron, TestProjection)
{	
    Polyhedron solid;
	int n=8;
	solid.num_cell0Ds = n;
	MatrixXd mat(n,3);
	mat << 0.0, 0.0, 1.224745,
	1.154701, 0.0, 0.4082483,
	-0.5773503, 1.0, 0.4082483,
	-0.5773503, -1.0, 0.4082483,
	0.5773503, 1.0, -0.4082483,
	0.5773503, -1.0, -0.4082483,
	-1.154701, 0.0, -0.4082483,		
	0.0, 0.0, -1.224745;
	solid.cell0Ds_coordinates = mat.transpose();
    project_points_onto_sphere(solid);
    Polyhedron projectedSolid;
	projectedSolid.num_cell0Ds = n;
	MatrixXd mat2(n,3);
	mat2 << 0.0, 0.0, 1.224745,
	1.154701, 0.0, 0.4082483,
	-0.5773503, 1.0, 0.4082483,
	-0.5773503, -1.0, 0.4082483,
	0.5773503, 1.0, -0.4082483,
	0.5773503, -1.0, -0.4082483,
	-1.154701, 0.0, -0.4082483,		
	0.0, 0.0, -1.224745;
	projectedSolid.cell0Ds_coordinates = mat2.transpose();
	bool Test = true;
	if (!solid.cell0Ds_coordinates.isApprox(projectedSolid.cell0Ds_coordinates))
		Test = false;
	EXPECT_EQ(Test, 0);
}


// tetraedro da normalizzare
TEST(TestPolyhedron, TestProjection2)
{	
    Polyhedron solid;
	int n=4;
	solid.num_cell0Ds = n;
	MatrixXd mat(n,3);
	mat << 2.0, 2.0, 2.0,
            -2.0, -2.0, 2.0,
            -2.0, 2.0, -2.0,
             2.0, -2.0, -2.0;
	solid.cell0Ds_coordinates = mat.transpose();
    project_points_onto_sphere(solid);
    Polyhedron projectedSolid;
	projectedSolid.num_cell0Ds = n;
	MatrixXd mat2(n,3);
	mat2 << 0.577, 0.577, 0.577,
	-0.577, -0.577, 0.577,
	-0.577, 0.577, -0.577,
	0.577, -0.577, -0.577;
	projectedSolid.cell0Ds_coordinates = mat2.transpose();
	bool Test = true;
	if (!solid.cell0Ds_coordinates.isApprox(projectedSolid.cell0Ds_coordinates))
		Test = false;
	EXPECT_EQ(Test, 0);
}}
























