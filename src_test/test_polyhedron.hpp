#pragma once

#include <iostream>
#include <vector>
#include <gtest/gtest.h>

#include "Utils.hpp"
#include "Eigen/Eigen"
#include "Polyhedron.hpp"
#include "Utils.hpp"


using namespace Eigen;
using namespace PolyhedronLibrary;

namespace PolyhedronLibrary {
	// test the correct projection of an already normalized solid
TEST(TestPolyhedron, TestProjection1)
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


	// test the correct projection of a not normalized solid
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
}


	// Test the dualization of the tetrahedron
TEST(TestPolyhedron, TestDualization1)
{	
	Polyhedron P;
	unsigned int p = 3;
	unsigned int q = 3;
	ASSERT_TRUE(Import_platonic_solid(p, q, P));
	Dualize(P);
	unsigned int n = P.num_cell0Ds;
	unsigned int m = P.num_cell1Ds;
	MatrixXd expected_coords(n,3);
	expected_coords <<  0.272166 , 0.471405 , 0.192450,
 -0.544331,  0.000000 , 0.192450,
0.272166 ,-0.471405  ,0.192450,
0.000000  ,0.000000 ,-0.577350;
	MatrixXi expected_edges(m,2);
	expected_edges << 0,1,
			1,2,
			0,2,
			2,3,
			0,3,
			1,3;
	bool coords_match = P.cell0Ds_coordinates.isApprox(expected_coords.transpose(),1e-6);
	bool edges_match = P.cell1Ds_extrema.isApprox(expected_edges.transpose());
	EXPECT_TRUE(coords_match);
	EXPECT_TRUE(edges_match);
}


	// Test the dualization of an octahedron
TEST(TestPolyhedron, TestDualization2)
{	
	Polyhedron P;
	unsigned int p = 3;
	unsigned int q = 4;
	ASSERT_TRUE(Import_platonic_solid(p, q, P));
	Dualize(P);
	unsigned int n = P.num_cell0Ds;
	unsigned int m = P.num_cell1Ds;
	MatrixXd expected_coords(n,3);
	expected_coords <<     0.471405,  0.471405,  0.471405,
-0.471405,  0.471405,  0.471405,
 -0.471405, -0.471405,  0.471405,
  0.471405, -0.471405,  0.471405,
  0.471405, -0.471405, -0.471405,
    0.471405,  0.471405, -0.471405,
 -0.471405,  0.471405, -0.471405,
 -0.471405, -0.471405, -0.471405;
	MatrixXi expected_edges(m,2);
	expected_edges << 0,1,
					  1,2,
					  2,3,
					  0,3,
					  3,4,
					  4,5,
					  0,5,
					  1,6,
					  5,6,
					  2,7,
					  6,7,
					  4,7;		
	bool coords_match = P.cell0Ds_coordinates.isApprox(expected_coords.transpose(),1e-6);
	bool edges_match = P.cell1Ds_extrema.isApprox(expected_edges.transpose());
	EXPECT_TRUE(coords_match);
	EXPECT_TRUE(edges_match);
};


}























