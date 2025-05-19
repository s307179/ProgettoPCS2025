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
	mat <<  0.57735,  0.57735,  0.57735,
-0.57735,  0.57735,  0.57735,
 0.57735, -0.57735,  0.57735,
-0.57735, -0.57735,  0.57735,
 0.57735,  0.57735, -0.57735,
-0.57735,  0.57735, -0.57735,
 0.57735, -0.57735, -0.57735,
-0.57735, -0.57735, -0.57735;
	solid.cell0Ds_coordinates = mat.transpose();
	Polyhedron projectedSolid;
	projectedSolid.num_cell0Ds = n;
	projectedSolid.cell0Ds_coordinates = solid.cell0Ds_coordinates;
    project_points_onto_sphere(solid);
	bool coords_match = solid.cell0Ds_coordinates.isApprox(projectedSolid.cell0Ds_coordinates,1e-5);
	EXPECT_TRUE(coords_match);
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
	bool coords_match = solid.cell0Ds_coordinates.isApprox(projectedSolid.cell0Ds_coordinates,1e-3);
	EXPECT_TRUE(coords_match);
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
 -0.471405, -0.471405,- 0.471405;
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

// Test the classI of a tetrahedron (b=1,c=0)
TEST(TestPolyhedron, TestClassI0)
{	
	Polyhedron P;
	unsigned int p = 3;
	unsigned int q = 3;
	unsigned int b = 1;
	unsigned int c = 0;
	ASSERT_TRUE(Import_platonic_solid(p, q, P));
	Triangulate(P,b,c);
	unsigned int n = P.num_cell0Ds;
	unsigned int m = P.num_cell1Ds;
	MatrixXd expected_coords(n,3);
	expected_coords <<   -0.816497,1.41421,-0.57735,
							1.63299,0,-0.57735,
						0,0,1.73205,
						-0.816497,-1.41421,-0.57735;
	MatrixXi expected_edges(m,2);
	expected_edges << 	0,2,
						1,2,
						0,1,
						2,3,
						0,3,
						1,3;
	bool coords_match = P.cell0Ds_coordinates.isApprox(expected_coords.transpose(),1e-5);
	bool edges_match = P.cell1Ds_extrema.isApprox(expected_edges.transpose());
	EXPECT_TRUE(coords_match);
	EXPECT_TRUE(edges_match);
};

// Test the classI of a tetrahedron (b=0,c=1)
TEST(TestPolyhedron, TestClassI1)
{	
	Polyhedron P;
	unsigned int p = 3;
	unsigned int q = 3;
	unsigned int b = 0;
	unsigned int c = 1;
	ASSERT_TRUE(Import_platonic_solid(p, q, P));
	Triangulate(P,b,c);
	unsigned int n = P.num_cell0Ds;
	unsigned int m = P.num_cell1Ds;
	MatrixXd expected_coords(n,3);
	expected_coords <<   -0.816497,1.41421,-0.57735,
							1.63299,0,-0.57735,
						0,0,1.73205,
						-0.816497,-1.41421,-0.57735;
	MatrixXi expected_edges(m,2);
	expected_edges << 	0,2,
						1,2,
						0,1,
						2,3,
						0,3,
						1,3;
	bool coords_match = P.cell0Ds_coordinates.isApprox(expected_coords.transpose(),1e-5);
	bool edges_match = P.cell1Ds_extrema.isApprox(expected_edges.transpose());
	EXPECT_TRUE(coords_match);
	EXPECT_TRUE(edges_match);
};


// Test the ClassII of a tetrahedron (b=c=1)

TEST(TestPolyhedron, TestClassII)
{	
	Polyhedron P;
	unsigned int p = 3;
	unsigned int q = 3;
	unsigned int b = 1;
	unsigned int c = 1;
	ASSERT_TRUE(Import_platonic_solid(p, q, P));
	Triangulate(P,b,c);
	MatrixXd expected_coords(14,3);
	expected_coords <<-0.816497,1.41421,-0.57735,
						0,0,1.73205,
					1.63299,0,-0.57735,
					-0.408248,0.707107,0.57735,
				0.816496,0,0.57735,
			0.408248,0.707107,-0.57735,
			0.272165,0.471405,0.19245,
		-0.816497,-1.41421,-0.57735,
		-0.408248,-0.707107,0.57735,
		-0.816497,0,-0.57735,
		-0.544331,0,0.19245,
		0.408248,-0.707107,-0.57735,
		0.272165,-0.471405,0.19245,
		-6.66667e-08,0,-0.57735;
	MatrixXi expected_edges(36,2);
	expected_edges << 0,6,3,6,0,3,5,6,0,5,2,6,2,5,4,6,2,4,1,6,1,4,1,3,7,10,8,10,7,
	8,9,10,7,9,0,10,0,9,3,10,1,10,1,8,2,12,4,12,11,12,2,11,7,12,7,11,
	8,12,1,12,0,13,5,13,9,13,7,13,11,13,2,13;			  
	bool coords_match = P.cell0Ds_coordinates.isApprox(expected_coords.transpose(),1e-5);
	bool edges_match = P.cell1Ds_extrema.isApprox(expected_edges.transpose());
	EXPECT_TRUE(coords_match);
	EXPECT_TRUE(edges_match);
};



















}























