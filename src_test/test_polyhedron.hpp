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
TEST(TestPolyhedron, TestClassI1)
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
TEST(TestPolyhedron, TestClassI2)
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

TEST(TestPolyhedron, TestClassII1)
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

// Test the ClassII of a tetrahedron (b=c=2)

TEST(TestPolyhedron, TestClassII2)
{	
	Polyhedron P;
	unsigned int p = 3;
	unsigned int q = 3;
	unsigned int b = 2;
	unsigned int c = 2;
	ASSERT_TRUE(Import_platonic_solid(p, q, P));
	Triangulate(P,b,c);
	MatrixXd expected_coords(38,3);
	expected_coords <<-8.1649660000000002e-01,1.4142140000000001e+00,-5.7735029999999998e-01,
-4.0824830000000001e-01,7.0710700000000004e-01,5.7735035000000001e-01,
4.0824819999999995e-01,7.0710700000000004e-01,-5.7735029999999998e-01,
-6.1237245000000007e-01,1.0606605000000000e+00,2.5000000014596679e-08,
-2.0412420000000003e-01,1.0606605000000000e+00,-5.7735029999999998e-01,
-2.7216556666666675e-01,9.4280933333333339e-01,-1.9245008333333333e-01,
8.1649649999999996e-01,0.0000000000000000e+00,5.7735035000000001e-01,
2.7216546666666663e-01,4.7140466666666669e-01,1.9245013333333336e-01,
1.6329929999999999e+00,0.0000000000000000e+00,-5.7735029999999998e-01,
1.2247447499999999e+00,0.0000000000000000e+00,2.5000000014596679e-08,
1.0206206000000000e+00,3.5355350000000002e-01,-5.7735029999999998e-01,
9.5257923333333328e-01,2.3570233333333335e-01,-1.9245008333333333e-01,
0.0000000000000000e+00,0.0000000000000000e+00,1.7320510000000000e+00,
-2.0412415000000000e-01,3.5355350000000002e-01,1.1547006750000000e+00,
4.0824824999999998e-01,0.0000000000000000e+00,1.1547006750000000e+00,
1.3608273333333332e-01,2.3570233333333335e-01,9.6225056666666664e-01,
-8.1649660000000002e-01,-1.4142140000000001e+00,-5.7735029999999998e-01,
-4.0824830000000001e-01,-7.0710700000000004e-01,5.7735035000000001e-01,
-8.1649660000000002e-01,0.0000000000000000e+00,-5.7735029999999998e-01,
-6.1237245000000007e-01,-1.0606605000000000e+00,2.5000000014596679e-08,
-8.1649660000000002e-01,-7.0710700000000004e-01,-5.7735029999999998e-01,
-6.8041383333333327e-01,-7.0710700000000004e-01,-1.9245008333333333e-01,
-5.4433106666666664e-01,0.0000000000000000e+00,1.9245013333333336e-01,
-8.1649660000000002e-01,7.0710700000000004e-01,-5.7735029999999998e-01,
-6.8041383333333327e-01,7.0710700000000004e-01,-1.9245008333333333e-01,
-2.0412415000000000e-01,-3.5355350000000002e-01,1.1547006750000000e+00,
-2.7216553333333332e-01,0.0000000000000000e+00,9.6225056666666664e-01,
4.0824819999999995e-01,-7.0710700000000004e-01,-5.7735029999999998e-01,
1.0206206000000000e+00,-3.5355350000000002e-01,-5.7735029999999998e-01,
9.5257923333333328e-01,-2.3570233333333335e-01,-1.9245008333333333e-01,
2.7216546666666663e-01,-4.7140466666666669e-01,1.9245013333333336e-01,
-2.0412420000000003e-01,-1.0606605000000000e+00,-5.7735029999999998e-01,
-2.7216556666666669e-01,-9.4280933333333339e-01,-1.9245008333333333e-01,
1.3608273333333332e-01,-2.3570233333333335e-01,9.6225056666666664e-01,
-4.0824833333333332e-01,7.0710700000000004e-01,-5.7735029999999998e-01,
-6.6666666705591141e-08,0.0000000000000000e+00,-5.7735029999999998e-01,
-4.0824833333333332e-01,-7.0710700000000004e-01,-5.7735029999999998e-01,
8.1649646666666664e-01,0.0000000000000000e+00,-5.7735029999999998e-01;

	MatrixXi expected_edges(108,2);
	expected_edges << 0,5,3,5,0,3,4,5,0,4,2,5,2,4,1,5,1,3,2,11,10,11,2,10,8,11,8,10,9,11,8,9,6,11,6,9,1,15,13,15,1,
	13,6,15,14,15,6,14,12,15,12,14,12,13,6,7,7,11,2,7,7,15,1,7,5,7,16,21,19,21,16,19,20,21,16,20,18,21,18,20,17,21,
	17,19,18,24,23,24,18,23,0,24,0,23,3,24,1,24,17,26,25,26,17,25,1,26,13,26,12,26,12,25,1,22,22,24,18,22,22,26,17,
	22,21,22,8,29,9,29,28,29,8,28,27,29,27,28,6,29,27,32,31,32,27,31,16,32,16,31,19,32,17,32,6,33,14,33,17,33,25,33,
	12,33,17,30,30,32,27,30,30,33,6,30,29,30,0,34,4,34,23,34,18,34,2,34,18,36,20,36,16,36,31,36,27,36,2,37,10,37,27,
	37,28,37,8,37,27,35,35,36,18,35,35,37,2,35,34,35;  
	bool coords_match = P.cell0Ds_coordinates.isApprox(expected_coords.transpose(),1e-4);
	bool edges_match = P.cell1Ds_extrema.isApprox(expected_edges.transpose());
	EXPECT_TRUE(coords_match);
	EXPECT_TRUE(edges_match);
};






// Test the shortest path algorithm 1

 TEST(TestPolyhedron, TestShortestPath1){
	Polyhedron P;
	unsigned int p = 3;
	unsigned int q = 4;
	bool t = Import_platonic_solid(p, q, P);
	if(!t){
		cerr <<"Could not import the platonic solid"<<endl;
	}
	Dualize(P);
	project_points_onto_sphere(P);
	unsigned int id_D = 0;
	unsigned int id_A = 7;
	vector<unsigned int> path = Short_path(P, id_D, id_A);

	double length = 0.0;
    for(unsigned int i=0; i < path.size() - 1; i++){
        unsigned int id_U = path[i];
        unsigned int id_V = path[i+1];
        Eigen::Vector3d U = P.cell0Ds_coordinates.col(id_U);
        Eigen::Vector3d V = P.cell0Ds_coordinates.col(id_V);
        length += (U-V).norm();

    }
	EXPECT_DOUBLE_EQ(length, 6/(sqrt(3)));
 };
	//Test 2 of the shortest path

TEST(TestPolyhedron, TestShortestPath2){
	Polyhedron Q;
	unsigned int p2 = 3;
	unsigned int q2 = 4;
	bool t2 = Import_platonic_solid(p2, q2, Q);
	if(!t2){
		cerr <<"Could not import the platonic solid"<<endl;
	}
	project_points_onto_sphere(Q);
	unsigned int id_D2 = 0;
	unsigned int id_A2 = 5;
	vector<unsigned int> path2 = Short_path(Q, id_D2, id_A2);
	

	double length2 = 0.0;
    for(unsigned int j=0; j < path2.size() - 1; j++){
        unsigned int id_U = path2[j];
        unsigned int id_V = path2[j+1];
        Eigen::Vector3d U = Q.cell0Ds_coordinates.col(id_U);
        Eigen::Vector3d V = Q.cell0Ds_coordinates.col(id_V);
        length2 += (U-V).norm();

    }
	EXPECT_DOUBLE_EQ(length2, 2*sqrt(2));
};





















}























