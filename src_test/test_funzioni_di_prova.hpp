#pragma once

#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include "funzioni_di_prova.cpp"
#include "Polyhedron"
namespace PolyhedronLibrary {
	// cubo già normalizzato
TEST(TestPolyhedron, TestProjection)
{	
    Polyhedron solid;
	int n=9;
	solid.numcell0ds = n;
	solid.cell0Ds_coordinates << {0.0, 0.0, 1.224745,
	1.154701, 0.0, 0.4082483,
	-0.5773503, 1.0, 0.4082483,
	-0.5773503, -1.0, 0.4082483,
	0.5773503, 1.0, -0.4082483,
	0.5773503, -1.0, -0.4082483,
	-1.154701, 0.0, -0.4082483,		
	0.0, 0.0, -1.224745}
    project_points_onto_sphere(solid);
    Polyhedron projectedSolid;
	projectedSolid.numcell0ds = n;
	projectedSolid.cell0Ds_coordinates << {0.0, 0.0, 1.224745,
	1.154701, 0.0, 0.4082483,
	-0.5773503, 1.0, 0.4082483,
	-0.5773503, -1.0, 0.4082483,
	0.5773503, 1.0, -0.4082483,
	0.5773503, -1.0, -0.4082483,
	-1.154701, 0.0, -0.4082483,		
	0.0, 0.0, -1.224745};
	EXPECT_EQ(solid, projectedSolid);
}
// tetraedro da normalizzare
TEST(TestPolyhedron, TestProjection)
{
	Polyhedron solid;
	int n=4;
	solid.numcell0ds = n;
	solid.cell0Ds_coordinates << {2.0, 2.0, 2.0,
            -2.0, -2.0, 2.0,
            -2.0, 2.0, -2.0,
             2.0, -2.0, -2.0};
    project_points_onto_sphere(solid);
    Polyhedron projectedSolid;
	projectedSolid.numcell0ds = n;
	projectedSolid.cell0Ds_coordinates<< {0.577, 0.577, 0.577,
	-0.577, -0.577, 0.577,
	-0.577, 0.577, -0.577,
	0.577, -0.577, -0.577};
	EXPECT_EQ(solid, projectedSolid);
}