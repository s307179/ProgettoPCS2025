#pragma once

#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include "Utils.hpp"

namespace PolyhedronLibrary {
	// cubo già normalizzato
TEST(TestPolyhedron, TestProjection)
{	
    std::Matrix3d points; 
	points << {0.0, 0.0, 1.224745,
	1.154701, 0.0, 0.4082483,
	-0.5773503, 1.0, 0.4082483,
	-0.5773503, -1.0, 0.4082483,
	0.5773503, 1.0, -0.4082483,
	0.5773503, -1.0, -0.4082483,
	-1.154701, 0.0, -0.4082483,		
	0.0, 0.0, -1.224745}
    project_points_onto_sphere(points);
    std::Matrix3d projectedPoints;
	projectedPoints << {0.0, 0.0, 1.224745,
	1.154701, 0.0, 0.4082483,
	-0.5773503, 1.0, 0.4082483,
	-0.5773503, -1.0, 0.4082483,
	0.5773503, 1.0, -0.4082483,
	0.5773503, -1.0, -0.4082483,
	-1.154701, 0.0, -0.4082483,		
	0.0, 0.0, -1.224745};
	EXPECT_EQ(points, projectedPoints);
}
// tetraedro da normalizzare
TEST(TestPolyhedron, TestProjection)
{
    std::Matrix3d points; 
	points << {2.0, 2.0, 2.0,
            -2.0, -2.0, 2.0,
            -2.0, 2.0, -2.0,
             2.0, -2.0, -2.0};
    project_points_onto_sphere(points);
    std::Matrix3d projectedPoints;
	projectedPoints << {0.577; 0.577; 0.577;
	-0.577, -0.577, 0.577;
	-0.577, 0.577, -0.577;
	0.577, -0.577, -0.577}
	EXPECT_EQ(points, projectedPoints);
}


}