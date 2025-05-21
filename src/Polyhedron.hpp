#pragma once

#include <iostream>
#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolyhedronLibrary{

struct Polyhedron
{
    unsigned int num_cell0Ds;
    unsigned int num_cell1Ds;
    unsigned int num_cell2Ds;

    vector<unsigned int> cell0Ds_id;
    vector<unsigned int> cell1Ds_id;
    vector<unsigned int> cell2Ds_id;

    MatrixXd cell0Ds_coordinates;
    MatrixXi cell1Ds_extrema;

    vector<vector<unsigned int>> cell2Ds_vertices;
    vector<vector<unsigned int>> cell2Ds_edges;
};

}