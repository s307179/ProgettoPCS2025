#pragma once

#include <iostream>
#include "Polyhedron.hpp"

using namespace std;

namespace PolyhedronLibrary{

/*Function to import a platonic solid*/
bool Import_solid(string &solid_name, Polyhedron &P);

/*Import the Cell0D properties from Cell0Ds.txt file*/
bool import_cell0Ds(string &solid_name, Polyhedron &P);

/*Import the Cell1D properties from Cell1Ds.txt file*/
bool import_cell1Ds(string &solid_name, Polyhedron &P);

/*Import the Cell2D properties from Cell2Ds.txt file*/
bool import_cell2Ds(string &solid_name, Polyhedron &P);

}



