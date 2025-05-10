#pragma once

#include <iostream>
#include "Polyhedron.hpp"

using namespace std;

namespace PolyhedronLibrary{

/*Function to import a platonic solid*/
bool Import_platonic_solid(const unsigned int &p, const unsigned int &q, Polyhedron &P);

/*Function to export the polyhedron properties and it will create 4 files: Cell0Ds.txt, Cell1Ds.txt, Cell2Ds.txt, Cell3Ds.txt*/
bool Export_polyhedron(Polyhedron &P);

}



