#pragma once

#include <iostream>
#include "Polyhedron.hpp"

using namespace std;

namespace PolyhedronLibrary{

/*Function to import a platonic solid*/
bool Import_platonic_solid(unsigned int &p, unsigned int &q, Polyhedron &P);

bool Export_polyhedron(Polyhedron &P);

}



