#pragma once

#include <iostream>
#include "Polyhedron.hpp"

using namespace std;

namespace PolyhedronLibrary{

/*Function to import a platonic solid
p,q: Schlafli's Simbol {p, q}
P: a Polyhedron struct
return the resul of import, true if is success, false otherwise */
bool Import_platonic_solid(unsigned int p, unsigned int q, Polyhedron &P);

/*Function to export the polyhedron properties and it will create 4 files: Cell0Ds.txt, Cell1Ds.txt, Cell2Ds.txt, Cell3Ds.txt
P: a Polyhedron struct */
void Export_polyhedron(Polyhedron &P);

/*Function to visulaize every properties of a input Polyhedron struct
P: a Polyhedron struct*/
void Visualize_polyhedron(Polyhedron &P);

/////////////////////////////////////////////////////////////
/*Function that compute the triangulation by input triangle ABC
*/
pair<vector<Eigen::Vector3d>, vector<Eigen::Vector3i>> Triangulation_basic_step(const Eigen::Vector3d &A, const Eigen::Vector3d &B, const Eigen::Vector3d &C, const unsigned int b);


void ClassI_polyhedron(Polyhedron &P, const unsigned int b, unsigned int p, unsigned int q);



}



