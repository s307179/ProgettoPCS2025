#pragma once

#include <iostream>
#include "Polyhedron.hpp"

using namespace std;

namespace PolyhedronLibrary{

/*Function to import a platonic solid
p,q: unsigned int for the Schlafli's Simbol {p, q}
P: a Polyhedron struct
return: the resul of import, true if is success, false otherwise */
bool Import_platonic_solid(unsigned int p, unsigned int q, Polyhedron &P);


/*Function to finish to fill the polyhedron struct using cell0Ds_coordinates and cell2Ds_vertices
P: a Polyhedron struct 
BE CAREFUL, should be used only after fill in the struct the numbers of cell, the coordinates matrix and cell2Ds_vertices*/
void finish_to_fill_struct(Polyhedron &P);


/*Function to export the polyhedron properties and it will create 4 files: Cell0Ds.txt, Cell1Ds.txt, Cell2Ds.txt, Cell3Ds.txt
P: a Polyhedron struct */
void Export_polyhedron(Polyhedron &P);


/*Function to visulaize every properties of a input Polyhedron struct
P: a Polyhedron struct*/
void Visualize_polyhedron(Polyhedron &P);


/*Function that compute the triangulation of ClassI by input triangle ABC
A: Eigen::Vector3d is one of the vertex of the input triangle
B: Eigen::Vector3d is one of the vertex of the input triangle
C: Eigen::Vector3d is one of the vertex of the input triangle
b: an unsigned int that lead the triangulation
return: a pair where the first element is a vector that contains the generated vertices, 
        while the second one is a vector that contains the generated faces (with the local id vertices that make it)*/
pair<vector<Eigen::Vector3d>, vector<Eigen::Vector3i>> classI_basic_step(const Eigen::Vector3d &A, const Eigen::Vector3d &B, const Eigen::Vector3d &C, const unsigned int b);


/*Function that compute the triangulation of ClassI by input polyhedron using the function "Triangulation_basic_step"
P: a Polyhedron struct
b: an unsigned int that lead the triangulation
p,q: unsigned int for the Schlafli's Simbol {p, q} */
void Triangulate(Polyhedron &P, const unsigned int b,const unsigned int c);


/* Function to project each point of the polyhedron onto the sphere in the origin with unitary radius
P: a Polyhedron struct */
void project_points_onto_sphere(Polyhedron &P);


/*Function to dualize an input polyhedron
P: a Polyhedron struct*/
void Dualize(Polyhedron &P);


/*Function to make cyclic an input face
face_new: id of the face to order
coord: the matrix with the coordinates of every point
return: a cyclic version of face_new, with ordered ids*/
vector<unsigned int> cycled_face_for_dual(vector<unsigned int> &face_new, Eigen::MatrixXd &coord);


/*Function that compute the triangulation of ClassII by input triangle ABC
A: Eigen::Vector3d is one of the vertex of the input triangle
B: Eigen::Vector3d is one of the vertex of the input triangle
C: Eigen::Vector3d is one of the vertex of the input triangle
b: an unsigned int that lead the triangulation
return: a pair where the first element is a vector that contains the generated vertices, 
        while the second one is a vector that contains the generated faces (with the local id vertices that make it)*/
pair<vector<Eigen::Vector3d>, vector<Eigen::Vector3i>> classII_basic_step(const Eigen::Vector3d &A, const Eigen::Vector3d &B, const Eigen::Vector3d &C, const unsigned int b);

/*Function that computes the shortest path between the given vertices id_D and id_A
and print the properties of the path (e.g length)
P: a Polyhedron struct
id_D, id_A: id of two vertices of the polyhedron P*/
vector<unsigned int> Short_path(Polyhedron &P, unsigned int id_D, unsigned int id_A);


}



