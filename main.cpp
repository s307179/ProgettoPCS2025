#include <iostream>
#include "Polyhedron.hpp"
#include "Utils.hpp"
#include "UCDutilities.hpp"

using namespace std;
using namespace PolyhedronLibrary;

int main(int argc, char* argv[])
{   
    if(argc !=5 && argc!=7)
    {
        cerr<<"usage: "<<argv[0]<<" <p> <q> <b> <c> [<begin id vertex> <end id vertex>]"<<endl;
        return 1;
    }
    
    const unsigned int p = stoi(argv[1]);
    const unsigned int q = stoi(argv[2]);
    const unsigned int b = stoi(argv[3]);
    const unsigned int c = stoi(argv[4]);

    Polyhedron P;

    //To build the platonic solid
    if(!Import_platonic_solid(p, q, P))
    {
        cerr<<"Error: the platonic solid {p,q} could not be imported, check the values of p and q"<<endl;
        return 2;
    }

    
    //To triangulate the polyhedron
    if((b < 1 && c < 1))
    {
        cerr<<"Error: the polyhedron could not be triangulated, check the values of b and c"<<endl;
        return 3;
    }

    if((b >= 1 && c == 0) || (b == 0 && c >=1) || (b == c && b != 0)) Triangulate(P, b, c);
    

    //To dualize the polyhedron
    if(p != 3 && q == 3) Dualize(P);

    //To project the polyhedron in the unitary sphere
    project_points_onto_sphere(P);

    //To visulize by terminal a polyhedron
    Visualize_polyhedron(P);

    //To create the CellXs.txt files
    Export_polyhedron(P); 
    
	//TODO aggiungi check sugli id e gestisci gli id fa command line
    Short_path(P, 0, 7);

    //To export the polyhedron in Paraview 
    Gedim::UCDUtilities utilities;
    utilities.ExportPoints("./Cell0Ds.inp",
                           P.cell0Ds_coordinates);

    utilities.ExportSegments("./Cell1Ds.inp",
                             P.cell0Ds_coordinates,
                             P.cell1Ds_extrema);
    
	
    

    return 0;
}