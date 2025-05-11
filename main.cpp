#include <iostream>
#include "Polyhedron.hpp"
#include "Utils.hpp"
#include "UCDutilities.hpp"

using namespace std;
using namespace PolyhedronLibrary;

int main(int argc, char* argv[])
{   
    if(argc < 5)
    {
        cerr<<"usage: "<<argv[0]<<" <p> <q> <b> <c> [<begin id vertex> <end id vertex>]"<<endl;
        return 1;
    }
    
    const unsigned int p = stoi(argv[1]);
    const unsigned int q = stoi(argv[2]);
    const unsigned int b = stoi(argv[3]);
    const unsigned int c = stoi(argv[4]);

    Polyhedron P;
    if(!Import_platonic_solid(p, q, P))
    {
        cerr<<"Error: the platonic solid {p,q} could not be imported, check the value of p and q"<<endl;
        return 2;
    }
    
    //To visulize by terminal a polyhedron 
    Visualize_polyhedron(P);

    //To create the CellXs.txt files
    Export_polyhedron(P); 
    
     //To export the polyhedron in Paraview
    Gedim::UCDUtilities utilities;
    utilities.ExportPoints("./Cell0Ds.inp",
                           P.cell0Ds_coordinates);

    utilities.ExportSegments("./Cell1Ds.inp",
                             P.cell0Ds_coordinates,
                             P.cell1Ds_extrema);

    

    return 0;
}