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
    if(! Import_platonic_solid(p, q, P))
    {
        cerr<<"Error: the platonic solid {p,q} could not be imported"<<endl;
        return 2;
    }
    
    Eigen::MatrixXd &A = P.cell0Ds_coordinates;
    Eigen::MatrixXi &B = P.cell1Ds_extrema;

    Eigen::MatrixXi &C = P.cell2Ds_vertices;
    Eigen::MatrixXi &D = P.cell2Ds_edges;

    cout<<A<<endl;
    cout<<endl;

    cout<<B<<endl;
    cout<<endl;
    
    cout<<C<<endl;
    cout<<endl;

    cout<<D<<endl;
    cout<<endl;

    

    return 0;
}