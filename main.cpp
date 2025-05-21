#include <iostream>
#include "Polyhedron.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"
#include <unordered_set>

using namespace std;
using namespace PolyhedronLibrary;

int main(int argc, char* argv[])
{   
    if(argc !=5 && argc!=7)
    {
        cerr<<"usage: "<<argv[0]<<" <p> <q> <b> <c> [<id_source> <id_destination>]"<<endl;
        return 1;
    }


    const int p = stoi(argv[1]);
    const int q = stoi(argv[2]);
    const int b = stoi(argv[3]);
    const int c = stoi(argv[4]);
    int id_source = -1;
    int id_destination = -1;

    Polyhedron P;

    //To build the platonic solid
    if(!Build_platonic_solid(p, q, P))
    {
        cerr<<"Error: the platonic solid {p,q} could not be imported, check the values of p and q"<<endl;
        return 2;
    }
    
    //To triangulate the polyhedron
    if((b < 1 && c == 0) || (b == 0 && c < 1) || (c <= 0 && b <= 0) || (b < 0 && c > 0) || (b > 0 && c < 0))
    {
        cerr<<"Error: the polyhedron could not be triangulated, check the values of b and c"<<endl;
        return 3;
    }
    else Triangulate(P, b, c);

    //To dualize the polyhedron
    if(p != 3 && q == 3) Dualize(P);

    //To project the polyhedron in the unitary sphere
    project_points_onto_sphere(P);

    //To visulize by terminal a polyhedron
    Visualize_polyhedron(P);

    //To create the CellXs.txt files
    Export_polyhedron(P); 
    
    

    //ShortPath property
    vector<Gedim::UCDProperty<double>> points_properties;
    vector<Gedim::UCDProperty<double>> segments_properties;
    points_properties.reserve(1); //We have only one Property
    segments_properties.reserve(1);

    vector<double> prop_vert(P.num_cell0Ds, 0.0);
    vector<double> prop_edges(P.num_cell1Ds, 0.0);

    //Fill the struct points_properties
    Gedim::UCDProperty<double> pointP;
    pointP.Label = "ShortPath";
    pointP.NumComponents = 1;
    pointP.Data = prop_vert.data();
    points_properties.push_back(pointP);
    
    //Fill the struct segments_properties
    Gedim::UCDProperty<double> edgeP;
    edgeP.Label = "ShortPath";
    edgeP.NumComponents = 1;
    edgeP.Data = prop_edges.data();
    segments_properties.push_back(edgeP);

    //Check if the id source and destination are valid
    if(argc == 7){
        id_source = stoi(argv[5]);
        id_destination = stoi(argv[6]);
        if((id_source < 0 || id_source >= P.num_cell0Ds) || (id_destination < 0 || id_destination >= P.num_cell0Ds))
        {
            cerr<<"Error: the shortes path between "<<id_source<<" and "<<id_destination<<" could not be imported, check the values of id_source and id_destination"<<endl;
            return 4;
        }
        else
        {
            vector<unsigned int> path = Short_path(P, id_source, id_destination);
        
            //Fill data prop_vert
            for(unsigned int id_v : path) prop_vert[id_v] = 1.0;

            //Fill data prop_edges
            for(unsigned int i=0; i < path.size() - 1; i++){
                unsigned int u = path[i];
                unsigned int v = path[i+1];
                std::unordered_set<unsigned int> edge_ref = {u, v};

                for(unsigned int j=0; j < P.cell1Ds_extrema.cols(); j++){
                    Eigen::Vector2i edge = P.cell1Ds_extrema.col(j);
                    unsigned int a = edge[0]; 
                    unsigned int b = edge[1];
                    if(edge_ref.count(a) + edge_ref.count(b) == 2) prop_edges[j] = 1.0;
                }
            }
        }
    }



    //To export the polyhedron in Paraview 
    Gedim::UCDUtilities utilities;
    utilities.ExportPoints("./Cell0Ds.inp",
                        P.cell0Ds_coordinates,
                        points_properties);

    utilities.ExportSegments("./Cell1Ds.inp",
                            P.cell0Ds_coordinates,
                            P.cell1Ds_extrema,
                            points_properties,
                        segments_properties);

    cout<<"Cell0Ds.inp, Cell1Ds.inp files created"<<endl;

    return 0;
}