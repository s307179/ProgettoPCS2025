#include "Utils.hpp"

#include <iostream>
#include <fstream>


namespace PolyhedronLibrary{

bool Import_platonic_solid(const unsigned int &p, const unsigned int &q, Polyhedron &P)
{
    //To establish the corresponding platonic solid {p,q}
    string poly_name = "dodecahedron";
    unsigned int V = 20;
    unsigned int E = 30;
    unsigned int F = 12;
    if(p == 3 && q == 3) 
    {
        poly_name = "tetrahedron";
        V = 4;
        E = 6;
        F = 4;   
    }
    else if (p == 3 && q == 4) 
    {
        poly_name = "octahedron";
        V = 6;
        E = 12;
        F = 8;    
    }
    else if (p == 3 && q == 5) 
    {
        poly_name = "icosahedron";
        V = 12;
        E = 30;
        F = 20;    
    }
    else if (p == 4 && q == 3) 
    {
        poly_name = "cube";
        V = 8;
        E = 12;
        F = 6;    
    }
    
    P.num_cell0Ds = V;
    P.num_cell1Ds = E;
    P.num_cell2Ds = F;
    P.num_cell3Ds = 1;
    
    P.cell0Ds_id.reserve(V);
    P.cell1Ds_id.reserve(E);
    P.cell2Ds_id.reserve(F);
    P.cell3D_id = 0;

    ifstream ifile("./" + poly_name + ".txt");
    if(ifile.fail()) return false;

    //Cell0Ds properties (vertices)
    //Fill the matrix cell0Ds_coordinates <--> matrix A
    MatrixXd &A = P.cell0Ds_coordinates; 
    A = MatrixXd::Zero(3, V);
    for(size_t i=0; i < V; i++)
    {
        P.cell0Ds_id.push_back(i);
        
        char trash;
        double x, y, z;
        ifile >> trash >> x >> y >> z;
        A(0,i) = x;
        A(1,i) = y;
        A(2,i) = z;    
    }

    //Cell2Ds properties (faces)
    //Fill the matrix cell2Ds_vertices <--> matrix B 
    MatrixXi &B = P.cell2Ds_vertices;
    MatrixXi &C = P.cell2Ds_edges;
    B = MatrixXi::Zero(p, F);
    C = MatrixXi::Zero(p, F);
    for(size_t i=0; i < F; i++)
    {
        P.cell2Ds_id.push_back(i);

        char trash;
        ifile >> trash;
        for(size_t j=0; j < p; j++)
        {
            unsigned int u;
            ifile >> u;
            B(j,i) = u;
        }
    }

    
    //Fill the matrix cell1Ds_extrema <--> matrix D
    MatrixXi &D = P.cell1Ds_extrema;
    D = MatrixXi::Zero(2, E);

    VectorXi check_valence;
    check_valence = VectorXi::Zero(P.num_cell0Ds);
    unsigned int id_edge = 0;
    for(unsigned int i=0; i < B.cols(); i++)
    {   
        VectorXi face = B.col(i);
        for(unsigned int k=0; k < p; k++)
        {
            unsigned int origin = face[k];
            unsigned int end = face[(k+1)%p];

            if(check_valence[origin] < q && check_valence[end] < q)
            {   
                bool not_added = true;

                for(unsigned int j=0; j < id_edge; j++)
                {
                    VectorXi edge = D.col(j);
                    unsigned int old_origin = edge[0];
                    unsigned int old_end = edge[1];

                    if(old_origin == end && old_end == origin)
                    {
                        not_added = false;
                    }
                }

                if(not_added)
                {
                    P.cell1Ds_id.push_back(id_edge);

                    D(0,id_edge) = origin;
                    D(1,id_edge) = end;
                    id_edge++;

                    check_valence[origin]++;
                    check_valence[end]++;
                }
                
            }
        }
    }

    //Fill the matrix cell2Ds_edges <--> matrix C
    for(unsigned int i=0; i < B.cols(); i++)
    {
        VectorXi face = B.col(i); //matrix B --> cell2Ds_vertices
        for(unsigned int j=0; j < p; j++)
        {
            unsigned int u = face[j];
            unsigned int v = face[(j+1)%p];

            for(unsigned int k=0; k < D.cols(); k++) //matrix D --> cell1Ds_extrema
            {
                VectorXi edge = D.col(k);
                unsigned int origin = edge[0];
                unsigned int end = edge[1];

                if ((u == origin && v == end) || (u == end && v == origin)) 
                {
                    C(j,i) = k;
                }    
            }
        }
    }
    
    return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Export_polyhedron(Polyhedron &P)
{
    return true;
}

}