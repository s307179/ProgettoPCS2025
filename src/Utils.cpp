#include "Utils.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>

namespace PolyhedronLibrary{

bool Import_platonic_solid(const unsigned int &p, const unsigned int &q, Polyhedron &P)
{
    //To establish the corresponding platonic solid {p,q}
    string poly_name;
    unsigned int V, E, F;
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
    else if(p == 5 && q == 3)
    {
        poly_name = "dodecahedron";
        V = 20;
        E = 30;
        F = 12;
    }
    else return false;
    
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
void Export_polyhedron(Polyhedron &P)
{
    //Cell0Ds.txt
    ofstream ofile1("Cell0Ds.txt");
    ofile1 << "Id;ShortPath;X;Y;Z\n"; //header 
    
    const MatrixXd &A = P.cell0Ds_coordinates;
    for(unsigned int id=0; id < P.num_cell0Ds; id++)
    {
        ofile1 << defaultfloat << id << ';' << 0 << ';' << scientific << setprecision(16) << A(0,id) << ';' << A(1,id) << ';' << A(2,id) << '\n';
    }
    ofile1.close();

    //Cell1Ds.txt
    ofstream ofile2("Cell1Ds.txt");
    ofile2 << "Id;ShortPath;Origin;End\n"; //header

    const MatrixXi &B = P.cell1Ds_extrema;
    for(unsigned int id=0; id < P.num_cell1Ds; id++)
    {
        ofile2 << id << ';' << 0 << ';' << B(0,id) << ';' << B(1,id) << '\n';
    }
    ofile2.close();

    //Cell2Ds.txt
    ofstream ofile3("Cell2Ds.txt");
    ofile3 << "Id;NumVertices;Vertices;NumEdges;Edges\n"; //header
    
    /*We have to respect the sequential rule described in the pdf project, 
    but during the filling of the 2 matrix for cell2Ds properties, we already check this,
    we can easly read the matrix: cell2Ds_vertices and cell2Ds_edges*/
    const MatrixXi &V = P.cell2Ds_vertices;
    const MatrixXi &E = P.cell2Ds_edges;
    for(unsigned int id=0; id < P.num_cell2Ds; id++)
    {
        ofile3 << id << ';' << V.rows();
        for(unsigned int i=0; i < V.rows(); i++)
        {
            ofile3 << ';' << V(i,id);
        }

        ofile3 << ';' << E.rows();
        for(unsigned int j=0; j < E.rows(); j++)
        {
            ofile3 << ';' << E(j,id);
        }
        ofile3 << '\n';
    }
    ofile3.close();

    //Cell3Ds.txt
    ofstream ofile4("Cell3Ds.txt");
    ofile4 << "Id;NumVertices;Vertices;NumEdges;Edges;NumFaces;Faces\n";

    ofile4 << 0 << ';' << P.num_cell0Ds;
    for(unsigned int id_vert=0; id_vert < P.num_cell0Ds; id_vert++)
        ofile4 << ';' << id_vert;
    
    ofile4 << ';' << P.num_cell1Ds;
    for(unsigned int id_edge=0; id_edge < P.num_cell1Ds; id_edge++)
        ofile4 << ';' << id_edge;

    ofile4 << ';' << P.num_cell2Ds;
    for(unsigned int id_face=0; id_face < P.num_cell2Ds; id_face++)
        ofile4 << ';' << id_face;
    ofile4 << '\n';
    ofile4.close();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
void Visualize_polyhedron(Polyhedron &P)
{
    cout << "num_cell0Ds: " << P.num_cell0Ds << endl;
    cout << "num_cell1Ds: " << P.num_cell1Ds << endl;
    cout << "num_cell2Ds: " << P.num_cell2Ds << endl;
    cout << "num_cell3Ds: " << P.num_cell3Ds << endl;
    cout << endl;

    cout << "cell0Ds_id = [";
    for(unsigned int i=0; i < P.num_cell0Ds; i++)
        cout << ' ' << P.cell0Ds_id[i];
    cout << " ]" << endl;

    cout << "cell1Ds_id = [";
    for(unsigned int i=0; i < P.num_cell1Ds; i++)
        cout << ' ' << P.cell1Ds_id[i];
    cout << " ]" << endl;

    cout << "cell2Ds_id = [";
    for(unsigned int i=0; i < P.num_cell2Ds; i++)
        cout << ' ' << P.cell02Ds_id[i];
    cout << " ]" << endl;

    cout << "cell2Ds_id = [ " << P.cell3D_id << " ]" <<endl;
    cout << endl;

    Eigen::MatrixXd &A = P.cell0Ds_coordinates;
    Eigen::MatrixXi &B = P.cell1Ds_extrema;
    Eigen::MatrixXi &C = P.cell2Ds_vertices;
    Eigen::MatrixXi &D = P.cell2Ds_edges;

    cout<<"Cell0Ds_coordinates: "<<endl;
    cout<<A<<endl;
    cout<<endl;

    cout<<"Cell1Ds_extrema: "<<endl;
    cout<<B<<endl;
    cout<<endl;
    
    cout<<"Cell2Ds_vertices: "<<endl;
    cout<<C<<endl;
    cout<<endl;

    cout<<"Cell2Ds_edges: "<<endl;
    cout<<D<<endl;
    cout<<endl;
}




}