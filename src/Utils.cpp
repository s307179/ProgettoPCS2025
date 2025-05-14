#include "Utils.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <set> 

namespace PolyhedronLibrary{

bool Import_platonic_solid(const unsigned int p, const unsigned int q, Polyhedron &P)
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
    cout << "num_cell0Ds: " << P.num_cell0Ds << " (vertices V)" << endl;
    cout << "num_cell1Ds: " << P.num_cell1Ds << " (edges E)" << endl;
    cout << "num_cell2Ds: " << P.num_cell2Ds << " (faces F)" << endl;
    cout << "num_cell3Ds: " << P.num_cell3Ds << " (polyhedron P)" << endl;
    cout << endl;

    cout << "cell0Ds_id = [";
    for(unsigned int i=0; i < P.num_cell0Ds; i++)
        cout << ' ' << P.cell0Ds_id[i];
    cout << " ]" << " (vertices id)" << endl;

    cout << "cell1Ds_id = [";
    for(unsigned int i=0; i < P.num_cell1Ds; i++)
        cout << ' ' << P.cell1Ds_id[i];
    cout << " ]" << " (edges id)" << endl;

    cout << "cell2Ds_id = [";
    for(unsigned int i=0; i < P.num_cell2Ds; i++)
        cout << ' ' << P.cell2Ds_id[i];
    cout << " ]" << " (faces id)" << endl;

    cout << "cell3D_id = [ " << P.cell3D_id << " ]" << " (polyhedron id)" <<endl;
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


////////////////////////////////////////////////////////////////////////////////////////////////////////
pair<vector<Eigen::Vector3d>, vector<Eigen::Vector3i>> Triangulation_basic_step(const Eigen::Vector3d &A, const Eigen::Vector3d &B, const Eigen::Vector3d &C, const unsigned int b)
{
    vector<Vector3d> vertici;
    vector<Vector3i> triangoli;

    // 1. Genera i vertici in coordinate baricentriche
    for (int u = 0; u <= b; u++) {
        for (int v = 0; v <= (b - u); v++) {
            int w = b - u - v;
            Vector3d punto = (u * A + v * B + w * C) / b; 
            vertici.push_back(punto); //C viene appeso per primo, B viene appeso in posizione b, A viene appeso in ultima posizone
        }
    }

    // 2. Genera i triangoli
    for (int u = 0; u < b; u++) 
    {
        for (int v = 0; v < (b - u); v++) 
        {
            // Calcola gli indici dei vertici per questa cella
            int current_row_start = u * (b + 1) - (u * (u - 1)) / 2;
            int next_row_start = (u + 1) * (b + 1) - (u * (u + 1)) / 2;

            int idx0 = current_row_start + v;          //(u, v)
            int idx1 = next_row_start + v;             //(u+1, v)
            int idx2 = current_row_start + v + 1;      //(u, v+1)
            
            // Aggiungi il triangolo "superiore"
            triangoli.emplace_back(idx0, idx1, idx2);

            // Aggiungi il triangolo "inferiore" se non siamo al bordo
            if (v < (b - u - 1)) 
            {
                int idx3 = next_row_start + v + 1;     //(u+1, v+1)
                triangoli.emplace_back(idx1, idx3, idx2);
            }
        }
    }

    return {vertici, triangoli};
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
void ClassI_polyhedron(Polyhedron &P, const unsigned int b, const unsigned int q)
{   
    vector<Eigen::Vector3i> new_faces;

    std::map<unsigned int, Eigen::Vector3d> map_id_vert_global; //mappa {id, vector} --> {id : vertex}
    unsigned int new_id = 0;

    for(unsigned int i=0; i < P.num_cell2Ds; i++)
    {
        Eigen::Vector3i face = P.cell2Ds_vertices.col(i);
        unsigned int id_A = face[0];
        unsigned int id_B = face[1];
        unsigned int id_C = face[2];
        
        Eigen::Vector3d A = P.cell0Ds_coordinates.col(id_A);
        Eigen::Vector3d B = P.cell0Ds_coordinates.col(id_B);
        Eigen::Vector3d C = P.cell0Ds_coordinates.col(id_C);


        std::pair [vertices, triangles] = Triangulation_basic_step(A, B, C, b);

        //build of map from local id to global id
        for(unsigned int v=0; v < vertices.size(); v++)
        {
            Eigen::Vector3d &vec = vertices[v];
            bool to_add = true;
            for(const auto &[key, val] : map_id_vert_global)
                if(vec.isApprox(val)) to_add = false;
            
            if(to_add)
            {
                map_id_vert_global.insert({new_id, vec});
                new_id++; 
            }  
        }

        //reindex the local face using the map
        for(unsigned int t=0; t < triangles.size(); t++)
        {
            Eigen::Vector3i &tri = triangles[t];
            unsigned int idx0, idx1, idx2; //id vertices of the triangle
            idx0 = tri[0];
            idx1 = tri[1];
            idx2 = tri[2];

            Eigen::Vector3d v1, v2, v3;
            v1 = vertices[idx0];
            v2 = vertices[idx1];
            v3 = vertices[idx2];

            for(const auto &[key, val] : map_id_vert_global)
            {
                if(v1.isApprox(val)) tri[0] = key;
                if(v2.isApprox(val)) tri[1] = key;
                if(v3.isApprox(val)) tri[2] = key;
            }
        }

        //store the faces in new_faces
        for(unsigned int t=0; t < triangles.size(); t++)
            new_faces.push_back(triangles[t]);
    }

    
    //overload the polyhedron struct with the new data
    unsigned int T = b*b;
    unsigned int V, E, F;
    if(q == 3)
    {
        V = 2*T + 2;
        E = 6*T;
        F = 4*T;
    }
    else if(q == 4)
    {
        V = 4*T + 2;
        E = 12*T;
        F = 8*T;
    }
    else
    {
        V = 10*T + 2;
        E = 30*T;
        F = 20*T;
    }

    P.num_cell0Ds = V;
    P.num_cell1Ds = E;
    P.num_cell2Ds = F;

    P.cell0Ds_id.reserve(V);
    P.cell1Ds_id.reserve(E);
    P.cell2Ds_id.reserve(F);
    P.cell0Ds_id.clear();
    P.cell1Ds_id.clear();
    P.cell2Ds_id.clear();
    for(unsigned int i=0; i < V ; i++) P.cell0Ds_id.push_back(i);
    for(unsigned int i=0; i < E; i++) P.cell1Ds_id.push_back(i);
    for(unsigned int i=0; i < F; i++) P.cell2Ds_id.push_back(i);

    //Fill cell0Ds_coordinates
    Eigen::MatrixXd &A = P.cell0Ds_coordinates;
    A = MatrixXd::Zero(3, V);
    for(auto &[k, v] : map_id_vert_global) A.col(k) = v;

    //Fill cell2Ds_vertices
    Eigen::MatrixXi &B = P.cell2Ds_vertices;
    B = MatrixXi::Zero(3, F);
    for(unsigned int i=0; i < new_faces.size(); i++) B.col(i) = new_faces[i];

    //Fill cell1Ds_extrema
    Eigen::MatrixXi &C = P.cell1Ds_extrema;
    C = MatrixXi::Zero(2, E);
    map<unsigned int, set<unsigned int>> check_map; // {edge_id, set{origin, end}}
    unsigned int id_edge = 0;
    for(size_t f=0; f < F; f++)
    {
        Eigen::Vector3i &face = new_faces[f];
        for(unsigned int i=0; i < 3; i++)
        {
            unsigned int origin = face[i];
            unsigned int end = face[(i+1)%3];

            bool to_add = true;
            for(auto &[key, val] : check_map)
            {
                if(val.count(origin) + val.count(end) > 1) to_add = false;
            }

            if(to_add)
            {   
                set<unsigned int> s = {origin, end};
                check_map.insert({id_edge, s});
                id_edge++;
            }
        }
    }

    for(auto &[id, o_and_e] : check_map)
    {
        unsigned int origin = *(o_and_e.begin());
        unsigned int end = *(std::next(o_and_e.begin()));
        C(0,id) = origin;
        C(1,id) = end;
    }
    
    //Fill MatrixXi cell2Ds_edges
    Eigen::MatrixXi &D = P.cell2Ds_edges;
    D = MatrixXi::Zero(3, F);
    for(unsigned int i=0; i < B.cols(); i++)
    {
        VectorXi face = B.col(i); //matrix B --> cell2Ds_vertices
        for(unsigned int j=0; j < 3; j++)
        {
            unsigned int u = face[j];
            unsigned int v = face[(j+1)%3];

            for(unsigned int k=0; k < C.cols(); k++) //matrix C --> cell1Ds_extrema
            {
                VectorXi edge = C.col(k);
                unsigned int origin = edge[0];
                unsigned int end = edge[1];

                if ((u == origin && v == end) || (u == end && v == origin)) 
                {
                    D(j,i) = k;
                }    
            }
        }
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////

void project_points_onto_sphere(Polyhedron &P)
	 {
	// necessary info of Polyhedron saved
	MatrixXd &S = P.cell0Ds_coordinates;
	int n = P.num_cell0Ds; 
	// Calculate barycenter coordinated 
	VectorXd barycenter =  VectorXd::Zero(3);
	for (int i =0; i<n; i++)
	{for (int j =0; j<3; j++)
		barycenter(j) += S(j,i)/n;
	 }
	 // move points with respect to the barycenter
	 for(int h = 0; h<n;h++)
	 {
		 for (int j =0; j<3; j++)
		 {S(j,h) -=barycenter(j);}
 
	}
	
	// project on the sphere
	for (int i = 0; i<n; i++)
	{double norm = sqrt(S(0,i)*S(0,i)+S(1,i)*S(1,i)+S(2,i)*S(2,i));
		for (int j=0;j<3;j++)
		{S(j,i) /= abs(norm);}
	}	;
	};



}