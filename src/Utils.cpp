#include "Utils.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <algorithm> // std::max std::min 
#include <limits>
#include <cmath>
#include<unordered_set>
#include <queue>

namespace PolyhedronLibrary{

void finish_to_fill_struct(Polyhedron &P)
{
    vector<vector<unsigned int>> &B = P.cell2Ds_vertices;
    //Fill cell1Ds_extrema
    //map for sorted pair of vertices --> Id egde 
    std::map<std::pair<unsigned int, unsigned int>, unsigned int> edge_map;
    Eigen::MatrixXi &C = P.cell1Ds_extrema;
    C = MatrixXi::Zero(2, P.num_cell1Ds);
    unsigned int id_edge = 0;
    for (unsigned int f = 0; f < P.num_cell2Ds; ++f) {
        vector<unsigned int> &face = B[f];
        for (unsigned int i = 0; i < face.size(); ++i) {
            unsigned int origin = face[i];
            unsigned int end = face[(i + 1) % face.size()];

            //I sort the vertices to avoid copy (es. (1,0) -> (0,1))
            unsigned int u = std::min(origin, end);
            unsigned int v = std::max(origin, end);
            std::pair<unsigned int, unsigned int> key(u, v);

            //Add the edge only if it already not exist 
            if (edge_map.find(key) == edge_map.end()) {
                edge_map[key] = id_edge;
                C(0, id_edge) = u;
                C(1, id_edge) = v;
                id_edge++;
            }
        }
    }

    //To verify if the number of generated edges is correct
    //assert(id_edge == P.num_cell1Ds);
    
    // Fill cell2Ds_edges
    vector<vector<unsigned int>> &D = P.cell2Ds_edges;
    D.reserve(P.num_cell2Ds);
    D.clear();
    for (unsigned int i = 0; i < B.size(); ++i) {
        vector<unsigned int> &face = B[i];
        vector<unsigned int> edges_f;
        edges_f.reserve(face.size());
        for (unsigned int j = 0; j < face.size(); ++j) {
            unsigned int origin = face[j];
            unsigned int end = face[(j + 1) % face.size()];

            //I sort the vertex to find the correct key in the edge_map
            unsigned int u = std::min(origin, end);
            unsigned int v = std::max(origin, end);
            std::pair<unsigned int, unsigned int> key(u, v);

            //I pick the id up from the edge_map
            edges_f.push_back(edge_map[key]);
        }
        D.push_back(edges_f);
    }

    //Fill cell0Ds_id, cell1Ds_id, cell2Ds_id
    vector<unsigned int> &cell0Ds_id = P.cell0Ds_id;
    vector<unsigned int> &cell1Ds_id = P.cell1Ds_id;
    vector<unsigned int> &cell2Ds_id = P.cell2Ds_id;
    cell0Ds_id.reserve(P.num_cell0Ds);
    cell1Ds_id.reserve(P.num_cell1Ds);
    cell2Ds_id.reserve(P.num_cell2Ds);
    P.cell0Ds_id.clear();
    P.cell1Ds_id.clear();
    P.cell2Ds_id.clear();

    for(unsigned int i=0; i < P.num_cell0Ds; i++) cell0Ds_id.push_back(i);
    for(unsigned int i=0; i < P.num_cell1Ds; i++) cell1Ds_id.push_back(i);
    for(unsigned int i=0; i < P.num_cell2Ds; i++) cell2Ds_id.push_back(i);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Build_platonic_solid(const int p, const int q, Polyhedron &P)
{
    //To establish the corresponding platonic solid {p,q}
    unsigned int V, E, F;
    if(p == 3 && q == 3) //tetrahedron  
    {
        V = 4;
        E = 6;
        F = 4;

        //Fill the matrix cell0Ds_coordinates <--> matrix A
        MatrixXd &A = P.cell0Ds_coordinates; 
        A = MatrixXd::Zero(3, V);
        Eigen::Vector3d v0(0.0, 0.0, 1.732051);
        Eigen::Vector3d v1(1.632993, 0.0, -0.5773503);
        Eigen::Vector3d v2(-0.8164966, 1.414214, -0.5773503);
        Eigen::Vector3d v3(-0.8164966, -1.414214, -0.5773503);
        A.col(0) = v0;
        A.col(1) = v1;
        A.col(2) = v2;
        A.col(3) = v3;

        //Fill the cell2Ds_vertices <--> B 
        vector<vector<unsigned int>> &B = P.cell2Ds_vertices;
        B.reserve(F);
        vector<unsigned int> f0 = {0, 1, 2};
        vector<unsigned int> f1 = {0, 2, 3};
        vector<unsigned int> f2 = {0, 3, 1};
        vector<unsigned int> f3 = {1, 3, 2};
        B.push_back(f0);
        B.push_back(f1);
        B.push_back(f2);
        B.push_back(f3);
    }
    else if ((p == 3 && q == 4) || (p == 4 && q == 3 )) //octahedron
    {
        V = 6;
        E = 12;
        F = 8;
        
        //Fill the matrix cell0Ds_coordinates <--> matrix A
        MatrixXd &A = P.cell0Ds_coordinates; 
        A = MatrixXd::Zero(3, V);
        Eigen::Vector3d v0(0.0, 0.0, 1.414214);
        Eigen::Vector3d v1(1.414214, 0.0, 0.0);
        Eigen::Vector3d v2(0.0, 1.414214, 0.0);
        Eigen::Vector3d v3(-1.414214, 0.0, 0.0);
        Eigen::Vector3d v4(0.0, -1.414214, 0.0);
        Eigen::Vector3d v5(0.0, 0.0, -1.414214);
        A.col(0) = v0;
        A.col(1) = v1; 
        A.col(2) = v2;
        A.col(3) = v3;
        A.col(4) = v4;
        A.col(5) = v5;

        //Fill the cell2Ds_vertices <--> B 
        vector<vector<unsigned int>> &B = P.cell2Ds_vertices;
        B.reserve(F);
        vector<unsigned int> f0 = {0, 1, 2};
        vector<unsigned int> f1 = {0, 2, 3};
        vector<unsigned int> f2 = {0, 3, 4};
        vector<unsigned int> f3 = {0, 4, 1};
        vector<unsigned int> f4 = {1, 4, 5};
        vector<unsigned int> f5 = {1, 5, 2};
        vector<unsigned int> f6 = {2, 5, 3};
        vector<unsigned int> f7 = {3, 5, 4};
        B.push_back(f0);
        B.push_back(f1);
        B.push_back(f2);
        B.push_back(f3);
        B.push_back(f4);
        B.push_back(f5);
        B.push_back(f6);
        B.push_back(f7);
    }
    else if ((p == 3 && q == 5) || (p == 5 && q == 3)) //icosahedron
    {
        V = 12;
        E = 30;
        F = 20;
        
        //Fill the matrix cell0Ds_coordinates <--> matrix A
        MatrixXd &A = P.cell0Ds_coordinates; 
        A = MatrixXd::Zero(3, V);
        Eigen::Vector3d v0(0.0, 0.0, 1.175571);
        Eigen::Vector3d v1(1.051462, 0.0, 0.5257311);
        Eigen::Vector3d v2(0.3249197, 1.0, 0.5257311);
        Eigen::Vector3d v3(-0.8506508, 0.618034, 0.5257311);
        Eigen::Vector3d v4(-0.8506508, -0.618034, 0.5257311);
        Eigen::Vector3d v5(0.3249197, -1.0, 0.5257311);
        Eigen::Vector3d v6(0.8506508, 0.618034, -0.5257311);
        Eigen::Vector3d v7(0.8506508, -0.618034, -0.5257311);
        Eigen::Vector3d v8(-0.3249197, 1.0, -0.5257311);
        Eigen::Vector3d v9(-1.051462, 0.0, -0.5257311);
        Eigen::Vector3d v10(-0.3249197, -1.0, -0.5257311);
        Eigen::Vector3d v11(0.0, 0.0, -1.175571);
        A.col(0) = v0;
        A.col(1) = v1;
        A.col(2) = v2;
        A.col(3) = v3;
        A.col(4) = v4;
        A.col(5) = v5;
        A.col(6) = v6;
        A.col(7) = v7;
        A.col(8) = v8;
        A.col(9) = v9;
        A.col(10) = v10;
        A.col(11) = v11;

        //Fill the cell2Ds_vertices <--> B 
        vector<vector<unsigned int>> &B = P.cell2Ds_vertices;
        B.reserve(F);
        vector<unsigned int> f0 = {0, 1, 2};
        vector<unsigned int> f1 = {0, 2, 3};
        vector<unsigned int> f2 = {0, 3, 4};
        vector<unsigned int> f3 = {0, 4, 5};
        vector<unsigned int> f4 = {0, 5, 1};
        vector<unsigned int> f5 = {1, 5, 7};
        vector<unsigned int> f6 = {1, 7, 6};
        vector<unsigned int> f7 = {1, 6, 2};
        vector<unsigned int> f8 = {2, 6, 8};
        vector<unsigned int> f9 = {2, 8, 3};
        vector<unsigned int> f10 = {3, 8, 9};
        vector<unsigned int> f11 = {3, 9, 4};
        vector<unsigned int> f12 = {4, 9, 10};
        vector<unsigned int> f13 = {4, 10, 5};
        vector<unsigned int> f14 = {5, 10, 7};
        vector<unsigned int> f15 = {6, 7, 11};
        vector<unsigned int> f16 = {6, 11, 8};
        vector<unsigned int> f17 = {7, 10, 11};
        vector<unsigned int> f18 = {8, 11, 9};
        vector<unsigned int> f19 = {9, 11, 10};
        B.push_back(f0);
        B.push_back(f1);
        B.push_back(f2);
        B.push_back(f3);
        B.push_back(f4);
        B.push_back(f5);
        B.push_back(f6);
        B.push_back(f7);
        B.push_back(f8);
        B.push_back(f9);
        B.push_back(f10);
        B.push_back(f11);
        B.push_back(f12);
        B.push_back(f13);
        B.push_back(f14);
        B.push_back(f15);
        B.push_back(f16);
        B.push_back(f17);
        B.push_back(f18);
        B.push_back(f19);

    }
    else return false;
    
    P.num_cell0Ds = V;
    P.num_cell1Ds = E;
    P.num_cell2Ds = F;

    //Call the function to finish to fill the polyhedron struct  
    finish_to_fill_struct(P);

    return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
void Export_polyhedron(const Polyhedron &P)
{
    //Cell0Ds.txt
    ofstream ofile1("Cell0Ds.txt");
    ofile1 << "Id;X;Y;Z\n"; //header 
    
    const MatrixXd &A = P.cell0Ds_coordinates;
    for(unsigned int id=0; id < P.num_cell0Ds; id++)
        ofile1 << defaultfloat << id << ';' << scientific << setprecision(16) << A(0,id) << ';' << A(1,id) << ';' << A(2,id) << '\n';
    ofile1.close();

    //Cell1Ds.txt
    ofstream ofile2("Cell1Ds.txt");
    ofile2 << "Id;Origin;End\n"; //header

    const MatrixXi &B = P.cell1Ds_extrema;
    for(unsigned int id=0; id < P.num_cell1Ds; id++)
        ofile2 << id << ';' << B(0,id) << ';' << B(1,id) << '\n';
    ofile2.close();

    //Cell2Ds.txt
    ofstream ofile3("Cell2Ds.txt");
    ofile3 << "Id;NumVertices;Vertices;NumEdges;Edges\n"; //header
    
    /*We have to respect the sequential rule described in the pdf project, 
    but during the filling of the 2 matrix for cell2Ds properties, we already check this,
    we can easly read the matrix: cell2Ds_vertices and cell2Ds_edges*/
    const vector<vector<unsigned int>> &V = P.cell2Ds_vertices;
    const vector<vector<unsigned int>> &E = P.cell2Ds_edges;
    for(unsigned int id=0; id < P.num_cell2Ds; id++)
    {
        ofile3 << id << ';' << V[id].size();
        for(unsigned int i=0; i < V[id].size(); i++)
            ofile3 << ';' << V[id][i];
        
        ofile3 << ';' << E[id].size();
        for(unsigned int j=0; j < E[id].size(); j++)
            ofile3 << ';' << E[id][j];

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

    cout<<"Cell0Ds.txt, Cell1Ds.txt, Cell2Ds.txt, Cell3Ds.txt files created"<<endl;
    cout<<endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
void Visualize_polyhedron(const Polyhedron &P)
{
    cout << "num_cell0Ds: " << P.num_cell0Ds << " (vertices V)" << endl;
    cout << "num_cell1Ds: " << P.num_cell1Ds << " (edges E)" << endl;
    cout << "num_cell2Ds: " << P.num_cell2Ds << " (faces F)" << endl;
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


    const Eigen::MatrixXd &A = P.cell0Ds_coordinates;
    const Eigen::MatrixXi &B = P.cell1Ds_extrema;
    const vector<vector<unsigned int>> &C = P.cell2Ds_vertices;
    const vector<vector<unsigned int>> &D = P.cell2Ds_edges;

    cout<<"Cell0Ds_coordinates: "<<endl;
    cout<<A<<endl;
    cout<<endl;

    cout<<"Cell1Ds_extrema: "<<endl;
    cout<<B<<endl;
    cout<<endl;
    
    cout<<"Cell2Ds_vertices: "<<endl;
    for(unsigned int i=0; i < C.size(); i++){
        cout<<'f'<<i<<": ";
        const vector<unsigned int> &face = C[i];
        for(unsigned int j=0; j < face.size(); j++){
            cout<<face[j]<<' ';
        }
        cout<<endl;
    }
    cout<<endl;

    cout<<"Cell2Ds_edges: "<<endl;
    for(unsigned int i=0; i < D.size(); i++){
        cout<<'f'<<i<<": ";
        const vector<unsigned int> &face = D[i];
        for(unsigned int j=0; j < face.size(); j++){
            cout<<face[j]<<' ';
        }
        cout<<endl;
    }
    cout<<endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
pair<vector<Eigen::Vector3d>, vector<Eigen::Vector3i>> classI_basic_step(const Eigen::Vector3d &A, const Eigen::Vector3d &B, const Eigen::Vector3d &C, const int b)
{
    vector<Vector3d> vertices;
    vector<Vector3i> triangles;

    //1. Cause the vertices in baricentric coordinates
    for (int u = 0; u <= b; u++) {
        for (int v = 0; v <= (b - u); v++) {
            int w = b - u - v;
            Vector3d punto = (u * A + v * B + w * C) / b; 
            vertices.push_back(punto); //C is the first one, B is at position b, A is the last one
        }
    }

    //2. Cause the triangles
    for (int u = 0; u < b; u++) 
    {
        for (int v = 0; v < (b - u); v++) 
        {
            // Compute the index of the vertices for this cell 
            int current_row_start = u * (b + 1) - (u * (u - 1)) / 2;
            int next_row_start = (u + 1) * (b + 1) - (u * (u + 1)) / 2;

            int idx0 = current_row_start + v;          //(u,v)
            int idx1 = next_row_start + v;             //(u+1,v)
            int idx2 = current_row_start + v + 1;      //(u,v+1)
            
            // Add the "superior" triangle
            triangles.emplace_back(idx0, idx1, idx2);

            // Add the "inferior" triangle if we are not at the edge
            if (v < (b - u - 1)) 
            {
                int idx3 = next_row_start + v + 1;     //(u+1, v+1)
                triangles.emplace_back(idx1, idx3, idx2);
            }
        }
    }

    return {vertices, triangles};
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
void Triangulate(Polyhedron &P, const int b, const int c)
{   
    vector<Eigen::Vector3i> new_faces;

    std::map<unsigned int, Eigen::Vector3d> map_id_vert_global; //map {key : val} --> {id : vertex}
    unsigned int new_id = 0;

    for(unsigned int i=0; i < P.num_cell2Ds; i++)
    {
        vector<unsigned int> &face = P.cell2Ds_vertices[i];
        unsigned int id_A = face[0];
        unsigned int id_B = face[1];
        unsigned int id_C = face[2];
        
        Eigen::Vector3d A = P.cell0Ds_coordinates.col(id_A);
        Eigen::Vector3d B = P.cell0Ds_coordinates.col(id_B);
        Eigen::Vector3d C = P.cell0Ds_coordinates.col(id_C);

        vector<Eigen::Vector3d> vertices;
        vector<Eigen::Vector3i> triangles;
        if((b >= 1 && c == 0) || (b == 0 && c >=1)) //Class I (geodetic polyhedron)
        {
            if(b != 0){
                std::pair<vector<Eigen::Vector3d>, vector<Eigen::Vector3i>> result = classI_basic_step(A, B, C, b);
                vertices = result.first;
                triangles = result.second;
            } 
            else{
                std::pair<vector<Eigen::Vector3d>, vector<Eigen::Vector3i>> result = classI_basic_step(A, B, C, c);
                vertices = result.first;
                triangles = result.second;
            } 
        }
        else if (b == c) //Class II
        {
            std::pair result = classII_basic_step(A, B, C, b);
            vertices = result.first;
            triangles = result.second;
        }
        

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
    unsigned int V = map_id_vert_global.size();
    unsigned int F = new_faces.size();
    unsigned int E = V + F - 2; //Eulero's formula
    
    P.num_cell0Ds = V;
    P.num_cell1Ds = E;
    P.num_cell2Ds = F;

    //Fill cell0Ds_coordinates
    Eigen::MatrixXd &A = P.cell0Ds_coordinates;
    A = MatrixXd::Zero(3, V);
    for(auto &[k, v] : map_id_vert_global) A.col(k) = v;

    //Fill cell2Ds_vertices
    vector<vector<unsigned int>> &B = P.cell2Ds_vertices;
    B.reserve(F);
    B.clear();
    for(unsigned int i=0; i < new_faces.size(); i++) 
    {
        Eigen::Vector3i &f = new_faces[i];
        vector<unsigned int> face;
        face.reserve(f.size());
        for(unsigned int j=0; j < f.size(); j++) face.push_back(f[j]);
        B.push_back(face);
    }
    
    finish_to_fill_struct(P);
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
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<unsigned int> cycled_face_for_dual(vector<unsigned int>& face_new, const Eigen::MatrixXd& coord) {
    vector<unsigned int> ordered;
    if (face_new.empty()) return ordered;

    //Compute the baricenter of the face
    Vector3d baricenter = Vector3d::Zero();
    for (auto id : face_new) {
        baricenter += coord.col(id);
    }
    baricenter /= face_new.size();

    //Calculate the normal to the face plane using three points
    Vector3d p0 = coord.col(face_new[0]);
    Vector3d p1 = coord.col(face_new[1]);
    Vector3d p2 = coord.col(face_new[2]);

    Vector3d v1 = p1 - p0;
    Vector3d v2 = p2 - p0;
    Vector3d normal = v1.cross(v2).normalized();

    //Reference direction from the center of gravity to the first point
    Vector3d ref_dir = (p0 - baricenter).normalized();

    //Local coordinate system on the plane
    Vector3d u_axis = ref_dir;
    Vector3d v_axis = normal.cross(u_axis).normalized();

    //Calculate the angle for each point-centroid vector projected onto the plane
    vector<pair<double, unsigned int>> angle_id_pairs;
    for (auto id : face_new) {
        Vector3d vec = coord.col(id) - baricenter;
        double x = vec.dot(u_axis);
        double y = vec.dot(v_axis);
        double angle = atan2(y, x); //arctangent gives me an angle in radians
        
        //Correct the angle to be between 0 and 2pi
        if (angle < 0) angle += 2 * M_PI; //M_PI = pi
        
        angle_id_pairs.emplace_back(angle, id); //using this I can forget to specify that I have to add a pair to the vector it's automatic
    }

    //Sort points based on calculated angle (std::sort --> growing order)
    sort(angle_id_pairs.begin(), angle_id_pairs.end());

    //I extract the ordered IDs
    for (const auto& pair : angle_id_pairs) {
        ordered.push_back(pair.second);
    }

    return ordered;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
void Dualize(Polyhedron &P)
{   
    unsigned int old_F = P.num_cell2Ds;
    unsigned int old_V = P.num_cell0Ds;
    P.num_cell0Ds = old_F;
    P.num_cell2Ds = old_V;


    //mapping the old_faces to new_vertices
    Eigen::MatrixXd new_vertices = MatrixXd::Zero(3, old_F); 

    Eigen::MatrixXd old_vertices = P.cell0Ds_coordinates;
    vector<vector<unsigned int>> &old_faces = P.cell2Ds_vertices;
    for(unsigned int i=0; i < old_F; i++)
    {
        vector<unsigned int> &face = old_faces[i];

        //compute the baricenter of the face
        Eigen::Vector3d S = Vector3d::Zero();
        for(unsigned int j=0; j < face.size(); j++)
        {
            unsigned int id = face[j];
            Eigen::Vector3d z = old_vertices.col(id);
            S += z;   
        }
        Eigen::Vector3d baricenter = S / face.size();

        new_vertices.col(i) = baricenter;
    } 

    //Overload and fill cell0Ds_coordinates
    Eigen::MatrixXd &A = P.cell0Ds_coordinates;
    A = MatrixXd::Zero(3, old_F);
    A = new_vertices;
    
    
    //mapping the old_vertices to new_faces
    vector<vector<unsigned int>> new_faces_by_vert;
    new_faces_by_vert.reserve(old_V);
    
    for(unsigned int id=0; id < old_V; id++)
    {   
        vector<unsigned int> candidate_face;
        for(unsigned int j=0; j < old_F; j++)
        {
            vector<unsigned int> &f = P.cell2Ds_vertices[j];
            std::set<unsigned int> face_set;
            for(unsigned int k=0; k < f.size(); k++) face_set.insert(f[k]); //fill teh set

            if(face_set.count(id) > 0) candidate_face.push_back(j);
        }
        //the candidate face must be cyclic
        vector<unsigned int> cyclyc = cycled_face_for_dual(candidate_face, P.cell0Ds_coordinates);
        new_faces_by_vert.push_back(cyclyc);  
    }

    //Overload and fill cell2Ds_vertices
    vector<vector<unsigned int>> &B = P.cell2Ds_vertices;
    B.reserve(old_V);
    B.clear();
    B = new_faces_by_vert;

    finish_to_fill_struct(P);   
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
pair<vector<Eigen::Vector3d>, vector<Eigen::Vector3i>> classII_basic_step(const Eigen::Vector3d &A, const Eigen::Vector3d &B, const Eigen::Vector3d &C, const int b)
{
    vector<Vector3d> new_vertices;
    vector<Vector3i> new_triangles;

    auto [base_vertices, base_triangles] = classI_basic_step(A, B, C, b);
    
    map<unsigned int, Vector3d> id_points_map;
    unordered_set<unsigned int> midpoints_ids; // To track midpoint IDs
    unsigned int id = 0;
    for(unsigned int i=0; i < base_triangles.size(); i++)
    {
        Vector3i T = base_triangles[i];

        Vector3d a = base_vertices[T[0]];
        Vector3d b_vert = base_vertices[T[1]];
        Vector3d c = base_vertices[T[2]];

        Vector3d M_ab = 0.5 * (a + b_vert);
        Vector3d M_bc = 0.5 * (b_vert + c);
        Vector3d M_ca = 0.5 * (c + a);

        Vector3d centroid = (a + b_vert + c) / 3;

        array<Vector3d, 7> local_points = {a, b_vert, c, M_ab, M_bc, M_ca, centroid};
        array<int, 7> id_local_points = {-1, -1, -1, -1, -1, -1, -1};
        for(unsigned int p=0; p < 7; p++)
        {
            Vector3d &P = local_points[p];
            bool to_add = true;
            for(auto &[key, val] : id_points_map)
                if(P.isApprox(val)){
                    to_add = false;
                    id_local_points[p] = key;
                }
            if(to_add){
                id_points_map[id] = P;
                id_local_points[p] = id;
                // Check if this point is a midpoint (p=3,4,5)
                if (p >=3 && p <=5) {
                    midpoints_ids.insert(id);
                }
                id++;
            } else {
                // Check if existing point is a midpoint (p=3,4,5)
                if (p >=3 && p <=5) {
                    midpoints_ids.insert(id_local_points[p]);
                }
            }
        }

        int id_a = id_local_points[0];
        int id_b = id_local_points[1];
        int id_c = id_local_points[2];
        int id_M_ab = id_local_points[3];
        int id_M_bc = id_local_points[4];
        int id_M_ca = id_local_points[5];
        int id_centroid = id_local_points[6];
        
        Vector3i t1 {id_a, id_centroid, id_M_ab};
        Vector3i t2 {id_a, id_centroid, id_M_ca};
        Vector3i t3 {id_c, id_centroid, id_M_ca};
        Vector3i t4 {id_c, id_centroid, id_M_bc};
        Vector3i t5 {id_b, id_centroid, id_M_bc};
        Vector3i t6 {id_b, id_centroid, id_M_ab}; 
        
        new_triangles.emplace_back(t1);
        new_triangles.emplace_back(t2);
        new_triangles.emplace_back(t3);
        new_triangles.emplace_back(t4);
        new_triangles.emplace_back(t5);
        new_triangles.emplace_back(t6);
    }

    //Fill new_vertices
    new_vertices.reserve(id_points_map.size());
    for(auto &[k, v]: id_points_map) new_vertices.push_back(v);

    //Now I have to remove some vertices
    //Process midpoints to replace four triangles with two --> I replace four triangles if they have a vertex in common and 
    //the vertex it is not a baricenter of the triangles belong to class 1
    unordered_map<unsigned int, vector<unsigned int>> midpoint_to_triangles;
    for (unsigned int i = 0; i < new_triangles.size(); ++i) {
        const Vector3i& tri = new_triangles[i];
        for (int j = 0; j < 3; ++j) {
            unsigned int v = tri[j];
            if (midpoints_ids.count(v)) {
                midpoint_to_triangles[v].push_back(i);
            }
        }
    }

    unordered_set<unsigned int> triangles_to_remove;
    vector<Vector3i> triangles_to_add;

    for (const auto& entry : midpoint_to_triangles) {
        unsigned int midpoint_id = entry.first;
        const vector<unsigned int>& triangles = entry.second;
        if (triangles.size() != 4) continue;

        unordered_set<unsigned int> other_vertices;
        for (auto tri_idx : triangles) {
            const Vector3i& tri = new_triangles[tri_idx];
            for (int j = 0; j < 3; ++j) {
                unsigned int v = tri[j];
                if (v != midpoint_id) {
                    other_vertices.insert(v);
                }
            }
        }
        if (other_vertices.size() != 4) continue;

        vector<unsigned int> original_vs, centroids;
        for (auto v : other_vertices) {
            const Vector3d& vertex = new_vertices[v];
            bool is_original = false;
            for (const auto& bv : base_vertices) {
                if (vertex.isApprox(bv)) {
                    is_original = true;
                    break;
                }
            }
            if (is_original) {
                original_vs.push_back(v);
            } else {
                centroids.push_back(v);
            }
        }

        if (original_vs.size() != 2 || centroids.size() != 2) continue;

        unsigned int a = original_vs[0];
        unsigned int b = original_vs[1];
        unsigned int C1 = centroids[0];
        unsigned int C2 = centroids[1];

        Vector3i new_tri1(C1, a, C2);
        Vector3i new_tri2(C2, b, C1);
        triangles_to_add.push_back(new_tri1);
        triangles_to_add.push_back(new_tri2);

        for (auto tri_idx : triangles) {
            triangles_to_remove.insert(tri_idx);
        }
    }

        vector<Vector3i> final_triangles;
    final_triangles.reserve(new_triangles.size() - triangles_to_remove.size() + triangles_to_add.size());
    for (unsigned int i = 0; i < new_triangles.size(); ++i) {
        if (!triangles_to_remove.count(i)) {
            final_triangles.push_back(new_triangles[i]);
        }
    }
    final_triangles.insert(final_triangles.end(), triangles_to_add.begin(), triangles_to_add.end());

    //Step 1: Identifying used vertices
    unordered_set<unsigned int> used_vertices;
    for(const auto& tri : final_triangles) {
        used_vertices.insert(tri[0]);
        used_vertices.insert(tri[1]);
        used_vertices.insert(tri[2]);
    }

    //Step 2: Reconstruction of clean vertices list
    vector<Vector3d> cleaned_vertices;
    unordered_map<unsigned int, unsigned int> id_remap;
    unsigned int new_id = 0;
    
    for(const auto& [old_id, vertex] : id_points_map) {
        if(used_vertices.count(old_id)) {
            id_remap[old_id] = new_id;
            cleaned_vertices.push_back(vertex);
            new_id++;
        }
    }

    //Step 3: Updating indices in triangles
    vector<Vector3i> cleaned_triangles;
    for(const auto& tri : final_triangles) {
        cleaned_triangles.emplace_back(
            id_remap[tri[0]], 
            id_remap[tri[1]], 
            id_remap[tri[2]]
        );
    }

    return {cleaned_vertices, cleaned_triangles};
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<unsigned int> Short_path(const Polyhedron &P, const int id_D, const int id_A)
{
    //Make the adjaceny list of the graph
    const Eigen::MatrixXi &edges = P.cell1Ds_extrema;
    
    vector<vector<unsigned int>> adj_list;
    adj_list.reserve(P.num_cell0Ds);
    for(unsigned int id=0; id < P.num_cell0Ds; id++){
        vector<unsigned int> adj_to_id;
        for(unsigned int j=0; j < edges.cols(); j++){
            Eigen::Vector2i edge = edges.col(j);
            unsigned int origin = edge[0];
            unsigned int end = edge[1];

            
            if(id == origin) adj_to_id.push_back(end);
            if(id == end) adj_to_id.push_back(origin);
        }
        adj_list.push_back(adj_to_id);
    }


    //Implementation of the BFS to look for short path between id_D and id_A
    vector<unsigned int> path;
    //Trivial case
    if(id_D == id_A){
        path.reserve(1);
        path.push_back(id_D);
        return path;
    }

    //Other cases
    const unsigned int N = P.num_cell0Ds;
    vector<bool> visited(N, false);
    vector<unsigned int> parent(N, N); //N = pivot value 
    std::queue<unsigned int> q;

    q.push(id_D);
    visited[id_D] = true;
    parent[id_D] = N;

    bool stop = false;
    while(!q.empty() || !stop){
        unsigned int current = q.front();
        q.pop();

        //Explore the neighbours of current
        for(unsigned int neighbour : adj_list[current]){
            if(!visited[neighbour]){
                visited[neighbour] = true;
                parent[neighbour] = current;
                q.push(neighbour);

                //Stop if neighbours is id_A (destination)
                if(neighbour == id_A)
                {
                    stop = true;
                    break;
                }     
            }
        }
    }

    //Reconstruct the path
    unsigned int node = id_A;

    while(node != N){
        path.push_back(node);
        node = parent[node];
    }

    std::reverse(path.begin(), path.end());

    ////////
    cout<<"Path: ";
    for(int i=0; i<path.size(); i++){
        cout<<path[i]<<' ';
    }
    cout<<endl;

    //Required output
    cout<<"The shortest path that links "<<id_D<<" and "<<id_A<<" is "<<(path.size()-1)<<" sides long"<<endl;

    double length = 0.0;
    for(unsigned int i=0; i < path.size() - 1; i++){
        unsigned int id_U = path[i];
        unsigned int id_V = path[i+1];

        Eigen::Vector3d U = P.cell0Ds_coordinates.col(id_U);
        Eigen::Vector3d V = P.cell0Ds_coordinates.col(id_V);
        length += (U-V).norm();

    }
    cout<<"The shortest path that links "<<id_D<<" and "<<id_A<<" is "<<length<<" long"<<endl<<endl;
    return path;
}








}

