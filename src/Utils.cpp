#include "Utils.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <algorithm> // std::max std::min 
#include <limits>
#include <cmath>

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
    assert(id_edge == P.num_cell1Ds);
    
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
bool Import_platonic_solid(unsigned int p, unsigned int q, Polyhedron &P)
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
    else if ((p == 3 && q == 4) || (p == 4 && q == 3 )) 
    {
        poly_name = "octahedron";
        V = 6;
        E = 12;
        F = 8;
        
        p = 3;
        q = 4;
    }
    else if ((p == 3 && q == 5) || (p == 5 && q == 3)) 
    {
        poly_name = "icosahedron";
        V = 12;
        E = 30;
        F = 20;
        
        p = 3;
        q = 5;
    }
    else return false;
    
    P.num_cell0Ds = V;
    P.num_cell1Ds = E;
    P.num_cell2Ds = F;

    ifstream ifile("./" + poly_name + ".txt");
    if(ifile.fail()) return false;

    //Cell0Ds properties (vertices)
    //Fill the matrix cell0Ds_coordinates <--> matrix A
    MatrixXd &A = P.cell0Ds_coordinates; 
    A = MatrixXd::Zero(3, V);
    for(unsigned int i=0; i < V; i++)
    {
        char trash;
        double x, y, z;
        ifile >> trash >> x >> y >> z;
        A(0,i) = x;
        A(1,i) = y;
        A(2,i) = z;    
    }

    //Cell2Ds properties (faces)
    //Fill the matrix cell2Ds_vertices <--> matrix B 
    vector<vector<unsigned int>> &B = P.cell2Ds_vertices;
    B.reserve(F);
    for(unsigned int i=0; i < F; i++)
    {
        char trash;
        ifile >> trash;
        vector<unsigned int> face;
        face.reserve(p);
        for(size_t j=0; j < p; j++)
        {
            unsigned int u;
            ifile >> u;
            face.push_back(u);
        }
        B.push_back(face);
    }

    //Call the function to finish to fill the polyhedron struct  
    finish_to_fill_struct(P);

    return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
void Export_polyhedron(Polyhedron &P)
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
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
void Visualize_polyhedron(Polyhedron &P)
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


    Eigen::MatrixXd &A = P.cell0Ds_coordinates;
    Eigen::MatrixXi &B = P.cell1Ds_extrema;
    vector<vector<unsigned int>> &C = P.cell2Ds_vertices;
    vector<vector<unsigned int>> &D = P.cell2Ds_edges;

    cout<<"Cell0Ds_coordinates: "<<endl;
    cout<<A<<endl;
    cout<<endl;

    cout<<"Cell1Ds_extrema: "<<endl;
    cout<<B<<endl;
    cout<<endl;
    
    cout<<"Cell2Ds_vertices: "<<endl;
    for(unsigned int i=0; i < C.size(); i++){
        cout<<'f'<<i<<": ";
        vector<unsigned int> &face = C[i];
        for(unsigned int j=0; j < face.size(); j++){
            cout<<face[j]<<' ';
        }
        cout<<endl;
    }
    cout<<endl;

    cout<<"Cell2Ds_edges: "<<endl;
    for(unsigned int i=0; i < D.size(); i++){
        cout<<'f'<<i<<": ";
        vector<unsigned int> &face = D[i];
        for(unsigned int j=0; j < face.size(); j++){
            cout<<face[j]<<' ';
        }
        cout<<endl;
    }
    cout<<endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
pair<vector<Eigen::Vector3d>, vector<Eigen::Vector3i>> Triangulation_basic_step(const Eigen::Vector3d &A, const Eigen::Vector3d &B, const Eigen::Vector3d &C, const unsigned int b)
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
void ClassI_polyhedron(Polyhedron &P, const unsigned int b, unsigned int p, unsigned int q)
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
    if(p == 3 && q == 3)
    {
        V = 2*T + 2;
        E = 6*T;
        F = 4*T;
    }
    else if((p == 3 &&  q == 4) || (p == 4 && q == 3))
    {
        V = 4*T + 2;
        E = 12*T;
        F = 8*T;

        q = 4;
    }
    else if((p == 3 && q == 5) || (p == 5 && q == 3))
    {
        V = 10*T + 2;
        E = 30*T;
        F = 20*T;

        q = 5;
    }

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
};


////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<unsigned int> cycled_face_for_dual(vector<unsigned int>& face_new, MatrixXd& coord) {
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





}