// prova funzione di proiezione del solido su una sfera centrata nell'origine e di raggio unitario

#include "Polyhedron.hpp"
#include <Eigen/Dense>

namespace PolyhedronLibrary{
 void project_points_onto_sphere(Polyhedron &P)
	 {
	// salvo il numero di vertici
	MatrixXd &S = P.cell0Ds_coordinates;
	int n = P.num_cell0Ds; 
	//calcolo il baricentro del solido
	VectorXd barycenter =  VectorXd::Zero(3);
	for (int i =0; i<n; i++)
	{{for (int j =0; j<3; j++)
	barycenter(j) += S(i,j)/n;}
	 }
	 // traslo rispetto al baricentro
	 for(size_t h = 0; h<n;h++)
	 {
		 for (int j =0; j<3; j++)
		 {S(h,j) -=barycenter(j);}
 
	}
	
	// proietto
	for (int i = 0; i<n; i++)
	{double norm = sqrt(S(i,0)*S(i,0)+S(i,1)*S(i,1)+S(i,2)*S(i,2));
		for (int j=0;j<3;j++)
		{S(i,j) /= abs(norm);}
	}	;
	};
}
