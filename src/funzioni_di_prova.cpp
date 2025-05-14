// prova funzione di proiezione del solido su una sfera centrata nell'origine e di raggio unitario

#include "Polyhedron.hpp"
#include <Eigen/Eigen>

namespace PolyhedronLibrary{
 void project_points_onto_sphere(Polyhedron &P)
	 {
	// salvo il numero di vertici
	MatrixXd &S = P.cell0Ds_coordinates;
	int n = P.num_cell0Ds; 
	//calcolo il baricentro del solido
	VectorXd barycenter =  VectorXd::Zero(3);
	for (int i =0; i<n; i++)
	{for (int j =0; j<3; j++)
		barycenter(j) += S(j,i)/n;
	 }
	 // traslo rispetto al baricentro
	 for(int h = 0; h<n;h++)
	 {
		 for (int j =0; j<3; j++)
		 {S(j,h) -=barycenter(j);}
 
	}
	
	// proietto
	for (int i = 0; i<n; i++)
	{double norm = sqrt(S(0,i)*S(0,i)+S(1,i)*S(1,i)+S(2,i)*S(2,i));
		for (int j=0;j<3;j++)
		{S(j,i) /= abs(norm);}
	}	;
	};
}
