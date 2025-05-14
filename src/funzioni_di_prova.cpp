// prova funzione di proiezione del solido su una sfera centrata nell'origine e di raggio unitario


 void project_points_onto_sphere(Polyhedron &P)
	 {
	// salvo il numero di vertici
	MatrixXd &S = P.cell0Ds_coordinates;
	double &n = P.numcell0ds; 
	//calcolo il baricentro del solido
	barycenter = VectorXi: Zero(3);
	for (i =0; i<n; i++)
	{{for (j =0; j<3; j++)
	barycenter(j) += S(i,j)/n;}
	 }
	 // traslo rispetto al baricentro
	 for(size_t h = 0; h<n;h++)
	 {{for (j =0; j<3; j++)
	 S(h,j) -=barycenter(j);}};
	
	// proietto
	for (i = 0; i<n; i++)
	{double norm = sqrt(S(i,0)^2+S(i,1)^2+S(i,2)^2);
		for (j=0;j<3;j++)
		{S(i,j) /= abs(norm)}}	;
	};
 
