
void cell_dihedral_angles(Cell &it){
	vector<double> dangles(6);
	Point v0 = it->vertex(0)->point();
	Point v1 = it->vertex(1)->point();
	Point v2 = it->vertex(2)->point();
	Point v3 = it->vertex(3)->point();
	dangles[0] = dihedral_angle(v0,v1,v2,v3);
	dangles[1] = dihedral_angle(v1,v2,v0,v3);
	dangles[2] = dihedral_angle(v0,v2,v1,v3);
	dangles[3] = dihedral_angle(v0,v3,v2,v1);
	dangles[4] = dihedral_angle(v3,v1,v2,v0);
	dangles[5] = dihedral_angle(v2,v3,v0,v1);
	it->info().m_dangles = dangles;
}

void quality_histogram(Delaunay & dt, const char* filename){
	double avgAngleTriangle=0;
	int r=5,index,nro=0;
	double sum=0;
	int sizeVector=ceil(180/r);
	vector<int> hist(sizeVector, 0);
	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		Cell c = it;
		cell_dihedral_angles(c);
		for(unsigned int i = 0; i < c->info().m_dangles.size() ; ++i){
			index=ceil(c->info().m_dangles[i]/r)-1;
			hist[index]++;
			sum+=c->info().m_dangles[i];
			++nro;
		}		
	}
	avgAngleTriangle=sum/nro;
	ofstream o(filename);
	for(int i=0;i<sizeVector;i++){
		o<<hist[i]<<endl;
	}
	o<<avgAngleTriangle<<endl;
}

void mesh_boundary_histogram(Delaunay & dt, const char* filename){
	vector<Vertex> boundaryVertices;
	boundary_vertices(dt,boundaryVertices);
	Vertex v;

	double avgAngleBoundary=0;
	int r=40,index;
	double sum=0;
	int sizeVector=ceil(360/r);
	vector<int> hist(sizeVector, 0);
	int nro=0;

	for(unsigned int i = 0; i < boundaryVertices.size(); ++i){
		v = boundaryVertices[i];
		//if( v->info().m_boundary || v->info().has_more_than_two_regions() || v->info().is_singular() )		continue;
		if( v->info().m_boundary || v->info().is_singular() )		continue;

		if( v->info().m_regions_around.size() == 2){
			solid_angle_mesh_vertex(dt,v);
			if(v->info().m_s_lower_angle <= 0 || !IsNumber(v->info().m_s_lower_angle) ) continue;
			index=ceil(v->info().m_s_lower_angle/r)-1;
			hist[index]++;
			sum+=v->info().m_s_lower_angle;
			++nro;
		}
		else{
			vector<double> solidAngles;
			solid_angle_mesh_vertex2(dt,v,solidAngles);
			for( unsigned int j = 0; j < solidAngles.size(); ++j ){
				double angle = solidAngles[j];
				//cout<<angle<<endl;
				if( angle <= 0 || angle > 360  || !IsNumber(angle) ) continue;
				index=ceil( abs(angle)/r )-1;
				hist[index]++;
				sum += angle;
				++nro;	
			}
		}
	}
	avgAngleBoundary=sum/(nro);
	ofstream o(filename);
	for(int i=0;i<sizeVector;i++){
		o<<hist[i]<<endl;
	}
	o<<avgAngleBoundary<<endl;

	restart_vertex_state(dt);
}

void polyhedron_boundary_histogram(Delaunay & dt, const char* filename){
	vector<Vertex> boundaryVertices;
	boundary_vertices(dt,boundaryVertices);
	Vertex v;

	//initiate_process(dt,boundaryVertices);

	double avgAngleBoundary=0;
	int r=40,index;
	double sum=0;
	int sizeVector=ceil(360/r);
	vector<int> hist(sizeVector, 0);
	int nro=0;
	for(unsigned int i = 0; i < boundaryVertices.size(); ++i){
		v = boundaryVertices[i];
		//if( v->info().m_boundary || v->info().has_more_than_two_regions() || v->info().is_singular() )		continue;
		if( v->info().m_boundary || v->info().is_singular() )		continue;
		if( v->info().m_regions_around.size() == 2){
			solid_angle_polyhedron_vertex(v);
			if(v->info().m_s_lower_angle <= 0 || v->info().m_s_lower_angle >= 360  || !IsNumber(v->info().m_s_lower_angle) ) continue;
			index=ceil(abs(v->info().m_s_lower_angle)/r)-1;
			hist[index]++;
			sum+=v->info().m_s_lower_angle;
			++nro;
		}
		else{
			//cout<<"size node: "<< v->info().m_regions_around.size()<<endl;
			vector<double> solidAngles;
			solid_angle_polyhedron_vertex2(v,solidAngles);
			for( unsigned int j = 0; j < solidAngles.size(); ++j ){
				double angle = solidAngles[j];
				//cout<<angle<<endl;
				if( angle <= 0 || angle > 360  || !IsNumber(angle) ) continue;
				index=ceil( abs(angle)/r )-1;
				hist[index]++;
				sum += angle;
				++nro;				
			}
			//cout<<"end "<<endl;
		}
	}
	avgAngleBoundary=sum/(nro);
	ofstream o(filename);
	for(int i=0;i<sizeVector;++i){
		o<<hist[i]<<endl;
	}
	o<<avgAngleBoundary<<endl;

	restart_vertex_state(dt);
}

void mesh_curvature_histogram(Delaunay & dt, const char* filename){
	vector<Vertex> boundaryVertices;
	boundary_vertices(dt,boundaryVertices);
	Vertex v;

	double avgAngleBoundary=0;
	double r=0.01,index;
	double sum=0;

	double limit = 0.5;
	int sizeVector=ceil(limit/r);
	vector<int> hist(sizeVector, 0);
	int nro=0;
	double maxvalue = 0;
	for(unsigned int i = 0; i < boundaryVertices.size(); ++i){
		Vertex v = boundaryVertices[i];
		if( v->info().is_singular_2() || v->info().is_corner() ) continue;
		set_stars_vertex(dt,v);
		vector<double> curvatures;
		vertex_curvature(v,curvatures);
		for( unsigned int j = 0; j < curvatures.size(); ++j){
			double k = curvatures[j];
			if( k > maxvalue) maxvalue=k;
			if( k < 0 || k > limit || !IsNumber(k) ) continue;
			if(k == 0)	k = 0.0001;		
			index=ceil( abs(k)/r )-1;
			hist[index]++;
			sum += k;
			++nro;	
		}
	}
	avgAngleBoundary=sum/(nro);
	ofstream o(filename);
	for(int i=0;i<sizeVector;i++){
		o<<hist[i]<<endl;
	}
	o<<avgAngleBoundary<<endl;
	o<<sum<<endl;
	o<<maxvalue<<endl;
	o.close();
	restart_vertex_state(dt);
}