
double sum_norm_square(Delaunay &dt,Facet &f){
	Vertex v0 = f.first->vertex( (f.second) );
	double sum=0;
	for(int i = 1; i < 4 ; ++i){
		//sum+=normSquare(rest( f.first->vertex( (f.second+i)%4 )->point() , v0->point() ));
		Vector vec = f.first->vertex( (f.second+i)%4 )->point() - v0->point();
		sum += vec*vec;
	}
	return sum;
}

double volume_around_vertex(Delaunay &dt,Vertex &v){
	list<Cell> incidentCells;
	dt.finite_incident_cells(v,back_inserter(incidentCells));
	double sum=0;
	for(std::list<Cell>::iterator it = incidentCells.begin(); it != incidentCells.end(); it++)
		sum+=abs( dt.tetrahedron(*it).volume() ) ;
	return sum;
}

Point inner_smoothed_point(Delaunay &dt, Vertex &v){
	list<Cell> incidentCells;
	dt.finite_incident_cells(v,back_inserter(incidentCells));
	Vector total=Vector(0,0,0);
	//cout<<"type : "<<v->info().m_type<<endl;

	
	/*cout<<"old point: "<<v->point()<<endl;
	cout<<"new point: "<<v->info().m_smoothed_point<<endl;

	vector<Point> points;
	points.push_back(v->point());
	points.push_back(v->info().m_smoothed_point);
	drawMesh(dt,points);
	getchar();*/

	for(std::list<Cell>::iterator it = incidentCells.begin(); it != incidentCells.end(); it++){ 
		Cell c = *it;
		Facet f = Facet(c,c->index(v));
		//cout<<"dt.triangle(f): "<<dt.triangle(f)<<endl;
		double sumNS= sum_norm_square(dt,f);
		//cout<<"sumNS: "<<sumNS<<endl;
		double S = sqrt( dt.triangle(f).squared_area() );
		//cout<<"S: "<<S<<endl;
		Vector n = unit_normal_facet(dt,f);
		//cout<<"n: "<<n<<endl;
		total = total + sumNS * S *n;
	}
	double ohmio = volume_around_vertex(dt,v);
	//cout<<"ohmio: "<<ohmio<<endl;
	Point x = v->point() - total/(6*ohmio);
	return x;
}

bool inner_smooth(Delaunay &dt, Vertex &v){
	Point p = inner_smoothed_point(dt,v);
	/*cout<<"new point: "<<p<<endl;
	cout<<"old point: "<<v->point()<<endl;

	vector<Point> points;
	points.push_back(p);
	points.push_back(v->point());
	drawMesh(dt,points);
	getchar();*/

	if( is_embedded_after_displacement(dt,v,p) ){
		v->set_point(p);
		cout<<"setted point"<<endl;
		return true;
	}
	cout<<"not setted point"<<endl;
	return false;
}

Point displacement(Delaunay &dt, Vertex &v){
	double t = 1, t_aux;	
	list<Cell> incidentCells;
	dt.finite_incident_cells(v,back_inserter(incidentCells));
	Vector A = v->point() - CGAL::ORIGIN;
	Vector B = v->info().m_smoothed_point - CGAL::ORIGIN;
	Vector X,Y,n;
	for(std::list<Cell>::iterator it = incidentCells.begin(); it != incidentCells.end(); it++){ 
		Cell c = *it;
		int i = c->index(v);
		X = c->vertex( dt.vertex_triple_index( i , 1 ) )->point() - c->vertex( dt.vertex_triple_index( i , 0 ) )->point();
		Y = c->vertex( dt.vertex_triple_index( i , 2 ) )->point() - c->vertex( dt.vertex_triple_index( i , 0 ) )->point();
		n = CGAL::cross_product(X,Y);
		if( n*(B-A) != 0){ // not parallel
			double d = n*Y;
			t_aux = ( d - n*A )/ ( n*(B-A) );
			if( t_aux >= 0 && t_aux <= 1){// intersection
				if( t_aux < t)
					t = t_aux;
			}
		}
	}
	Vector intersection;
	if(t == 1)
		intersection = A + t*(B - A);
	else
		intersection = A + t*0.98*(B - A);

	return CGAL::ORIGIN + intersection;
}

bool boundary_smooth( Delaunay &dt, Vertex &v ){
	/*cout<<"point: "<<v->point()<<endl;
	cout<<"smooth point: "<<v->info().m_smoothed_point <<endl;
	cout<<"regions: "<<v->info().m_regions_around.size() <<endl;*/
	if(v->info().m_smoothed_point == v->point() ) return false;
	Point p = displacement(dt,v);


	/*cout<<"new point: "<<p<<endl;
	cout<<"old point: "<<v->point()<<endl;

	vector<Point> points;
	points.push_back(v->point());
	points.push_back(p);	
	drawMesh(dt,points);
	getchar();*/

	if( is_embedded_after_displacement(dt,v,p) ){
		v->set_point(p);
		cout<<"setted boundary point "<<endl;
		return true;
	}
	cout<<"not setted boundary point"<<endl;
	return false;
}

void ODT(Delaunay &dt, int num){
	mark_boundary(dt);
	vector<Vertex> boundary;
	boundary_vertices(dt,boundary);
	for( int i = 0 ; i < num ; ++i){
		for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) {

			Vertex v = vi;
			//cout<<"v point "<<v->point()<<endl;
			if( v->info().is_corner() ){
				//cout<<"corner "<<endl;
			}
			else if( v->info().is_boundary() ){
				//cout<<"boundary "<<endl;
				boundary_smooth(dt,v);
			}		
			else{//inner
				//cout<<"internal "<<endl;
				inner_smooth(dt,v);
			}
		}
	}
}