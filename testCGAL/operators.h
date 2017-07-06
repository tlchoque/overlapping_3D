
bool is_obtuse(Vector a , Vector b){
	return a.x()*b.y() - a.y()*b.x() < 0;
}

double angle(Point &i, Point &v, Point &j){
	Vector a = i - v;
	Vector b = j - v;
	double th = acos( (a*b)/( norm(a)*norm(b) ) );

	return th;
}

double angle2(Vector x, Vector y){
	double th = acos( ( x*y )/ ( norm( x ) * norm( y ) ) );
	if( is_obtuse(x,y) ){
		//cout<<"signed"<<endl;
		th = -th;
	}
	return th;
}


double cotan(Vertex v0, Vertex v1, Vertex v2){
	Vector u = v1->point() - v0->point();
	Vector v = v2->point() - v0->point();
	double uv = u*v;
	double denom = sqrt( ( u*u  )*( v*v ) - uv*uv );
	if( denom == 0 ) return 0;
	return uv/denom;		
}

double length(Point &a, Point &b){
	Vector d = a - b;
	return norm(d);
}

double cot(double a){
	return cos(a)/sin(a);
}


double voronoi_area(Point &i, Point &v, Point &j){
	return (  CGAL::squared_distance(v,i)*cot( angle(v,i,j) ) + CGAL::squared_distance(v,j)*cot( angle(v,i,j) ) )/ 8;
}

double triangle_area(Point &v1, Point &v2, Point &v3){
	Vector a = v2 - v1;
	Vector b = v3 - v1;
	return norm( CGAL::cross_product(a,b) )/2;
}

//for optimal delaunay
Vector unit_normal_triangle(Point &v1, Point &v2, Point &v3){
	Vector a = v2 - v1;
	Vector b = v3 - v1;
	Vector cp = CGAL::cross_product( a , b );
	cp = unit_vector(cp);
	return cp;
}

Vector norm_tangent_plane(vector<Vertex> &vertices, Vertex &v){
	Vector sum = Vector(0,0,0);
	for(unsigned int i = 0; i < vertices.size(); ++i){
		int next_i = (i+1)%vertices.size();
		double S = triangle_area( v->point(),vertices[i]->point(),vertices[next_i]->point() );
		Vector n = unit_normal_triangle( v->point(),vertices[i]->point(),vertices[next_i]->point() );
		sum = sum + S*n;
	}
	sum = unit_vector(sum);
	return sum;
}



double mixed_area(vector<Vertex> &vertices, Vertex &v){
	double mixed_area = 0;
	for(unsigned int i = 0; i < vertices.size(); ++i){
		int next_i = ( i + 1 )%vertices.size();
		if( CGAL::angle(vertices[i]->point(), v->point(), vertices[next_i]->point()) == CGAL::OBTUSE )
			mixed_area += triangle_area(vertices[i]->point(), v->point(), vertices[next_i]->point() ) / 2;
		else if( CGAL::angle(v->point(), vertices[i]->point(),  vertices[next_i]->point()) == CGAL::OBTUSE  || 
				 CGAL::angle(v->point(),  vertices[next_i]->point() , vertices[i]->point()) == CGAL::OBTUSE	)
			mixed_area += triangle_area(vertices[i]->point(), v->point(), vertices[next_i]->point() ) / 4;
		else
			mixed_area += voronoi_area( vertices[i]->point(), v->point(), vertices[next_i]->point() );
	}
	return mixed_area;
}

double mean_curvature(vector<Vertex> &vertices, Vertex &v){
	Vector k = Vector(0,0,0);
	for(unsigned int i = 0; i < vertices.size(); ++i){
		int next_i = ( i + 1 )%vertices.size();
		size_t prev_i = ( i + vertices.size() - 1 )%vertices.size();
		Vector cen = v->point() - vertices[i]->point();

		Vector prev_t = vertices[prev_i]->point() - v->point();
		Vector t = vertices[i]->point() - vertices[prev_i]->point() ;
		Vector next_t = vertices[next_i]->point() - vertices[i]->point();
		Vector aux = v->point() - vertices[next_i]->point();

		//cen = ( cot( angle2( prev_t,t ) ) + cot( angle2( next_t,aux ) ) )*cen;
		cen = ( cotan( vertices[prev_i],v,vertices[i] ) + cotan( vertices[next_i],v,vertices[i]  ) )*cen;

		k = k + cen;
	}

	//vector<Point> pts;
	//pts.push_back(v->point() );
	//pts.push_back( CGAL::ORIGIN + (200*k + ( v->point() - CGAL::ORIGIN ) ) ); // for grey

	//drawMesh(dt,pts);
	//getchar();

	//cout<<"normal curvature: "<<k<<endl;
	double Am = mixed_area(vertices,v);
	k = k / (Am*2);
	return norm(k)/2;
}


double mean_curvature2(vector<Vertex> &vertices, Vertex &v){
	Vector k = Vector(0,0,0);
	Vector sum = Vector(0,0,0);
	for(unsigned int i = 0; i < vertices.size(); ++i){
		int next_i = ( i + 1 )%vertices.size();
		size_t prev_i = ( i + vertices.size() - 1 )%vertices.size();
		//Vector cen = vertices[i]->point() - v->point();
		Vector cen = v->point() - vertices[i]->point();
		cen = ( cotan( vertices[prev_i],v,vertices[i] ) + cotan( vertices[next_i],v,vertices[i]  ) )*cen;
		k = k + cen;
		//cout<<"cen and k "<<cen<<" "<<k<<endl;

		//normal vector
		double S = triangle_area( v->point(),vertices[i]->point(),vertices[next_i]->point() );
		Vector n = unit_normal_triangle( v->point(),vertices[i]->point(),vertices[next_i]->point() );
		sum = sum + S*n;
	}
	
	double Am = mixed_area(vertices,v);
	//cout<<"mixed area and k "<<Am<<" "<<k<<endl;
	if(Am == 0) return 0;

	//k = k / (Am*2);
	k = (1/(Am*2))*k;
	double mc = norm(k)/2;
	//cout<<"k sum mc "<<k<<" "<<sum<<" "<<mc<<endl;
	if( sum*k < 0 )
		return -mc;
	return mc;
}

double gaussian_curvature(vector<Vertex> &vertices, Vertex &v){
	double G=0;
	for(unsigned int i = 0; i < vertices.size(); ++i){
		int next_i = ( i + 1 )%vertices.size();
		Vector a = vertices[i]->point() - v->point() ;
		Vector b = v->point() - vertices[next_i]->point();

		//double theta = angle( vertices[i]->point(),v->point(),vertices[next_i]->point() );

		double theta = angle2(a, b);		
		G += theta;
	}
	double Am = mixed_area(vertices,v);
	cout<<G<<endl;
	G = 2*PI - G;
	G = G / (Am);
	return G;
}

void principal_curvatures(vector<Vertex> &vertices, Vertex &v, double &k1,double &k2){
	double kH = mean_curvature(vertices,v);
	double kG = gaussian_curvature(vertices,v);
	double d = kH*kH - kG;	
	if( d > 0){
		double difV = sqrt( d );
		k1 = kH + difV;
		k2 = kH - difV;
		cout<<" principal: "<<k1<<" * "<<k2<<endl;
	}	
}

//for histograms and 

double dihedral_angle(Point &a, Point &b, Point &c, Point &d){
	Vector b1 = c - a ;
	Vector b2 = b - a ;
	Vector b3 = d - a ;
	Vector c1 = CGAL::cross_product(b1,b2);
	Vector c2 = CGAL::cross_product(b2,b3);
	return acos(- ( c1*( 1/norm(c1) ) * c2*( 1/norm(c2) )  ))* 180.0 / PI;
}

Vector middle_point(vector<Vertex> &star){
	Vector m = Vector(0,0,0);
	for(unsigned int i = 0; i < star.size(); ++i)
		m = m + (star[i]->point() - CGAL::ORIGIN );
	return m/star.size();
}

Vector average_normal(vector<Vertex> &star, Vector &p){
	Vector q,x1,x2,normal(0,0,0);
	for(unsigned int i = 0; i < star.size(); ++i){
		int next_i = ( i + 1 )%star.size();
		x1 = ( star[i]->point() - CGAL::ORIGIN ) - p;
		x2 = ( star[next_i]->point() - CGAL::ORIGIN ) - p;
		q = CGAL::cross_product(x1,x2);
		normal = normal + q;
	}
	normal = unit_vector(normal);
	return normal;
}

Vector normal_edges(vector<Vertex> &star, Vector &p){
	Vector q,x1,x2,normal(0,0,0);
	for(unsigned int i = 0; i < star.size(); ++i){
		int next_i = ( i + 1 )%star.size();
		x1 = ( star[i]->point() - CGAL::ORIGIN ) - p;
		x2 = ( star[next_i]->point() - CGAL::ORIGIN ) - p;
		q = CGAL::cross_product(x1,x2);
		normal = normal + q;
	}
	return normal;
}


double solid_angle(Point &o, Point &a, Point &b, Point &c){
	Vector ac = a - o;
	Vector bc = b - o;
	Vector cc = c - o;
	double triple = abs( CGAL::cross_product(ac,bc)*cc );
	double A = norm(ac);
	double B = norm(bc);
	double C = norm(cc);
	return 2*atan( triple/abs( A*B*C + C*(ac*bc) + B*(ac*cc) + A*(bc*cc) ) ) * 180 / PI;	
}

void solid_angle_polyhedron_vertex(Vertex &v){
	int label = v->info().m_regions_around[0];
	vector<Vertex> star = v->info().m_s_stars[0];
	Vector p = middle_point(star);
	double total=0;

	/*cout<<"label: "<<label<<" size: "<<star.size()<<" star 0:"<< star[0]->point()<<" middle:"<< CGAL::ORIGIN + p<<endl;
	for(unsigned int i = 0; i < star.size(); ++i)
		cout<<star[i]->point()<<endl;*/

	for(unsigned int i = 0; i < star.size(); ++i){
		int next_i = ( i + 1 )%star.size();
		total += solid_angle( v->point() , CGAL::ORIGIN + p, star[i]->point(),star[next_i]->point());
	}

	//change this for  future
	total = abs(total);

	if(total > 360)		v->info().m_s_lower_angle = 720 - total ;
	else	v->info().m_s_lower_angle = total;
}

double solid_angle_for_star(vector<Vertex> star, Vertex &v){
	Vector p = middle_point(star);
	double total=0;
	for(unsigned int i = 0; i < star.size(); ++i){
		int next_i = ( i + 1 )%star.size();
		total += solid_angle( v->point() , CGAL::ORIGIN + p, star[i]->point(),star[next_i]->point());
	}
	//change this for  future
	total = abs(total);
	return total;
}

void solid_angle_polyhedron_vertex2(Vertex &v, vector<double> &solidAngles){
	//cout<<"stars size: "<<v->info().m_s_stars.size()<<endl;
	for( unsigned int i = 0 ; i < v->info().m_s_stars.size(); ++i){
		solidAngles.push_back( solid_angle_for_star( v->info().m_s_stars[i], v ) );
	}
}

void solid_angle_mesh_vertex(Delaunay &dt,Vertex v){
	list<Cell> cells;
	dt.finite_incident_cells(v,back_inserter(cells));
	int label = v->info().m_regions_around[0];
	Vertex v1,v2,v3;
	double total = 0;

	for(std::list<Cell>::iterator  it = cells.begin(); it != cells.end(); it++){ 	
		if( (*it)->info().m_label != label ) continue;
		Cell cell = *it;
		int indx = cell->index(v);	
		v1 = cell->vertex( dt.vertex_triple_index(indx,0) );
		v2 = cell->vertex( dt.vertex_triple_index(indx,1) );
		v3 = cell->vertex( dt.vertex_triple_index(indx,2) );
		total += solid_angle( cell->vertex(indx)->point(),v1->point(),v2->point(),v3->point() );
	}	
	if(total > 360)		v->info().m_s_lower_angle = 720 - total ;
	else	v->info().m_s_lower_angle = total;
}

void solid_angle_mesh_vertex2(Delaunay &dt,Vertex v,  vector<double> &solidAngles){
	list<Cell> cells;
	dt.finite_incident_cells(v,back_inserter(cells));
	int label = v->info().m_regions_around[0];
	Vertex v1,v2,v3;
	double total = 0;

	for(unsigned int i = 0; i < v->info().m_regions_around.size(); ++i){
		solidAngles.push_back(0);
	}

	for(std::list<Cell>::iterator  it = cells.begin(); it != cells.end(); it++){ 	
		int label = (*it)->info().m_label;
		vector<int>::iterator ite = std::find(v->info().m_regions_around.begin(), v->info().m_regions_around.end(), label );
		if( ite != v->info().m_regions_around.end() ){// found
			int position = std::distance( v->info().m_regions_around.begin(), ite );	
			Cell cell = *it;
			int indx = cell->index(v);	
			v1 = cell->vertex( dt.vertex_triple_index(indx,0) );
			v2 = cell->vertex( dt.vertex_triple_index(indx,1) );
			v3 = cell->vertex( dt.vertex_triple_index(indx,2) );
			solidAngles[position] += solid_angle( cell->vertex(indx)->point(),v1->point(),v2->point(),v3->point() );
		}
	}			
}

//for ODT

Vector ODT_middle_point(vector<Vertex> &star){
	Vector m = Vector(0,0,0);
	for(unsigned int i = 0; i < star.size(); ++i)
		m = m + (star[i]->info().m_smoothed_point - CGAL::ORIGIN );
	return m/star.size();
}

Vector ODT_average_normal(vector<Vertex> &star, Vector &p){
	Vector q,x1,x2,normal(0,0,0);
	for(unsigned int i = 0; i < star.size(); ++i){
		int next_i = ( i + 1 )%star.size();
		x1 = ( star[i]->info().m_smoothed_point - CGAL::ORIGIN ) - p;
		x2 = ( star[next_i]->info().m_smoothed_point - CGAL::ORIGIN ) - p;
		q = CGAL::cross_product(x1,x2);
		normal = normal + q;
	}
	normal = unit_vector(normal);
	return normal;
}

//for new repairing

Vector normal_curvature(vector<Vertex> &vertices, Vertex &v){
	Vector n = norm_tangent_plane(vertices,v);
	double c = mean_curvature(vertices,v);
	return c*n;
}

double vertex_curvature(Vertex &v, vector<double> &curvatures){
	double total=0;
	for( unsigned int i = 0 ; i < v->info().m_sa_stars.size(); ++i){
		if( v->info().m_regions_around[i] == HOLE ) continue;
		curvatures.push_back( mean_curvature( v->info().m_sa_stars[i], v ) );
	}
	for(unsigned int i = 0; i < curvatures.size(); ++i){
		total += curvatures[i];
	}
	return total;
}