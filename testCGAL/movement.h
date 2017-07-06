#include <algorithm> 
#include <limits>

double squaredHalfAnnulusWidth(Point &v1, Point &v2, Point &v3, Point &v4, Point &v5){
	double sa = pow(v1.x(),2) + pow(v1.y(),2) + pow(v1.z(),2);
	double sb = pow(v2.x(),2) + pow(v2.y(),2) + pow(v2.z(),2);
	double sc = pow(v3.x(),2) + pow(v3.y(),2) + pow(v3.z(),2);
	double sd = pow(v4.x(),2) + pow(v4.y(),2) + pow(v4.z(),2);
	double se = pow(v5.x(),2) + pow(v5.y(),2) + pow(v5.z(),2);

	double a1 = 2*( v1.x() - v5.x() );
	double b1 = 2*( v1.y() - v5.y() );
	double c1 = 2*( v1.z() - v5.z() );
	double d1 = sa - se;

	double a2 = 2*( v2.x() - v4.x() );
	double b2 = 2*( v2.y() - v4.y() );
	double c2 = 2*( v2.z() - v4.z() );
	double d2 = sb - sd;
	
	double a3 = 2*( v4.x() - v3.x() );
	double b3 = 2*( v4.y() - v3.y() );
	double c3 = 2*( v4.z() - v3.z() );
	double d3 = sd - sc;

	double m01 = a1*b2 - a2*b1;
	double m02 = a1*b3 - a3*b1;
	double m12 = a2*b3 - a3*b2;
	double x01 = b1*c2 - b2*c1;
	double x02 = b1*c3 - b3*c1;
	double x12 = b2*c3 - b3*c2;
	double y01 = c1*a2 - c2*a1;
	double y02 = c1*a3 - c3*a1;
	double y12 = c2*a3 - c3*a2;
	double z01 = a1*b2 - a2*b1;
	double z02 = a1*b3 - a3*b1;
	double z12 = a2*b3 - a3*b2;
	double g = m01*c3 - m02*c2 + m12*c1;
	double x012 = x01*d3 - x02*d2 + x12*d1;
	double y012 = y01*d3 - y02*d2 + y12*d1;
	double z012 = z01*d3 - z02*d2 + z12*d1;
	double den = 1/g;

	double xannulus = x012*den;
	double yannulus = y012*den;
	double zannulus = z012*den;

	double r2 = pow( v2.x() - xannulus,2 ) + pow( v2.y() - yannulus,2 ) + pow( v2.z() - zannulus,2 );
	double R2 = pow( v1.x() - xannulus,2 ) + pow( v1.y() - yannulus,2 ) + pow( v1.z() - zannulus,2 );

	double w = ( r2 + R2 - 2*sqrt(r2*R2) )/4;
	return w;
}

int getIndex(Delaunay &dt,Cell &c1,Cell &c2){
	if( c2->neighbor(0) == c1)
		return 0;
	else if( c2->neighbor(1) == c1)
		return 1;
	else if( c2->neighbor(2) == c1)
		return 2;
	else 
		return 3;
}

void getTolerance(Delaunay &dt, Vertex v){
	list<Facet> facets;
	dt.finite_incident_facets(v,back_inserter(facets));
	for(std::list<Facet>::iterator it = facets.begin(); it != facets.end(); it++){ 	
		Cell c1 = it->first;
		Cell c2 = c1->neighbor(it->second);
		if( dt.is_infinite(c1) || dt.is_infinite(c2)) continue;
		// vertex of facet
		Vertex v1 = c1->vertex( it->second );
		Vertex v2 = c1->vertex( (it->second+1)%4 );
		Vertex v3 = c1->vertex( (it->second+2)%4 );
		Vertex v4 = c1->vertex( (it->second+3)%4 );		
		Vertex v5 = c2->vertex( getIndex( dt,c1,c2 ) );

		double tol = squaredHalfAnnulusWidth( v1->point(), v2->point(), v3->point(), v4->point(), v5->point() );
		if(tol < v->info().m_tolerance)
			v->info().m_tolerance = tol;
	}

	list<Cell> cells;
	dt.finite_incident_cells(v,back_inserter(cells));
	for(std::list<Cell>::iterator it = cells.begin(); it != cells.end(); it++){ 
		Cell c1 = *it;
		Cell c2 = c1->neighbor( c1->index(v) );
		if(dt.is_infinite(c2)) continue;

		int indexv = c1->index(v);
		Vertex v1 = v;
		Vertex v2 = c1->vertex( (indexv+1)%4 );
		Vertex v3 = c1->vertex( (indexv+2)%4 );
		Vertex v4 = c1->vertex( (indexv+3)%4 );		
		Vertex v5 = c2->vertex( getIndex( dt,c1,c2 ) );
		double tol = squaredHalfAnnulusWidth( v1->point(), v2->point(), v3->point(), v4->point(), v5->point() );
		if(tol < v->info().m_tolerance )
			v->info().m_tolerance = tol;
	}
}

void sphere_properties(Point &v1, Point &v2, Point &v3, Point &v4, Point &v5, Point &c, double &r){
	double sa = pow(v1.x(),2) + pow(v1.y(),2) + pow(v1.z(),2);
	double sb = pow(v2.x(),2) + pow(v2.y(),2) + pow(v2.z(),2);
	double sc = pow(v3.x(),2) + pow(v3.y(),2) + pow(v3.z(),2);
	double sd = pow(v4.x(),2) + pow(v4.y(),2) + pow(v4.z(),2);
	double se = pow(v5.x(),2) + pow(v5.y(),2) + pow(v5.z(),2);

	double a1 = 2*( v1.x() - v5.x() );
	double b1 = 2*( v1.y() - v5.y() );
	double c1 = 2*( v1.z() - v5.z() );
	double d1 = sa - se;

	double a2 = 2*( v2.x() - v4.x() );
	double b2 = 2*( v2.y() - v4.y() );
	double c2 = 2*( v2.z() - v4.z() );
	double d2 = sb - sd;
	
	double a3 = 2*( v4.x() - v3.x() );
	double b3 = 2*( v4.y() - v3.y() );
	double c3 = 2*( v4.z() - v3.z() );
	double d3 = sd - sc;

	double m01 = a1*b2 - a2*b1;
	double m02 = a1*b3 - a3*b1;
	double m12 = a2*b3 - a3*b2;
	double x01 = b1*c2 - b2*c1;
	double x02 = b1*c3 - b3*c1;
	double x12 = b2*c3 - b3*c2;
	double y01 = c1*a2 - c2*a1;
	double y02 = c1*a3 - c3*a1;
	double y12 = c2*a3 - c3*a2;
	double z01 = a1*b2 - a2*b1;
	double z02 = a1*b3 - a3*b1;
	double z12 = a2*b3 - a3*b2;
	double g = m01*c3 - m02*c2 + m12*c1;
	double x012 = x01*d3 - x02*d2 + x12*d1;
	double y012 = y01*d3 - y02*d2 + y12*d1;
	double z012 = z01*d3 - z02*d2 + z12*d1;
	double den = 1/g;

	double xannulus = x012*den;
	double yannulus = y012*den;
	double zannulus = z012*den;

	c = Point(xannulus,yannulus,zannulus);

	double r2 = pow( v2.x() - xannulus,2 ) + pow( v2.y() - yannulus,2 ) + pow( v2.z() - zannulus,2 );
	double R2 = pow( v1.x() - xannulus,2 ) + pow( v1.y() - yannulus,2 ) + pow( v1.z() - zannulus,2 );
	double w = ( r2 + R2 - 2*sqrt(r2*R2) )/4;
	//cout<<"r: "<<sqrt(r2)<<endl;
	r = sqrt(r2)+sqrt(w);
}

double inline dot(Point &a,  Point &b){
	return  a.x()*b.x() + a.y()*b.y() +  a.z()*b.z();
}

void boundary_sphere_properties(Point &v1, Point &v2, Point &v3, Point &v4, Point &c, double &r){
	/*Point a = Point(v2.x() - v3.x(), v2.y() - v3.y(),v2.z() - v3.z());
	Point b = Point(v4.x() - v3.x(), v4.y() - v3.y(),v4.z() - v3.z());


	Point n = Point( a.y()*b.z() - a.z()*b.y(), a.z()*b.x() - a.x()*b.z(),  a.x()*b.y() - a.y()*b.x() );
	double normal = sqrt( n.x()*n.x() + n.y()*n.y() + n.z()*n.z() );
	n = Point(n.x()/normal,n.y()/normal,n.z()/normal);
	Point v = Point(v1.x() - v3.x(), v1.y() - v3.y(),v1.z() - v3.z());*/

	Vector va = v2 - v3;
	Vector vb = v4 - v3;
	Vector nr = CGAL::cross_product(va,vb);
	nr = unit_vector(nr);
	Vector vv = v1 - v3;
	c = CGAL::centroid(v2,v3,v4);
	r = (nr*vv)/2;
}

int getNearPoint(Point &p0,  Point &p1 , Point &pc, double &r , double &t){
	//cout<<"stop getNearPoint -1"<<endl;
	Vector d = p1 - p0;
	Vector f = p0 - pc;
	double a = d*d;
	double b = 2*f*d;
	double c = f*f- r*r;
	double discriminant = b*b-4*a*c;
	if( discriminant < 0 )
		return 0;
	else{
		//cout<<"stop getNearPoint 0"<<endl;
		discriminant = sqrt( discriminant );
		t = (-b - discriminant)/(2*a);
		//cout<<"stop getNearPoint 1"<<endl;
		if( t >= 0 && t <= 1 ){
			return 1;
		}
		return 0;
	}
}

int getDisplacedPoint(Delaunay &dt, Vertex &v, Point &p1){
	double t = 1000;
	list<Facet> facets;
	dt.finite_incident_facets(v,back_inserter(facets));
	for(std::list<Facet>::iterator it = facets.begin(); it != facets.end(); it++){ 	
		Cell c1 = it->first;
		Cell c2 = c1->neighbor(it->second);
		///if( dt.is_infinite(c1) || dt.is_infinite(c2)) continue;
		Vertex v1 = c1->vertex( it->second );
		Vertex v2 = c1->vertex( (it->second+1)%4 );
		Vertex v3 = c1->vertex( (it->second+2)%4 );
		Vertex v4 = c1->vertex( (it->second+3)%4 );		
		//Vertex v5 = c2->vertex( getIndex( dt,c1,c2 ) );
		Vertex v5 = dt.mirror_vertex(c1,it->second);

		double ta;Point center;double radio;
		if( dt.is_infinite(c1) )
			boundary_sphere_properties(v5->point(),v2->point(),v3->point(),v4->point(),center,radio);
		else if( dt.is_infinite(c2) )
			boundary_sphere_properties(v1->point(),v2->point(),v3->point(),v4->point(),center,radio);
		else
			sphere_properties(v1->point(), v2->point(), v3->point(), v4->point(), v5->point(), center, radio);

		int step = getNearPoint(v->point(),p1,center,radio,ta);
		if(step){
			/*cout<<center<<" "<<radio<<endl;
			cout<<ta<<endl;*/
			if(ta < t)  t = ta;
		}
	}

	list<Cell> cells;
	dt.finite_incident_cells(v,back_inserter(cells));
	for(std::list<Cell>::iterator it = cells.begin(); it != cells.end(); it++){ 
		Cell c1 = *it;
		Cell c2 = c1->neighbor( c1->index(v) );
		//if(dt.is_infinite(c2)) continue;

		int indexv = c1->index(v);
		Vertex v1 = v;
		Vertex v2 = c1->vertex( (indexv+1)%4 );
		Vertex v3 = c1->vertex( (indexv+2)%4 );
		Vertex v4 = c1->vertex( (indexv+3)%4 );		
		//Vertex v5 = c2->vertex( getIndex( dt,c1,c2 ) );
		Vertex v5 = dt.mirror_vertex(c1,indexv);
		double ta;Point center;double radio;

		if( dt.is_infinite(c2) )
			boundary_sphere_properties(v1->point(),v2->point(),v3->point(),v4->point(),center,radio);
		else
			sphere_properties( v1->point(), v2->point(), v3->point(), v4->point(), v5->point(), center, radio);

		int step = getNearPoint(v->point(),p1,center,radio,ta);
		if(step){
			//cout<<center<<" "<<radio<<endl;
			//cout<<ta<<endl;
			if(ta < t)  t = ta;
		}
	}
	//cout<<t<<endl;
	if(t > 0){		
		//t = 0.7*t;
		p1 = Point( v->point().x() + t*(p1.x() - v->point().x()) , v->point().y() + t*(p1.y() - v->point().y()), v->point().z() +t*(p1.z() - v->point().z()) );
		return 1;
	}
	else{
		p1 = v->point();
		return 0;
	}
}