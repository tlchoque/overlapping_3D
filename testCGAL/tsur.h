
int map_edge(map<pair<Vertex,Vertex>,int> &hash_edges, pair<Vertex,Vertex> &p){
    std::map<pair<Vertex,Vertex>,int>::const_iterator it = hash_edges.find( p );
    if (it == hash_edges.end())
        return -1;
    return it->second;
}

void edges_of_boundary(Delaunay &dt, vector<Edge> &edges){
	Vertex v1,v2,v3;
	map<pair<Vertex,Vertex>,int> hash_edges;
	for(Finite_facets_iterator fi=dt.finite_facets_begin();fi!=dt.finite_facets_end(); fi++){
		if( fi->first->info().m_label == fi->first->neighbor(fi->second)->info().m_label )		continue;
		if( dt.is_infinite( fi->first ) || dt.is_infinite( fi->first->neighbor(fi->second) ) )		continue;

		v1 = fi->first->vertex( dt.vertex_triple_index(fi->second,0) );
		v2 = fi->first->vertex( dt.vertex_triple_index(fi->second,1) );
		v3 = fi->first->vertex( dt.vertex_triple_index(fi->second,2) );

		Cell c = fi->first ;
		int i1 = c->index(v1);
		int i2 = c->index(v2);
		int i3 = c->index(v3);

		pair<Vertex,Vertex> e1,e2,e3,reverse_e1,reverse_e2,reverse_e3;
		e1 = make_pair(v1,v2);
		e2 = make_pair(v1,v3);
		e3 = make_pair(v2,v3);
		reverse_e1 = make_pair(v2,v1);
		reverse_e2 = make_pair(v3,v1);
		reverse_e3 = make_pair(v3,v2);
		
		if( map_edge( hash_edges,e1 ) == -1 &&  map_edge( hash_edges,reverse_e1 ) == -1 ){
			hash_edges[e1] = 1;
			hash_edges[reverse_e1] = 1;
			Edge e = Edge(c,i1,i2);
			edges.push_back(e);
		}
		if( map_edge( hash_edges,e2 ) == -1 &&  map_edge( hash_edges,reverse_e2 ) == -1 ){
			hash_edges[e2] = 1;
			hash_edges[reverse_e2] = 1;
			Edge e = Edge(c,i1,i3);
			edges.push_back(e);
		}
		if( map_edge( hash_edges,e3 ) == -1 &&  map_edge( hash_edges,reverse_e3 ) == -1 ){
			hash_edges[e3] = 1;
			hash_edges[reverse_e3] = 1;
			Edge e = Edge(c,i2,i3);
			edges.push_back(e);
		}
	}
}

int common_label(Vertex &v, Vertex &w){
	int label;
	if( v->info().m_regions_around[0] != 0	)	label = v->info().m_regions_around[0];
	else	label = v->info().m_regions_around[1];

	//label = v->info().m_regions_around[0];
	if( label == w->info().m_regions_around[0] || label == w->info().m_regions_around[1] )
		return label;
	else return -1;
}

int common_label2(Vertex &v, Vertex &w){
	for(unsigned int i = 0; i < v->info().m_regions_around.size(); ++i){
		for(unsigned int j = 0; j < w->info().m_regions_around.size(); ++j){
			if( w->info().m_regions_around[j] == v->info().m_regions_around[i] ){
				return w->info().m_regions_around[j];
			}
		}	
	}
	return -1;
}

void boundary_facets_vertex(Delaunay &dt,vector<Facet> &facets,Vertex &v,int label){
	list<Cell> incident;
	Cell vc,c;
	dt.finite_incident_cells(v,back_inserter(incident));
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		if( (*it)->info().m_label != label) continue;
		vc = *it;
		int index = vc->index(v);
		for(unsigned int j = 0; j < 3;++j){
			int index2 = dt.vertex_triple_index(index,j); 
			c = vc->neighbor( index2 );
			if(c->info().m_label != label){
				Facet f = Facet(vc,index2);
				facets.push_back(f);
			}
		}
	}
}

vector<Vertex> star_vertex(Delaunay &dt,Vertex &v,int label){
	map<Vertex,Vertex> vv;
	vector<Facet> facets;
	boundary_facets_vertex(dt,facets,v,label);
	for(unsigned int i=0;i < facets.size();++i){
		Facet f = facets[i];
		Cell c = f.first;		
		int index0 = f.second;
		int indexv = c->index(v);
		int index1 = dt.next_around_edge(index0,indexv);
		int index2 = dt.next_around_edge(index1,indexv);
		vv[c->vertex(index1)] = c->vertex(index2);
		if(i == 0)
			vv[v] = c->vertex(index1);
	}
	vector<Vertex> orderedVertices;// the first is repeated at last
	Vertex it = v;

	for(unsigned int i=0;i < facets.size();++i){
		it = vv[it];
		orderedVertices.push_back(it);
	}	
	return orderedVertices;
}

double star_volume( vector<Vertex> &star, Vector &v,Vector &p){
	Vector x1,x2,x3;
	double volTotal=0;
	for(unsigned int i = 0; i < star.size(); ++i){
		int next_i = ( i + 1 )%star.size();
		volTotal += CGAL::volume( CGAL::ORIGIN + v , star[i]->point(), star[next_i]->point(), CGAL::ORIGIN + p );
	}
	return volTotal;
}

bool is_point_in_polygon(Vector &p, vector<Vector> &points){
	Vector a(0,0,0),pv,vec1,vec2;
	double n1,n2,ip,theta,total=0;
	for(unsigned int i = 0; i < points.size();++i){
		int next_i = ( i + 1 )%points.size();
		vec1 = points[i]  -  p ;
		vec2 = points[next_i] - p ;
		pv=CGAL::cross_product(vec1,vec2);
		n1=norm(vec1);
		n2=norm(vec2);
		ip = vec1*vec2;
		theta = acos( ip/(n1*n2) );
		if(i==0)	a = Vector(pv.x(),pv.y(),pv.z());
		if( a*pv >= 0 )
			total+=theta;
		else
			total-=theta;
	}
	if ( total < 1e-2 ) return 0;
	else return 1;
}

void balance_star(Vertex &v, vector<Vertex> &star){
	double V=0,V1=0,h,hf;
	Vector p,n,m,x,Ver,newVer,point;
	vector<Vector> points;
	p = middle_point(star);
	n = average_normal(star,p);
	for(unsigned int i = 0; i < star.size(); ++i){
		h = ( (star[i]->point() - CGAL::ORIGIN) - p) * n;	
		point = (star[i]->point() - CGAL::ORIGIN)  - h*n;
		points.push_back(point);
		//cout<<"h: "<<h<<endl;
	}
	if(!is_point_in_polygon(p,points))	return;

	V = star_volume(star, p , v->point() - CGAL::ORIGIN);
	V1 = star_volume(star, p , p + n);
	hf = V/V1;
	newVer = p + hf*n;
	v->set_point(CGAL::ORIGIN + newVer);
} 

int index_label_in_star(Vertex v, int label){

	/*if( v->info().m_regions_around[0] == label)
		return 0;
	else return 1;*/

	for( unsigned int i = 0; i < v->info().m_regions_around.size(); ++i ){
		if( v->info().m_regions_around[i] == label )
			return i;
	}
	return -1;
}

void set_stars_to_vertex(Delaunay &dt, Vertex &v){
	map<Vertex,Vertex> vv;
	vector<Facet> facets;
	int label = v->info().m_regions_around[0];
	boundary_facets_vertex(dt,facets,v,label);
	for(unsigned int i=0;i < facets.size();++i){
		Facet f = facets[i];
		Cell c = f.first;		
		int index0 = f.second;
		int indexv = c->index(v);
		int index1 = dt.next_around_edge(index0,indexv);
		int index2 = dt.next_around_edge(index1,indexv);
		vv[c->vertex(index1)] = c->vertex(index2);
		if(i == 0)
			vv[v] = c->vertex(index1);
	}
	vector<Vertex> orderedVertices;// the first is repeated at last
	Vertex it = v;
	for(unsigned int i=0;i < facets.size();++i){
		it = vv[it];
		orderedVertices.push_back(it);
	}
	v->info().m_s_stars[0] = orderedVertices;
	std::reverse(orderedVertices.begin(), orderedVertices.end());
	v->info().m_s_stars[1] = orderedVertices;
}

void set_stars_to_vertex2(Delaunay &dt, Vertex &v){	
	/*cout<<"stop 0 set stars "<<v->info().m_regions_around.size() - 2<<endl;
	cout<<"size regions around: "<<v->info().m_regions_around.size() <<endl;
	cout<<"size stars: "<<v->info().m_s_stars.size() <<endl;*/
	for(unsigned int i = 0; i < v->info().m_regions_around.size() - 2;++i){
		vector<Vertex> star;
		v->info().m_s_stars.push_back(star);
	}

	/*cout<<"size regions around: "<<v->info().m_regions_around.size() <<endl;
	cout<<"size stars: "<<v->info().m_s_stars.size() <<endl;
*/
	map<int,int> label_to_index;
	for(unsigned int i = 0; i < v->info().m_regions_around.size(); ++i){
		label_to_index[ v->info().m_regions_around[i] ] = i;
	}

	vector<vector<Cell>> cellsSet( v->info().m_s_stars.size() );	

	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));	
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		cellsSet[ label_to_index[ (*it)->info().m_label ]  ].push_back(*it);
	}
	//cout<<"stop 1/2 set stars"<<endl;

	for(unsigned int i = 0; i < cellsSet.size(); ++i){
		//cout<<"item i: "<<i<<" cell size: "<<cellsSet[i].size()<<endl;
		v->info().m_s_stars[i] = ring_vertices( dt,cellsSet[i],v,v->info().m_regions_around[i] );
		/*for(unsigned int j = 0; j < v->info().m_s_stars[i].size(); ++j){
			cout<<v->info().m_s_stars[i][j]->point()<<endl;
		}
		cout<<endl;*/
	}
}

//new smoothing

void greater_edge_in_vertex(Delaunay &dt, Vertex &v){
	double max=0;
	list<Edge> incident;
	dt.finite_incident_edges(v,back_inserter(incident));
	for(std::list<Edge>::iterator it = incident.begin(); it != incident.end(); it++){
		double square = dt.segment(*it).squared_length();
		if( square  > max)
			max = square;
	}
	v->info().m_s_greater_edge = sqrt(max);
}

void initiate_process(Delaunay &dt, vector<Vertex> &boundaryVertices){
	for(unsigned int i = 0 ; i < boundaryVertices.size(); ++i ){
		Vertex v = boundaryVertices[i];
		greater_edge_in_vertex(dt,v);
		//cout<<"stop 0"<<endl;
		if( v->info().has_two_regions() )
			set_stars_to_vertex(dt,v);
		else{
			if( v->info().has_more_than_two_regions() )
				set_stars_to_vertex2(dt,v);
		}
		//cout<<"stop 1"<<endl;
	}
}

bool is_embedded_after_displacement(Delaunay &dt, Vertex v, const Point &p){
	if(v->point() == p) return true;
	Point ant = v->point();
	v->set_point(p);
	std::size_t size;
	// are incident cells well-oriented
	std::vector<Cell> cells;
	cells.reserve(64);
	dt.finite_incident_cells(v, std::back_inserter(cells));
	size = cells.size();

	for(std::size_t i=0; i<size; i++) {
		Cell c = cells[i];
		if(dt.is_infinite(c)) continue;
		if(dt.orientation(c->vertex(0)->point(), c->vertex(1)->point(), c->vertex(2)->point(), c->vertex(3)->point()) != CGAL::POSITIVE){
			v->set_point(ant);
			return false;
		}
	}
	v->set_point(ant);
	return true;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void TSUR3D( Delaunay &dt, vector<Edge> &edges){
	Vertex v1,v2;
	double h,h1,h2,h1f,h2f,vol1,vol2,rlx=0;
	Vector p1,p2,n1,n2,n,m,V1,V2,newV1,newV2,pp1,pp2;	
	int counter = 0;
	bool draw=0;
	for ( unsigned int i = 0; i < edges.size(); ++i){
		v1 = edges[i].first->vertex(edges[i].second);
		v2 = edges[i].first->vertex(edges[i].third);
		 
		vector<Point> test;
		test.push_back(v1->point());
		test.push_back(v2->point());
		
		if( ( v1->info().is_corner() || v1->info().has_more_than_two_regions() || v1->info().is_singular_2() ) || 
			( v2->info().is_corner() || v2->info().has_more_than_two_regions() || v1->info().is_singular_2() ) )		continue;


		/*if( ( v1->info().is_corner()  || v1->info().is_singular_2() ) || 
			( v2->info().is_corner()  || v1->info().is_singular_2() ) )		continue;*/

		//if(v1->info().has_more_than_two_regions() || v2->info().has_more_than_two_regions() )	draw = 1;
		if(draw){
			drawMesh(dt,test);
			getchar();
		}

		int label = common_label2(v1,v2);
		
		vector<Vertex> star1 = v1->info().m_s_stars[index_label_in_star(v1,label)];
		vector<Vertex> star2 = v2->info().m_s_stars[index_label_in_star(v2,label)];

		V1 = v1->point() - CGAL::ORIGIN;
		V2 = v2->point() - CGAL::ORIGIN;		
		p1 = middle_point(star1);
		p2 = middle_point(star2);
		pp1 = p1;
		pp2 = p2;
		n1 = average_normal(star1,p1);
		n2 = average_normal(star2,p2);
		n = n1 + n2;
		n = unit_vector(n);
		m = (V1 + V2)/2;

		h1 = (V1 - m) * n;
		h2 = (V2 - m) * n;
		if(h1 == 0 || h2 == 0) continue;//indeterminate

		p1 = V1 - h1 * n;
		p2 = V2 - h2 * n;	

		vol1 = star_volume(star1, V1, p1);
		vol2 = star_volume(star2, V2, p2);

		double divisor = (vol1/h1 + vol2/h2);
		if( divisor == 0) continue;//indeterminate

		h =  2*(vol1 + vol2)/divisor;
		h1f = rlx * h1 + (1.0 - rlx) * h;
		h2f = rlx * h2 + (1.0 - rlx) * h;
		newV1 = p1 + h1f*n;
		newV2 = p2 + h2f*n;

		//if( !IsNumber( newV1.x() ) || !IsNumber( newV1.y() ) || !IsNumber( newV1.z() ) || !IsNumber( newV2.x() ) || 
		//!IsNumber( newV2.y() ) || !IsNumber( newV2.z() ) ) draw=1;

		v1->set_point(CGAL::ORIGIN + newV1);
		v2->set_point(CGAL::ORIGIN + newV2);
		
		if(draw){
			for(unsigned int i = 0; i < star1.size(); ++i)
				cout<<star1[i]->point()<<endl;
			for(unsigned int i = 0; i < star2.size(); ++i)
				cout<<star2[i]->point()<<endl;

			double test1 = star_volume(star1, V1, newV1);
			double test2 = star_volume(star2, V2, newV2);
		
			cout<<"counter: "<<counter<<endl;
			cout<<"LABEL: "<<label<<endl;
			cout<<"v1 size labels : "<<v1->info().m_regions_around.size()<<endl;
			cout<<"v2 size labels : "<<v2->info().m_regions_around.size()<<endl;
			cout<<"h1: "<<h1<<" h2: "<<h2<<endl;
			cout<<"h: "<<h<<" n: "<<n<<endl;
			cout<<"vol1: "<<vol1<<" vol2: "<<vol2<<endl;
			cout<<"test1: "<<test1<<" test2: "<<test2<<endl;
			cout<<"sum input: "<<vol1+vol2<<endl;
			cout<<"sum output: "<<test1+test2<<endl;
			cout<<"V1: "<<V1<<" V2: "<<V2<<endl;
			cout<<"newV1: "<<newV1<<" newV2: "<<newV2<<endl;
			cout<<endl;

			test.clear();
			test.push_back(CGAL::ORIGIN + newV1);
			test.push_back(CGAL::ORIGIN + newV2);
			drawMesh(dt,test);
			getchar();
			draw = 0;
		}
		++counter;
	}
}

void balance3D( Delaunay &dt, vector<Vertex> &boundaryVertices ){
	for ( unsigned int i = 0 ; i < boundaryVertices.size() ; ++i){
		Vertex v = boundaryVertices[i] ;
		if( v->info().is_corner() || v->info().has_more_than_two_regions()  || v->info().is_singular_2() )		continue;
		vector<Vertex> star = v->info().m_s_stars[0];
		balance_star(v,star);
	}
}

void restricted_TSUR3D( Delaunay &dt, vector<Edge> &edges){
	Vertex v1,v2;
	double vol1,vol2,rlx=0;
	Vector p1,p2,n1,n2,n,m,V1,V2,newV1,newV2,pp1,pp2;	
	int counter = 0;
	for ( unsigned int i = 0; i < edges.size(); ++i){
		v1 = edges[i].first->vertex(edges[i].second);
		v2 = edges[i].first->vertex(edges[i].third);
		 
		/*vector<Point> test;
		test.push_back(v1->point());
		test.push_back(v2->point());*/

		/*drawMesh(dt,test);
		getchar();*/

		if( ( v1->info().is_corner() || v1->info().has_more_than_two_regions() || v1->info().is_singular_2() ) || 
			( v2->info().is_corner() || v2->info().has_more_than_two_regions() || v1->info().is_singular_2() ) )		continue;

		greater_edge_in_vertex(dt,v1);
		greater_edge_in_vertex(dt,v2);

		double d = v1->info().m_s_greater_edge/4;
		double cen = v2->info().m_s_greater_edge/4;
		if( d > cen )	d = cen;

		int label = common_label(v1,v2);
		
		vector<Vertex> star1 = v1->info().m_s_stars[index_label_in_star(v1,label)];
		vector<Vertex> star2 = v2->info().m_s_stars[index_label_in_star(v2,label)];

		V1 = v1->point() - CGAL::ORIGIN;
		V2 = v2->point() - CGAL::ORIGIN;		
		m = (V1 + V2)/2;	

		p1 = middle_point(star1);
		p2 = middle_point(star2);
		n1 = average_normal(star1,p1);
		n2 = average_normal(star2,p2);
		n = n1 + n2;
		n = unit_vector(n);
		
		int s1 = sgn( (V1-m)*n );
		int s2 = sgn( (V2-m)*n );

		p1 = V1 + s1 * n;
		p2 = V2 + s2 * n;	

		vol1 = star_volume(star1, V1, p1);
		vol2 = star_volume(star2, V2, p2);
		double D;
		if( abs(vol1) > abs(vol2) ){
			D = d*vol2/vol1;
			newV1 = V1 + D*s1*n;
			newV2 = V2 + d*s2*n;
			newV1 = V1 + D*n;
			newV2 = V2 + d*n;
		}
		else{
			D = d*vol1/vol2;
			newV1 = V1 + d*s1*n;
			newV2 = V2 + D*s2*n;
			newV1 = V1 + d*n;
			newV2 = V2 + D*n;
		}

		if( is_embedded_after_displacement(dt,v1,CGAL::ORIGIN + newV1) ){
			Point ant = v1->point();
			v1->set_point(CGAL::ORIGIN + newV1);
			if( is_embedded_after_displacement(dt,v2,CGAL::ORIGIN + newV2) ){
				//cout<<"mision complete"<<endl;
				//for drawing
				/*v1->set_point(ant);
				drawMesh(dt,test);
				getchar();
				v1->set_point(CGAL::ORIGIN + newV1);*/
				//
				v2->set_point(CGAL::ORIGIN + newV2);
			}
			else
				v1->set_point(ant);
		}

		/*cout<<"n: "<<n<<endl;
		cout<<"s1: "<<s1<<" s2: "<<s2<<endl;
		cout<<"p1: "<<p1<<" p2: "<<p2<<endl;
		cout<<"vol1: "<<vol1<<" vol2: "<<vol2<<endl;
		cout<<"d: "<<d<<" D: "<<D<<endl;
		cout<<"V1: "<<V1<<" V2: "<<V2<<endl;
		cout<<"newV1: "<<newV1<<" newV2: "<<newV2<<endl;*/

		//getchar();

		/*test.clear();
		test.push_back(CGAL::ORIGIN + newV1);
		test.push_back(CGAL::ORIGIN + newV2);
		drawMesh(dt,test);
		getchar();*/

		++counter;
	}
}

void restricted_balance_star(Delaunay &dt,Vertex &v, vector<Vertex> &star){
	double V=0,V1=0,h,hf;
	Vector p,n,m,x,Ver,newVer,point;
	vector<Vector> points;
	p = middle_point(star);
	n = average_normal(star,p);
	for(unsigned int i = 0; i < star.size(); ++i){
		h = ( (star[i]->point() - CGAL::ORIGIN) - p) * n;	
		point = (star[i]->point() - CGAL::ORIGIN)  - h*n;
		points.push_back(point);
		//cout<<"h: "<<h<<endl;
	}
	if(!is_point_in_polygon(p,points))	return;

	V = star_volume(star, p , v->point() - CGAL::ORIGIN);
	V1 = star_volume(star, p , p + n);
	hf = V/V1;
	newVer = p + hf*n;
	if( is_embedded_after_displacement(dt,v,CGAL::ORIGIN + newVer) )
		v->set_point(CGAL::ORIGIN + newVer);
} 

void restricted_balance3D( Delaunay &dt, vector<Vertex> &boundaryVertices ){
	for ( unsigned int i = 0 ; i < boundaryVertices.size() ; ++i){
		Vertex v = boundaryVertices[i] ;
		if( v->info().is_corner() || v->info().has_more_than_two_regions() || v->info().is_singular_2() )		continue;
		vector<Vertex> star = v->info().m_s_stars[0];
		restricted_balance_star(dt,v,star);
	}
}

void restricted_tsur(Delaunay &dt, int iteration){
	mark_boundary(dt);
	vector<Vertex> boundaryVertices;
	vector<Edge> edges;
	boundary_vertices(dt,boundaryVertices);
	edges_of_boundary(dt,edges);
	initiate_process(dt,boundaryVertices);
	
	for( int i = 0; i < iteration; ++i){
		restricted_balance3D(dt,boundaryVertices);
		cout<<"stop 1"<<endl;
		/*drawMesh(dt);
		getchar();*/
		
		restricted_TSUR3D(dt,edges);
		cout<<"stop 2"<<endl;
		/*drawMesh(dt);
		getchar();*/
	}
	restart_vertex_state(dt);
	cout<<"end tsur"<<endl;
}

void tsur(Delaunay &dt, int iteration){
	mark_boundary(dt);
	vector<Vertex> boundaryVertices;
	vector<Edge> edges;
	boundary_vertices(dt,boundaryVertices);

	edges_of_boundary(dt,edges);
	initiate_process(dt,boundaryVertices);

	clock_t begin = clock();
	
	for( int i = 0; i < iteration; ++i){
		balance3D(dt,boundaryVertices);
		cout<<"stop 1"<<endl;
		
		TSUR3D(dt,edges);		
		cout<<"stop 2"<<endl;

		/*drawMesh(dt);
		getchar();*/
	}
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<<"time tsur: "<<elapsed_secs<<endl;

	restart_vertex_state(dt);
	cout<<"end tsur"<<endl;
}

// for ODT

double ODT_star_volume( vector<Vertex> &star, Vector &v,Vector &p){
	Vector x1,x2,x3;
	double volTotal=0;
	for(unsigned int i = 0; i < star.size(); ++i){
		int next_i = ( i + 1 )%star.size();
		volTotal += CGAL::volume( CGAL::ORIGIN + v , star[i]->info().m_smoothed_point, star[next_i]->info().m_smoothed_point, CGAL::ORIGIN + p );
	}
	return volTotal;
}

void ODT_TSUR3D( Delaunay &dt, vector<Edge> &edges){
	Vertex v1,v2;
	double h,h1,h2,h1f,h2f,vol1,vol2,rlx=0;
	Vector p1,p2,n1,n2,n,m,V1,V2,newV1,newV2,pp1,pp2;	
	int counter = 0;
	bool draw=0;
	for ( unsigned int i = 0; i < edges.size(); ++i){
		v1 = edges[i].first->vertex(edges[i].second);
		v2 = edges[i].first->vertex(edges[i].third);
		 
		vector<Point> test;
		test.push_back(v1->point());
		test.push_back(v2->point());
		
		if( ( v1->info().is_corner() || v1->info().has_more_than_two_regions() || v1->info().is_singular_2() ) || 
			( v2->info().is_corner() || v2->info().has_more_than_two_regions() || v1->info().is_singular_2() ) )		continue;


		/*if( ( v1->info().is_corner()  || v1->info().is_singular_2() ) || 
			( v2->info().is_corner()  || v1->info().is_singular_2() ) )		continue;*/

		//if(v1->info().has_more_than_two_regions() || v2->info().has_more_than_two_regions() )	draw = 1;
		if(draw){
			drawMesh(dt,test);
			getchar();
		}

		int label = common_label2(v1,v2);
		
		vector<Vertex> star1 = v1->info().m_s_stars[index_label_in_star(v1,label)];
		vector<Vertex> star2 = v2->info().m_s_stars[index_label_in_star(v2,label)];

		V1 = v1->info().m_smoothed_point - CGAL::ORIGIN;
		V2 = v2->info().m_smoothed_point - CGAL::ORIGIN;		

		p1 = ODT_middle_point(star1);
		p2 = ODT_middle_point(star2);
		pp1 = p1;
		pp2 = p2;
		n1 = ODT_average_normal(star1,p1);
		n2 = ODT_average_normal(star2,p2);
		n = n1 + n2;
		n = unit_vector(n);
		m = (V1 + V2)/2;

		h1 = (V1 - m) * n;
		h2 = (V2 - m) * n;
		if(h1 == 0 || h2 == 0) continue;//indeterminate

		p1 = V1 - h1 * n;
		p2 = V2 - h2 * n;	

		vol1 = ODT_star_volume(star1, V1, p1);
		vol2 = ODT_star_volume(star2, V2, p2);

		double divisor = (vol1/h1 + vol2/h2);
		if( divisor == 0) continue;//indeterminate

		h =  2*(vol1 + vol2)/divisor;
		h1f = rlx * h1 + (1.0 - rlx) * h;
		h2f = rlx * h2 + (1.0 - rlx) * h;
		newV1 = p1 + h1f*n;
		newV2 = p2 + h2f*n;

		//if( !IsNumber( newV1.x() ) || !IsNumber( newV1.y() ) || !IsNumber( newV1.z() ) || !IsNumber( newV2.x() ) || 
		//!IsNumber( newV2.y() ) || !IsNumber( newV2.z() ) ) draw=1;

		/*v1->set_point(CGAL::ORIGIN + newV1);
		v2->set_point(CGAL::ORIGIN + newV2);*/

		v1->info().m_smoothed_point = CGAL::ORIGIN + newV1;
		v2->info().m_smoothed_point = CGAL::ORIGIN + newV2;
		
		/*if(draw){
			for(unsigned int i = 0; i < star1.size(); ++i)
				cout<<star1[i]->point()<<endl;
			for(unsigned int i = 0; i < star2.size(); ++i)
				cout<<star2[i]->point()<<endl;

			double test1 = star_volume(star1, V1, newV1);
			double test2 = star_volume(star2, V2, newV2);
		
			cout<<"counter: "<<counter<<endl;
			cout<<"LABEL: "<<label<<endl;
			cout<<"v1 size labels : "<<v1->info().m_regions_around.size()<<endl;
			cout<<"v2 size labels : "<<v2->info().m_regions_around.size()<<endl;
			cout<<"h1: "<<h1<<" h2: "<<h2<<endl;
			cout<<"h: "<<h<<" n: "<<n<<endl;
			cout<<"vol1: "<<vol1<<" vol2: "<<vol2<<endl;
			cout<<"test1: "<<test1<<" test2: "<<test2<<endl;
			cout<<"sum input: "<<vol1+vol2<<endl;
			cout<<"sum output: "<<test1+test2<<endl;
			cout<<"V1: "<<V1<<" V2: "<<V2<<endl;
			cout<<"newV1: "<<newV1<<" newV2: "<<newV2<<endl;
			cout<<endl;

			test.clear();
			test.push_back(CGAL::ORIGIN + newV1);
			test.push_back(CGAL::ORIGIN + newV2);
			drawMesh(dt,test);
			getchar();
			draw = 0;
		}
		++counter;*/
	}
}

void ODT_balance_star(Vertex &v, vector<Vertex> &star){
	double V=0,V1=0,h,hf;
	Vector p,n,m,x,Ver,newVer,point;
	vector<Vector> points;
	p = ODT_middle_point(star);
	n = ODT_average_normal(star,p);
	for(unsigned int i = 0; i < star.size(); ++i){
		h = ( (star[i]->info().m_smoothed_point - CGAL::ORIGIN) - p) * n;	
		point = (star[i]->info().m_smoothed_point - CGAL::ORIGIN)  - h*n;
		points.push_back(point);
	}
	if(!is_point_in_polygon(p,points))	return;

	V = ODT_star_volume(star, p , v->info().m_smoothed_point - CGAL::ORIGIN);
	V1 = ODT_star_volume(star, p , p + n);
	hf = V/V1;
	newVer = p + hf*n;
	//v->set_point(CGAL::ORIGIN + newVer);
	v->info().m_smoothed_point = CGAL::ORIGIN + newVer;
} 

void ODT_balance3D( Delaunay &dt, vector<Vertex> &boundaryVertices ){
	for ( unsigned int i = 0 ; i < boundaryVertices.size() ; ++i){
		Vertex v = boundaryVertices[i] ;
		if( v->info().is_corner() || v->info().has_more_than_two_regions()  || v->info().is_singular_2() )		continue;
		vector<Vertex> star = v->info().m_s_stars[0];
		ODT_balance_star(v,star);
	}
}

void initialize_smoothed_point(vector<Vertex> &boundary){
	for(unsigned int i = 0; i < boundary.size(); ++i){
		boundary[i]->info().m_smoothed_point = boundary[i]->point();
	}
}

void ODT_TSUR(Delaunay &dt, int iteration){
	mark_boundary(dt);
	vector<Vertex> boundaryVertices;
	vector<Edge> edges;
	boundary_vertices(dt,boundaryVertices);
	initialize_smoothed_point(boundaryVertices);

	edges_of_boundary(dt,edges);
	initiate_process(dt,boundaryVertices);

	for( int i = 0; i < iteration; ++i){
		ODT_balance3D(dt,boundaryVertices);
		cout<<"stop 1"<<endl;
		
		ODT_TSUR3D(dt,edges);		
		cout<<"stop 2"<<endl;
	}
	restart_vertex_state(dt);
	cout<<"end tsur"<<endl;
}