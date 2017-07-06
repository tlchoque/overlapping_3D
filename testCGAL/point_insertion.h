//new intersection

void boundary_facets_vertex_test(Delaunay &dt,vector<Facet> &facets,Vertex &v,vector<Cell> &boundary_cells){
	list<Cell> incident;
	Cell vc,c;
	int label = boundary_cells[0]->info().m_label;
	for(unsigned int i = 0; i < boundary_cells.size(); ++i){
		Cell vc = boundary_cells[i];
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

vector<Facet> ordered_facets(vector<Cell> &boundary_cells, Delaunay &dt,Vertex v){
	vector<Facet> ordered_facets;
	vector<Facet> facets;
	boundary_facets_vertex_test(dt,facets,v,boundary_cells);
	map<Vertex,Facet> vf;
	map<Vertex,Vertex> vv;
	for(unsigned int i=0;i < facets.size();++i){
		Facet f = facets[i];
		Cell c = f.first;		
		int index0 = f.second;
		int indexv = c->index(v);
		int index1 = dt.next_around_edge(index0,indexv);
		int index2 = dt.next_around_edge(index1,indexv);
		vv[c->vertex(index1)] = c->vertex(index2);
		vf[c->vertex(index1)] = f;
		if(i == 0){
			vv[v] = c->vertex(index1);
		}
	}
	Vertex it = v;
	for(unsigned int i=0;i < facets.size();++i){
		it = vv[it];
		Facet f = vf[it];
		ordered_facets.push_back(f);

	}
	return ordered_facets;
}

vector<Cell> order_boundary_cells(vector<Cell> &boundary_cells, Delaunay &dt,Vertex v){
	vector<Cell> ordered;
	for(unsigned int i = 0; i < boundary_cells.size(); ++i)
		boundary_cells[i]->info().m_boundary = true;

	vector<Facet> facets = ordered_facets(boundary_cells,dt,v);	
	for(unsigned int i = 0; i < facets.size(); ++i){
		Cell c = facets[i].first;
		if( c->info().m_boundary == true ){
			ordered.push_back(c);
			c->info().m_boundary = false;
		}	
	}
	return ordered;
}

//--------------

bool hole_cell(Cell &c,deque<Cell> &holes){
	c->info().m_label = HOLE;
	holes.push_back(c);
	//if( !c->info().m_hole_state ){// state point f the cell is in the deque		
	//	c->info().m_hole_state = true;
	//	holes.push_back(c);
	//	return true;
	//}		
	return false;
}

void decrease_decimals(Point &p){
	double dec= 1000;
	double x = floor(p.x()  * dec)/dec;
	double y = floor(p.y()  * dec)/dec;
	double z = floor(p.z()  * dec)/dec;
	p = Point(x,y,z);
}

int new_map(map<Point,int> &m, Point& p){
    std::map<Point,int>::const_iterator it = m.find( p );
    if (it == m.end())
        return -1;
    return it->second;
}

Cell adjacent_label(Delaunay &dt,Cell &c, Vertex &v,int label){
	Cell c1 = c->neighbor( dt.vertex_triple_index(c->index(v),0) );
	Cell c2 = c->neighbor( dt.vertex_triple_index(c->index(v),1) );
	Cell c3 = c->neighbor( dt.vertex_triple_index(c->index(v),2) );
	if( c1->info().m_label != -1 && c1->info().m_label != label)
		return c1;
	else if( c2->info().m_label != -1 && c2->info().m_label != label)
		return c2;
	else if( c3->info().m_label != -1 && c3->info().m_label != label)
		return c3;
	else return c;
}

Point average_neighbor_boundary_vertices(Delaunay &dt,vector<Cell> &cells, Vertex &v){
	Point p = v->point();
	vector<Vertex> vertices;
	vertices.push_back(v);
	Vertex v1;
	for(unsigned int i = 0; i < cells.size(); ++i){
		int idx = cells[i]->index(v);
		for(unsigned int j = 0; j < 3; ++j){
			v1 = cells[i]->vertex( dt.vertex_triple_index(idx,j) );
				if( v1->info().m_auxiliar == false ){
				p = Point(p.x() + v1->point().x(), p.y() + v1->point().y(),p.z() + v1->point().z() );
				v1->info().m_auxiliar = true;
				vertices.push_back(v1);
			}
		}
	}
	for(unsigned int i = 0; i < vertices.size(); ++i)
		vertices[i]->info().m_auxiliar = false;
	if(vertices.size() != 0 )
		p = Point( p.x() / vertices.size(),p.y() / vertices.size(),p.z() / vertices.size() );
	return p;
}

void relabel_with_neighbor(Delaunay &dt, deque<Cell> &cellsd,Vertex &v,Vertex &s, int label,deque<Cell> &holes){
	unsigned int sizec = 0;
	size_t cellSize = cellsd.size();	
	while(cellsd.size()!=0 && sizec < cellSize ){
		Cell main = cellsd.back();
		cellsd.pop_back();
		Cell c = main->neighbor(main->index(v));
		if(c->info().m_label != -1 && c->info().m_label != label){//get the oposite label	
			if( c->info().m_label == HOLE)		hole_cell(main,holes);
			else	main->info().m_label = c->info().m_label;				
			main->info().m_original_label = c->info().m_original_label;
		}
		else{// get the adjacent, but no oposite
			Cell idx = adjacent_label(dt,main,v,label);
			if(idx != main)	{
				if( idx->info().m_label == HOLE)	 hole_cell(main,holes);
				else	main->info().m_label = idx->info().m_label;
				main->info().m_original_label = idx->info().m_original_label;
			}
			else	{cellsd.push_front(main);++sizec;}
		}
	}
	//CHECK this step
	if( sizec >= cellSize){// there not exist and adjacent cel with different label 
		//cout<<" problem insert label"<<endl;
				
	/*	cout<<"no pos guau; label :"<<label<<endl;
		getchar();*/
		
		int alternateLabel=0;
		if(label == 0)	alternateLabel=1;
		for(unsigned int i = 1; i < s->info().m_star.size(); ++i){
			if( s->info().m_star[i].size() > s->info().m_star[alternateLabel].size()  && i!=label)
				alternateLabel = i;
		}	
		for(unsigned int i = 0; i <cellsd.size();++i){
			if( alternateLabel == HOLE) 	hole_cell(cellsd[i],holes);
			else  cellsd[i]->info().m_label = alternateLabel;
			cellsd[i]->info().m_original_label = alternateLabel;
		}
	}
}

void relabel_with_hash(Delaunay &dt,Vertex &v,Vertex &w, int label, map<Point,int> &pointLabel, deque<Vertex> &singularVertices, deque<Cell> &holes){
	deque<Cell> withoutLabel;
	list<Cell> incidentCellsW;
	dt.finite_incident_cells(w,back_inserter(incidentCellsW));
	for(std::list<Cell>::iterator it = incidentCellsW.begin(); it != incidentCellsW.end(); it++){
		Cell c = *it;
		Point cenc = CGAL::centroid(dt.tetrahedron(c));
		decrease_decimals(cenc);
		int mapLabel = new_map(pointLabel,cenc);
		if( mapLabel != -1){
			if( mapLabel == HOLE)	hole_cell(c,holes);
			else	c->info().m_label = mapLabel;
			c->info().m_original_label = mapLabel;
		}
		else
			withoutLabel.push_back(c);
	} 
	if(withoutLabel.size()!=0){
		relabel_with_neighbor(dt,withoutLabel,w,v,label,holes);
	}
}

void relabel_with_neighbor_repairing_1(Delaunay &dt, deque<Cell> &cellsd,Vertex &v,Vertex &s, int label,deque<Cell> &holes){
	for(unsigned int i = 0; i < cellsd.size(); ++i ){
		cellsd[i]->info().m_label = label;
	}
}

void relabel_with_hash_2(Delaunay &dt,Vertex &v,Vertex &w, int label, map<Point,int> &pointLabel, map<Point,int> &pointOriginalLabel, deque<Vertex> &singularVertices, deque<Cell> &holes,bool flag= false){
	deque<Cell> withoutLabel;
	list<Cell> incidentCellsW;
	dt.finite_incident_cells(w,back_inserter(incidentCellsW));
	for(std::list<Cell>::iterator it = incidentCellsW.begin(); it != incidentCellsW.end(); it++){
		Cell c = *it;
		Point cenc = CGAL::centroid(dt.tetrahedron(c));
		decrease_decimals(cenc);
		int mapLabel = new_map(pointLabel,cenc);
		if( mapLabel != -1){
			//cout<<"map label "<<mapLabel<<endl;
			if( mapLabel == HOLE)	hole_cell(c,holes);
			else	c->info().m_label = mapLabel;
			c->info().m_original_label = pointOriginalLabel[cenc];
		}
		else
			withoutLabel.push_back(c);
	} 
	if(withoutLabel.size()!=0){
		if(flag) {			
			relabel_with_neighbor_repairing_1(dt,withoutLabel,w,v,label,holes);
		}
		else	relabel_with_neighbor(dt,withoutLabel,w,v,label,holes);
	}
}

void hash_centroids_to_labels(Delaunay &dt, Vertex &v, vector<Cell> &cells, map<Point,int> &pointLabel,Point p){// cells doesn't have to be hashed
	list<Cell> incidentCells;
	dt.finite_incident_cells(v,back_inserter(incidentCells));
	Vertex v1,v2,v3;
	for(std::list<Cell>::iterator it = incidentCells.begin(); it != incidentCells.end(); it++){ 
		Cell c = *it;// if c is in cells continue
		if ( std::find(cells.begin(), cells.end(), c) != cells.end() )	continue; // c is in cells
		int idx = c->index(v);
		v1 = c->vertex( dt.vertex_triple_index(idx,0) );
		v2 = c->vertex( dt.vertex_triple_index(idx,1) );
		v3 = c->vertex( dt.vertex_triple_index(idx,2) );
		Point cen = CGAL::centroid(v1->point(),v2->point(),v3->point(),p);
		decrease_decimals(cen);
		pointLabel[cen] = c->info().m_label;
	}
}

void hash_centroids_to_labels_2(Delaunay &dt, Vertex &v, vector<Cell> &keep, map<Point,int> &pointLabel, map<Point,int> &pointOriginalLabel,Point p){// cells doesn't have to be hashed
	list<Cell> incidentCells;
	dt.finite_incident_cells(v,back_inserter(incidentCells));
	Vertex v1,v2,v3;
	for(std::list<Cell>::iterator it = incidentCells.begin(); it != incidentCells.end(); it++){ 
		Cell c = *it;// if c is in cells continue
		if ( std::find(keep.begin(), keep.end(), c) != keep.end() )	continue; // c is in cells
		int idx = c->index(v);
		v1 = c->vertex( dt.vertex_triple_index(idx,0) );
		v2 = c->vertex( dt.vertex_triple_index(idx,1) );
		v3 = c->vertex( dt.vertex_triple_index(idx,2) );
		Point cen = CGAL::centroid(v1->point(),v2->point(),v3->point(),p);
		decrease_decimals(cen);
		pointLabel[cen] = c->info().m_label;
		pointOriginalLabel[cen] = c->info().m_original_label;
	}
}

void check_if_cells_set_is_out(Delaunay &dt, vector<Cell> &cells, Vertex v,int label, deque<Vertex> &singularVertices){
	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		Cell c = *it;
		if ( std::find(cells.begin(), cells.end(), c) == cells.end() )	continue;// c is not in cells
		double a2 = 0;int lbl = 0;int idx=0; // there will always be a facet 
		for(unsigned int j = 0; j < 3;++j){
			idx = dt.vertex_triple_index(c->index(v),j);
			if( dt.triangle( Facet(c, idx) ).squared_area() > a2 && c->neighbor( idx )->info().m_label != label && c->neighbor( idx )->info().m_label != -1) {
				a2 = dt.triangle( Facet(c, idx) ).squared_area();
				lbl = c->neighbor( idx )->info().m_label;
			}
		}// at the end, we got the label & index			
		c->info().m_label = lbl;
		c->info().m_state = 1;	

		int index = c->index(v);
		for(unsigned int j = 0; j < 3;++j){
			Vertex vc = c->vertex( dt.vertex_triple_index(index,j) );
			update_regions_around_vertex( dt, vc );
			if( vc->info().m_state != 2  && vc->info().is_singular() ){
				singularVertices.push_front( vc );
				vc->info().m_state=2;
			}
		}	
	}
	// not update v because v is uncomplete
}

void detect_new_singularities(Delaunay &dt,Vertex w,deque<Vertex> &singularVertices){
	list<Vertex> verticesW;
	dt.finite_incident_vertices(w,back_inserter(verticesW));
	for(std::list<Vertex>::iterator it = verticesW.begin(); it != verticesW.end(); it++){ 
		Vertex v = *it;
		update_regions_around_vertex(dt, v );
		if( v->info().m_state != 2 && v->info().is_singular() ){
			v->info().m_state=2;
			singularVertices.push_back(v);
		}
	}
	update_regions_around_vertex(dt, w );
}

// new insertion point

vector<Cell> boundary_cells(Delaunay &dt,vector<Cell> &erode, Vertex &v){
	int label = erode[0]->info().m_label;
	vector<Cell> cells;
	Cell vc,c0;
	for( unsigned int i = 0; i < erode.size(); ++i){
		vc = erode[i];
		int index = vc->index(v);
		for(unsigned int i = 0; i < 3;++i){
			c0 =  vc->neighbor( dt.vertex_triple_index(index,i) );
			if( c0->info().m_label != label ){
				cells.push_back(vc);
				break;
			}
		}
	}
	return cells;
}

Point intersection_1(Delaunay &dt, vector<Cell> &erode, Vertex &v){// one cell
	int index = erode[0]->index(v);
	Cell c = erode[0];

	Point t= CGAL::centroid( c->vertex( dt.vertex_triple_index(index,0) )->point(),  
		c->vertex( dt.vertex_triple_index(index,1) )->point(), c->vertex( dt.vertex_triple_index(index,2) )->point() );

	return CGAL::ORIGIN + ( v->point() - CGAL::ORIGIN ) + 0.4*( t - v->point() );


	/*return CGAL::centroid( c->vertex( dt.vertex_triple_index(index,0) )->point(),  
		c->vertex( dt.vertex_triple_index(index,1) )->point(), c->vertex( dt.vertex_triple_index(index,2) )->point() );*/
}

Point intersection_2(Delaunay &dt, vector<Cell> &erode, Vertex &v){
	//cout<<"insrtion 2"<<endl;
	Cell c = erode[0],c0;
	int index = erode[0]->index(v);
	int label = erode[0]->info().m_label;
	Vector middle;
	for(unsigned int i = 0; i < 3;++i){
		int j = dt.vertex_triple_index(index,i);
		c0 =  c->neighbor( j );
		if( c0->info().m_label == label ){
			int index1 = dt.next_around_edge(j,index);
			int index2 = dt.next_around_edge(index1,index);
			middle = ( (c->vertex(index1)->point() - CGAL::ORIGIN) + ( c->vertex(index2)->point() - CGAL::ORIGIN) )/2;
			break;
		}
	}
	//middle = ( v->point() - CGAL::ORIGIN ) + 0.85*( ( CGAL::ORIGIN + middle ) - v->point() ); // change to 0.5
	middle = ( v->point() - CGAL::ORIGIN ) + 0.6*( ( CGAL::ORIGIN + middle ) - v->point() );
	return CGAL::ORIGIN + middle;
}

Point n_intersection_3(Delaunay &dt, vector<Cell> &erode, Vertex &v){
	//cout<<"size of vector to erode: "<<erode.size()<<endl;
	vector<Cell> boundary = boundary_cells(dt,erode,v);
	Point c,d,m,x;
	Vector cv,cd,cm,ce,vx,test;
	Point intersection;
	double distance = 10000000;

	for( unsigned int i = 0; i < boundary.size() - 1; ++i){
		c = boundary[i]->circumcenter();
		cv = v->point() - c;
		for( unsigned int j = i+1; j < boundary.size(); ++j){
			d = boundary[j]->circumcenter(); // d->circumcenter 
			cd = unit_vector(d - c);
			cm = ( cv*cd )*cd;		
			ce = cm-cv;// ce -> vector parallel to mx
			m = c + cm;// m-> middle point between v and x  HERE i can make short the computation
			x = m + ce;
			vx = x - v->point();
			if( vx*vx < distance){
				intersection = x;
				distance = vx*vx;
			}
		}
	}	
	//intersection = (CGAL::ORIGIN + ( v->point() - CGAL::ORIGIN ) + 0.4*(intersection - v->point() ) );	//originally 0.8 , for graphics 0.3  x-1.5 z +1s
	
	intersection = (CGAL::ORIGIN + ( v->point() - CGAL::ORIGIN ) + 0.8*(intersection - v->point() ) );
	
	for( unsigned int i = 0; i < erode.size(); ++i){
		Cell c = erode[i];
		//if( CGAL::side_of_bounded_sphere( c->vertex(0)->point(), c->vertex(1)->point(), 
		//	c->vertex(2)->point(), c->vertex(3)->point(), intersection ) != CGAL::ON_BOUNDED_SIDE ){
		//	//cout<<"not intersection: "<<endl;
		//	return v->point();
		//}

		if( dt.side_of_sphere(c,intersection,true) != CGAL::ON_BOUNDED_SIDE){
			/*cout<<" not on bounded side "<<dt.tetrahedron(c)<<endl;
			cout<<"label " <<c->info().m_label<<endl;*/
			/*vector<Cell> vt;
			vt.push_back(c);
			drawCells(dt,vt);
			getchar();*/

			
			/*vector<Point> vt;
			Point a = c->vertex(0)->point();
			Point b = c->vertex(1)->point();
			Point d = c->vertex(2)->point();
			Point e = c->vertex(3)->point();
			vt.push_back(intersection);
			vt.push_back(a);
			vt.push_back(b);
			vt.push_back(d);
			vt.push_back(e);
			drawMesh(dt,vt);
			getchar();*/
			
			return v->point();
		}
	}
	return intersection;
}

Point center_circle_intersection(Vector n, Point P, Point T){
	Vector d = T - P;
	Vector proy = (n * d )*n;
	Point Q = P + ( d - proy);
	return Q;
}

void sphere_sphere_intersection(Point a, Point b, Point v, Point &P,Vector &n){
	Vector dot = b - a;
	Vector r = a - v;
	Vector R = b - v;
	double h = 0.5  + 0.5*(r*r - R*R )/(dot*dot);
	P = a + h*(dot);
	n = unit_vector(b - P) ;
}

Point circle_circle_intersection(Point a, Point b, Point v){
	Vector dot = b - a;
	Vector r = a - v;
	Vector R = b - v;
	double h = 0.5  + 0.5*(r*r - R*R )/(dot*dot);
	Point P = a + h*(dot);
	return v + 2*(P - v) ;
}

Point intersection_3(Delaunay &dt, vector<Cell> &erode, Vertex &v){
	//if(erode.size() > 5 ) cout<<"greater than 5: "<<erode.size()<<endl;

	vector<Cell> boundary = boundary_cells(dt,erode,v);
	vector<Cell> ordered = order_boundary_cells(boundary,dt,v);	
	Point intersection = v->point();

	vector<Point> pts;	
	double distance = 10000000000000000;
	for(unsigned int i = 0; i < ordered.size(); ++i){
		int next_index;
		if( i == ordered.size() - 1)		next_index = 0;
		else	next_index = i + 1;
		Point a = ordered[i]->circumcenter();
		Point b = ordered[next_index]->circumcenter();
		Point P; Vector n;
		sphere_sphere_intersection(a,b,v->point(),P,n);
		//cout<<"iteration: "<<i<<endl;
		pts.clear();
		pts.push_back( v->point() );
		pts.push_back( a );
		pts.push_back( b );
		pts.push_back(P);
		/*cells.push_back( ordered[i] );
		cells.push_back( ordered[next_index] );*/
		vector<pair<Point,double>>  circles;
		double radius = norm(a - v->point() );
		double radius2 = norm(b - v->point() );
		circles.push_back( make_pair( a , radius ) );
		circles.push_back( make_pair( b , radius2 ) );

		/*drawMesh(dt,pts);
		draw_spheres(circles);
		getchar();*/

		bool cen = 0;
		double min_angle;		
		Vector plane_normal;
		Point x;
		for(unsigned int j = 0; j < ordered.size(); ++j){
			if( j == i || j == next_index)	continue;
			Cell t = ordered[j];
			Point T = t->circumcenter();			
			Point Q = center_circle_intersection(n,P,T); // proyection of T in plane of circle P 
			
			Point O = circle_circle_intersection(P,Q, v->point());// interesction between circle P and Q			
			//for cross angle sign

			Vector vp_normalized = unit_vector( v->point() - P );
			Vector op_normalized = unit_vector( O - P );
			Vector cr = CGAL::cross_product( vp_normalized , op_normalized );
			double angle = acos( vp_normalized*op_normalized );

			/*pts.push_back( T );
			pts.push_back( Q );*/
			pts.push_back( O );

			double rr = norm( T - v->point() );
			//circles.push_back( make_pair( T, rr) );

			/*drawMesh(dt,pts);
			getchar();*/

			if( cen == 0){
				min_angle = angle;
				x = O;
				plane_normal = unit_vector(cr);
				cen = 1;
			}
			else{
				if( cr*plane_normal < 0 )  {
					x = v->point();
					break;
				}
				if( angle < min_angle ){
					min_angle = angle;
					x = O;					
				}				
			}
		}
		if( x != v->point() ){
			Vector dis = x - v->point();
			double val = dis*dis;
			if( val < distance){
				intersection = x;
				distance =  val;
			}
		}
	}
	
	pts.push_back( intersection );
	
	/*drawMesh(dt,pts);
	getchar();*/

	intersection = (CGAL::ORIGIN + ( v->point() - CGAL::ORIGIN ) + 0.8*(intersection - v->point() ) );
	if( intersection == v->point() ) return intersection;

	//cout<<"v and intersection before:  "<<v->point()<<" * "<<intersection<<endl;

	for( unsigned int i = 0; i < erode.size(); ++i){
		Cell c = erode[i];
		if( dt.side_of_sphere(c,intersection,true) != CGAL::ON_BOUNDED_SIDE){			
			//cout<<"side_of_sphere error"<<endl;
			return v->point();
		}
	}
	//cout<<"v and intersection after:  "<<v->point()<<" * "<<intersection<<endl;
	/*cout<<"stop 3"<<endl;	
	drawMesh(dt,pts);
	getchar();*/
	
	return intersection;
}

Point intersection_path(Delaunay &dt, vector<Cell> &path, Vertex &v){
	Point intersection;
	Point c,d,m,x;
	Vector cv,cd,cm,ce,vx;

	c = path.front()->circumcenter();// think
	d = path.back()->circumcenter();	
	cv = v->point() - c;


	//d = path[1]->circumcenter();
	cd = unit_vector(d - c);	
	cm = ( cv*cd )*cd;		
	

	ce = cm-cv;// ce -> vector parallel to mx
	m = c + cm;// m-> middle point between v and x  HERE i can make short the computation
	x = m + ce;
	intersection = x;

	intersection = (CGAL::ORIGIN + ( v->point() - CGAL::ORIGIN ) + 0.5*(intersection - v->point() ) );
	//cout<<"intersection after: "<<intersection<<endl;

	for( unsigned int i = 0; i < path.size(); ++i){
		Cell c = path[i];
		if( CGAL::side_of_bounded_sphere( c->vertex(0)->point(), c->vertex(1)->point(), 
			c->vertex(2)->point(), c->vertex(3)->point(), intersection ) != CGAL::ON_BOUNDED_SIDE ){
			return v->point();
		}
	}
	return intersection;
}

bool circumsphere_intersection(Delaunay &dt,vector<Cell> &erode, Vertex v, Point &collision){
	if( erode.size() == 1 ){
		//cout<<"# of cells 1"<<endl;
		collision = intersection_1(dt,erode,v);
	}
	else if( erode.size() == 2  ){
		//cout<<"# of cells 2"<<endl;
		collision = intersection_2(dt,erode,v);
	}
	else{
		//cout<<"# of cells +2"<<endl;
		collision = intersection_3(dt,erode,v);
	}
	if(collision == v->point() ) return false;
	return true;
}

bool circumsphere_intersection_path(Delaunay &dt,vector<Cell> &erode, Vertex v, Point &collision){
	//cout<<"size path "<<erode.size()<<endl;
	if( erode.size() == 1 )
		collision = intersection_1(dt,erode,v);
	else if( erode.size() == 2  )
		collision = intersection_2(dt,erode,v);
	else
		collision = intersection_path(dt,erode,v);
	if(collision == v->point() ) return false;
	else return true;
}

void unvisited(Delaunay &dt, Vertex v, vector<Cell> &preserved_cells ){
	list<Cell> cells;
	dt.finite_incident_cells(v,back_inserter(cells));
	for(std::list<Cell>::iterator it = cells.begin(); it != cells.end(); it++){ 
		if( !(*it)->info().m_visited ){
			int index = (*it)->index(v);
			Cell c = (*it)->neighbor(index);
			if( !dt.is_infinite(c) )
				preserved_cells.push_back(c);
		}
		else
			(*it)->info().m_visited = false;
	}
}

vector<Cell> mark_visited_boundary_cells(Delaunay &dt,vector<Cell> &erode, Vertex &v){
	int label = erode[0]->info().m_label;
	vector<Cell> cells;
	Cell vc,c0;
	for( unsigned int i = 0; i < erode.size(); ++i){
		vc = erode[i];
		vc->info().m_visited=true;
		int index = vc->index(v);
		for(unsigned int i = 0; i < 3;++i){
			c0 =  vc->neighbor( dt.vertex_triple_index(index,i) );
			if( c0->info().m_label != label ){
				cells.push_back(vc);
				break;
			}
		}
	}
	return cells;
}

bool check_sphere_collision(Delaunay &dt, vector<Cell> &collision_cells, Point &collision,Vertex v){
	double t = 1,ta;
	for( unsigned int i = 0; i < collision_cells.size(); ++i){
		Cell c = collision_cells[i];
		Point center = c->circumcenter();
		Vector r = center - c->vertex(0)->point();
		double radio = sqrt( r * r );
		int step = getNearPoint(v->point(),collision,center,radio,ta);
		if(step){
			if(ta < t)  t = ta;
		}
		//cout<<"stop check_sphere_collision -1 4"<<endl;
	}
	if(t > 0 && t <= 1){	
		t = 0.98*t;
		collision = v->point() + t*(collision - v->point()); 
		//cout<<v->point()<<" * "<<collision<<" * "<<t<<endl;
		if( collision == v->point() )
			return 0;
		else	return 1;
	}
	else{
		collision = v->point();
		return 0;
	}
}

void check_cells_in_conflict(Delaunay &dt, vector<Cell> &collision_cells, int erodelabel, Vertex v, Point &collision, vector<Cell> &cell_out_star){
	for( unsigned int i = 0;i < collision_cells.size(); ++i ){
		//cout<<"stop conflict 0 "<<endl;
		Cell c = collision_cells[i];
		if( c->has_vertex(v) )  continue;
		int label = c->info().m_label;
		if( dt.side_of_sphere(c,collision,false) == CGAL::ON_BOUNDED_SIDE ){
			//cout<<"stop conflict 1 "<<endl;
			if( erodelabel == label ){			
				collision_cells.erase(collision_cells.begin() + i);
				--i;
				for(unsigned int j = 0; j < 4 ; ++j){	
					//cout<<"stop conflict 2 "<<endl;
					Cell vc = c->neighbor(j);
					if( !vc->has_vertex(v) && !dt.is_infinite(vc) ){
						//collision_cells.push_back(vc);
						cell_out_star.push_back(vc);
					}
				}
			}
		}
	}
}

bool check_sphere_collision2(Delaunay &dt, vector<Cell> &collision_cells, Point &collision,Vertex v,bool draw){
	double t = 1,ta;
	//cout<<"check_sphere_collision2 0 "<<endl;
	for( unsigned int i = 0; i < collision_cells.size(); ++i){
		//cout<<i <<" :: stop -1"<<endl;
		
		Cell c = collision_cells[i];
		//cout<<"tetra: "<<dt.tetrahedron(c)<<endl;
		//cout<<"stop 0"<<endl;
		
		if( dt.is_infinite(c)){
			cout<<"is infinite "<<endl;
			getchar();
		}

		/*if(draw && i == 4){
			cout<<"tetra: "<<dt.tetrahedron(c)<<endl;
			vector<Point> pts;
			pts.push_back( c->vertex(0)->point() );
			pts.push_back( c->vertex(1)->point() );
			pts.push_back( c->vertex(2)->point() );
			pts.push_back( c->vertex(3)->point() );
			vector<Cell> cells;
			cells.push_back(c);
			drawMesh(dt,pts);
			drawCells(dt,cells);
			getchar();
		}*/

		Point center = c->circumcenter();
		//cout<<"stop 1"<<endl;	
		Vector r = center - c->vertex(0)->point();
		//cout<<"stop 2"<<endl;
		double radio = sqrt( r * r );
		//cout<<"stop 3"<<endl;	
		
		int step = getNearPoint(v->point(),collision,center,radio,ta);

		//cout<<"stop 4"<<endl;
		//cout<<"check_sphere_collision2 2 collision: "<<collision<<endl;
		// begin the draw test
		/*vector<Point> pts;
		pts.push_back( v->point() );
		pts.push_back( v->point() + t*(collision - v->point()) );		*/	

		if(step)	if(ta < t)  t = ta;
		//cout<<"stop 5"<<endl;

		//cout<<"check_sphere_collision2 2 5 "<<endl;

		/*if( draw){			
			pts.push_back( v->point() + t*(collision - v->point()) );			
			pts.push_back( c->vertex(0)->point() );
			pts.push_back( c->vertex(1)->point() );
			pts.push_back( c->vertex(2)->point() );
			pts.push_back( c->vertex(3)->point() );
			pts.push_back( c->circumcenter() );
			cout<<"this is collision, label "<<c->info().m_label<<endl;
			drawMesh(dt,pts);
			getchar();
		}*/
	}

	//cout<<"check_sphere_collision2 3 "<<endl;

	if(t > 0 && t <= 1){	
		t = 0.98*t;
		collision = v->point() + t*(collision - v->point()); 
		if( collision == v->point() )
			return 0;
		else	return 1;
	}
	else{
		collision = v->point();
		return 0;
	}
}

bool check_preserved_cells(Delaunay &dt, vector<Cell> &preserved_cells, Point &collision,Vertex v){
	for( unsigned int i = 0; i < preserved_cells.size(); ++i){
		Cell c = preserved_cells[i];
		//if( CGAL::side_of_bounded_sphere( c->vertex(0)->point(), c->vertex(1)->point(), 
		//	c->vertex(2)->point(), c->vertex(3)->point(), collision ) != CGAL::ON_UNBOUNDED_SIDE ){
		//		//cout<<"preserved: "<<c->info().m_label<<endl;
		//	return false;
		//}

		//cout<<"tetra "<<dt.tetrahedron(c)<<endl;
		if( dt.side_of_sphere(c,collision,false) == CGAL::ON_BOUNDED_SIDE ){
			return false;
		}
		//else cout<<"preserved"<<endl;
	}
	return true;
}

bool check_preserved_cells2(Delaunay &dt, vector<Cell> &preserved_cells, Point &collision,Vertex v, bool draw){
	for( unsigned int i = 0; i < preserved_cells.size(); ++i){
		Cell c = preserved_cells[i];
		if( dt.side_of_sphere(c,collision,false) == CGAL::ON_BOUNDED_SIDE ){
			/*if( draw){
				cout<<"here"<<endl;
				vector<Point> cellpoints;
				cellpoints.push_back( collision );
				cellpoints.push_back( v->point() );
				cellpoints.push_back(c->vertex(0)->point() );
				cellpoints.push_back(c->vertex(1)->point() );
				cellpoints.push_back(c->vertex(2)->point() );
				cellpoints.push_back(c->vertex(3)->point() );
				cellpoints.push_back(c->circumcenter() );
				drawMesh(dt,cellpoints);
				getchar();
			}		*/	
			return false;
		}
		//else cout<<"preserved"<<endl;
	}
	return true;
}

void select_preserved_cells(Delaunay &dt, vector<Cell> &preserved,vector<Cell> &erode, Vertex v){
	for(unsigned int i = 0; i < preserved.size(); ++i){
		Cell c = preserved[i];
		int vindex = c->index(v);
		for( unsigned int j = 0; j < 3; ++j){			
			int index = dt.vertex_triple_index(vindex,j);
			Cell vc = c->neighbor(index);
			if ( std::find(erode.begin(), erode.end(), vc) != erode.end() )	{ // c is in erode, then take out vc from vector keep
				c->info().m_visited = false;
				preserved.erase( preserved.begin() + i ) ;
				--i;
				break;
			}
		}
	}
}

bool safe_point(Delaunay &dt, vector<Cell> &erode, vector<Cell> &keep, Vertex v, Point &collision){	
	//cout<<"safe stop 0 "<<endl;
	if( !circumsphere_intersection(dt,erode,v, collision) ) { /*cout<<"not intersection: "<<endl; */return false;}
	//cout<<"safe stop 1 "<<endl;
	vector<Cell> preserved_cells = mark_visited_boundary_cells(dt,keep,v);
	vector<Cell> collision_cells;
	unvisited(dt,v,collision_cells);
	//cout<<"safe stop 2 "<<endl;
	if( !check_sphere_collision(dt,collision_cells,collision,v) ) { /*cout<<"collision: "<<endl;*/ return false;}
	//cout<<"safe stop 3 "<<endl;

	if( !check_preserved_cells(dt,preserved_cells,collision,v) ) { /*cout<<"preserved: "<<endl;*/ return false;}

	//if( !check_preserved_cells2(dt,preserved_cells,collision,v,draw) ) { cout<<"preserved: "<<endl; return false;}
	//cout<<"safe stop 4 "<<endl;
	return true;
}

bool safe_point2(Delaunay &dt, vector<Cell> &erode, vector<Cell> &keep, Vertex v, Point &collision, vector<Cell> &cells_out_star,bool draw){	
	//cout<<"safe stop 0 "<<endl;
	cells_out_star.clear();
	if( !circumsphere_intersection(dt,erode,v, collision) ) {/*cout<<"not intersection: "<<endl;*/ return false;}
	//cout<<"safe stop 1 "<<endl;
	vector<Cell> preserved_cells = mark_visited_boundary_cells(dt,keep,v);
	vector<Cell> collision_cells;

	//select_preserved_cells(dt,preserved_cells,erode,v);
	unvisited(dt,v,collision_cells);

	//cout<<"safe stop 2 "<<endl;
	vector<Point> cellpoints;

	/*if(draw){
		cout<<"size: "<<erode.size()<<endl;
		cellpoints.push_back( collision );
		cellpoints.push_back( v->point() );
		drawMesh(dt,cellpoints);
		getchar();
	}*/

	check_cells_in_conflict(dt,collision_cells,erode[0]->info().m_label,v,collision,cells_out_star);

	//if( check_preserved_cells2(dt,preserved_cells,collision,v,draw) ) { cout<<"preserved: "<<endl;}

	if( !check_sphere_collision2(dt,collision_cells,collision,v,draw) ) { /*cout<<"collision: "<<endl;*/ return false;}
	if( !check_sphere_collision2(dt,cells_out_star,collision,v,draw) ) { /*cout<<"collision2 : "<<endl;*/ return false;}
	//cout<<"safe stop 3 "<<endl;

	if( !check_preserved_cells2(dt,preserved_cells,collision,v,draw) ) { /*cout<<"not preserved: "<<endl;*/ return false;}

	//cout<<"safe stop 4 "<<endl;
	//if( !check_preserved_cells2(dt,preserved_cells,collision,v,draw) ) { cout<<"preserved: "<<endl; return false;}
	//cout<<"safe stop 4 "<<endl;
	return true;
}

void detect_new_singularities_2(Delaunay &dt,Vertex w,deque<Vertex> &singularVertices){
	list<Vertex> verticesW;
	dt.finite_incident_vertices(w,back_inserter(verticesW));
	for(std::list<Vertex>::iterator it = verticesW.begin(); it != verticesW.end(); it++){ 
		Vertex v = *it;
		update_vertex_info(dt, v );
		if( v->info().m_state != 2 && v->info().is_singular_2() ){
			v->info().m_state=2;
			if( v->info().has_two_regions() )
				singularVertices.push_back(v);
			else
				singularVertices.push_front(v);
		}
	}
	update_vertex_info(dt, w );
	if( w->info().m_state != 2 && w->info().is_singular_2() ){
		w->info().m_state=2;
		//singularVertices.push_back(w);
		if( w->info().has_two_regions() )
			singularVertices.push_back(w);
		else
			singularVertices.push_front(w);
	}
}

bool insert_point_1(Delaunay &dt, deque<Vertex> &singularVertices, vector<vector<Cell>> &cellSets,Vertex &v,int label,int min,int max,deque<Cell> &holes,bool draw ){
	Point pv2 = average_neighbor_boundary_vertices(dt,cellSets[max],v);
	Point pv1 = average_neighbor_boundary_vertices(dt,cellSets[min],v);
	int step2 = getDisplacedPoint(dt,v,pv2);
	if(step2){
		map<Point,int> pointLabel;
		hash_centroids_to_labels(dt,v,cellSets[max],pointLabel,pv2);	
		Vertex w = dt.insert(pv2);	
		check_if_cells_set_is_out(dt,cellSets[min],v,label,singularVertices);		
		relabel_with_hash(dt,v,w,label,pointLabel,singularVertices,holes);
		detect_new_singularities_2(dt,w,singularVertices);
		update_vertex_info(dt, v );

		if(draw){
			aroundVertex(dt,v,w,label);
			getchar();
		}
		return true;
	}
	else if( getDisplacedPoint(dt,v,pv1) ){
		map<Point,int> pointLabel;
		hash_centroids_to_labels(dt,v,cellSets[min],pointLabel,pv1);	
		Vertex w = dt.insert(pv1);	
		check_if_cells_set_is_out(dt,cellSets[max],v,label,singularVertices);		
		relabel_with_hash(dt,v,w,label,pointLabel,singularVertices,holes);
		detect_new_singularities_2(dt,w,singularVertices);
		update_vertex_info(dt, v );

		if(draw){
			aroundVertex(dt,v,w,label);
			getchar();
		}
		return true;
	}
	else{
		//cout<<"not possible inserted"<<endl;
		return false;
	}
}

bool insert_point_2(Delaunay &dt, deque<Vertex> &singularVertices, vector<vector<Cell>> &cellSets,Vertex &v,int label,int min,int max,deque<Cell> &holes,bool draw , bool flag = false){
	Point p;
	/*cout<<"vrai insert_point_2"<<endl;
	if(flag)	cout<<"vrai flag"<<endl;*/
	if( safe_point(dt,cellSets[max],cellSets[min],v,p) ){
		//cout<<" point inserted 1"<<endl;
		map<Point,int> pointLabel;
		map<Point,int> pointOriginaLabel;
		hash_centroids_to_labels_2(dt,v,cellSets[min],pointLabel,pointOriginaLabel,p);	
		Vertex w = dt.insert(p, cellSets[max][0]);	
		//relabel_with_hash_2(dt,v,w,cellSets[min][0]->info().m_label,pointLabel,pointOriginaLabel,singularVertices,holes);
		relabel_with_hash_2(dt,v,w,label,pointLabel,pointOriginaLabel,singularVertices,holes,flag);
		detect_new_singularities_2(dt,w,singularVertices);
		update_vertex_info(dt, v );

		if(draw){
			/*vector<Vertex> vertices;
			vertices.push_back(v);
			vertices.push_back(w);
			save_vertex_star(dt,vertices,"identification/insertion/2comp","erode");*/

			aroundVertex(dt,v,w,label);
			//drawMesh_star(dt);
			//save_vertex_star(dt,"identification/insertion/1comp","dilate");
			getchar();
		}
		//cout<<" point inserted 1"<<endl;
		return true;
	}	
	else if( safe_point(dt,cellSets[min],cellSets[max],v,p) ){
		//cout<<" point inserted 2"<<endl;
		map<Point,int> pointLabel;
		map<Point,int> pointOriginaLabel;
		/*hash_centroids_to_labels(dt,v,cellSets[max],pointLabel,p);	
		Vertex w = dt.insert(p, cellSets[max][0]);	
		relabel_with_hash(dt,v,w,cellSets[max][0]->info().m_label,pointLabel,singularVertices,holes);*/
		hash_centroids_to_labels_2(dt,v,cellSets[max],pointLabel,pointOriginaLabel,p);	

		/*for(unsigned int i = 0; i < cellSets[max].size(); ++i){
			cout<<i<<" tetra:"<<dt.tetrahedron(cellSets[max][i])<<endl;
			cout<<" label:"<<cellSets[max][i]->info().m_label<<endl;
		}*/

		Vertex w = dt.insert(p, cellSets[min][0]);	
		//relabel_with_hash_2(dt,v,w,cellSets[max][0]->info().m_label,pointLabel,pointOriginaLabel,singularVertices,holes);
		relabel_with_hash_2(dt,v,w,label,pointLabel,pointOriginaLabel,singularVertices,holes,flag);

		/*for(unsigned int i = 0; i < cellSets[max].size(); ++i){
			cout<<i<<" tetra:"<<dt.tetrahedron(cellSets[max][i])<<endl;
			cout<<" label:"<<cellSets[max][i]->info().m_label<<endl;
		}*/

		detect_new_singularities_2(dt,w,singularVertices);
		update_vertex_info( dt, v );
		if(draw){
			/*
			vector<Vertex> vertices;
			vertices.push_back(v);
			vertices.push_back(w);
			save_vertex_star(dt,vertices,"identification/insertion/2comp","erode");*/

			aroundVertex(dt,v,w,label);
			//drawMesh_star(dt);
			getchar();
		}
		//cout<<" point inserted 2"<<endl;
		return true;
	}
	else{
		//cout<<"not possible inserted"<<endl;
		return false;
	}
}

// separate 1 region from all the others
bool insert_point_3(Delaunay &dt, deque<Vertex> &singularVertices, vector<vector<Cell>> &cellSets,Vertex &v,int label,int min,int max,deque<Cell> &holes,bool draw ){
	Point p;
	if( safe_point(dt,cellSets[min],cellSets[max],v,p) ){
		//cout<<"inserting point "<<endl;
		map<Point,int> pointLabel;
		hash_centroids_to_labels(dt,v,cellSets[max],pointLabel,p);	
		Vertex w = dt.insert(p, cellSets[max][0]);	
		relabel_with_hash(dt,v,w,label,pointLabel,singularVertices,holes);
		detect_new_singularities_2(dt,w,singularVertices);
		update_vertex_info(dt, v );
		if(draw){
			aroundVertex(dt,v,w,label);
			getchar();
		}
		return true;
	}
	else{
		//cout<<"not possible inserted"<<endl;
		return false;
	}
}

// point insertion to dilate
vector<Cell> complement_cells(Delaunay &dt,vector<Cell> &erode, Vertex &v){
	vector<Cell> complement;
	list<Cell> cells;
	int label = erode[0]->info().m_label;
	dt.finite_incident_cells(v,back_inserter(cells));
	for(std::list<Cell>::iterator it = cells.begin(); it != cells.end(); it++){ 
		Cell c = *it;
		if( c->info().m_label == label){
			if( std::find( erode.begin(), erode.end(), c) == erode.end() ){//not found
				complement.push_back(c);
			}
		}		
	}
	return complement;
}

vector<Cell> mark_visited_boundary_cells_path(Delaunay &dt,vector<Cell> &erode, Vertex &v){
	int label = erode[0]->info().m_label;
	vector<Cell> cells;
	Cell vc,c0;
	for( unsigned int i = 0; i < erode.size(); ++i){
		vc = erode[i];
		//vc->info().m_visited=true;
		int index = vc->index(v);
		for(unsigned int j = 0; j < 3;++j){
			c0 =  vc->neighbor( dt.vertex_triple_index(index,j) );
			if( c0->info().m_label != label ){				
				vc->info().m_visited = true;
				cells.push_back(vc);
				break;
			}
		}
	}
	return cells;
}

bool safe_point_to_dilate(Delaunay &dt, vector<Cell> &erode, vector<Cell> &keep, Vertex v, Point &collision, vector<Cell> &cell_out_star,bool draw){	
	//cout<<"safe_point_to_dilate stop 0 "<<endl;
	if( !circumsphere_intersection_path(dt,erode,v, collision) ) { /*cout<<"not intersection: "<<endl;*/ return false;}
	vector<Cell> preserved_cells = mark_visited_boundary_cells_path(dt,keep,v);// take all the complement cells  and take the boundary
	vector<Cell> collision_cells;

	//cout<<"safe_point_to_dilate stop 1 "<<endl;
	//test this step
	select_preserved_cells(dt,preserved_cells,erode,v);
	// end test
	//cout<<"safe_point_to_dilate stop 2 "<<endl;

	unvisited(dt,v,collision_cells);// it is ok

	for(unsigned int i = 0; i < collision_cells.size(); ++i){
		if( dt.is_infinite(collision_cells[i] ) ){
			cout<<"is infinite in safe"<<endl;
			getchar();
		}		
	}
	//cout<<"safe_point_to_dilate stop 3 "<<endl;

	//test for collision
	//vector<Cell> cell_out_star;
	check_cells_in_conflict(dt,collision_cells,erode[0]->info().m_label,v,collision,cell_out_star);
	//end test 

	//cout<<"safe_point_to_dilate stop 4 "<<endl;

	/*vector<Point> cellpoints;
	if(draw){
		cellpoints.push_back( collision );
		cellpoints.push_back( v->point() );
		drawMesh(dt,cellpoints);
		getchar();
	}*/

	//if( check_preserved_cells2(dt,preserved_cells,collision,v,draw) ) { cout<<"preserved: "<<endl;}

	if( !check_sphere_collision2(dt,collision_cells,collision,v,draw) ) { /*cout<<"collision: "<<endl;*/ return false;}
	if( !check_sphere_collision2(dt,cell_out_star,collision,v,draw) ) { /*cout<<"collision2 : "<<endl;*/ return false;}
	
	//cout<<"safe_point_to_dilate stop 5"<<endl;

	/*if(draw){
		cellpoints.push_back( collision );
		drawMesh(dt,cellpoints);
		getchar();
	}*/

	if( !check_preserved_cells2(dt,preserved_cells,collision,v,draw) ) { /*cout<<"not preserved: "<<endl;*/ return false;}
	//cout<<"safe_point_to_dilate stop 6 "<<endl;
	return true;
}

void hash_centroids_path(Delaunay &dt, Vertex &v, vector<Cell> &keep, map<Point,int> &pointLabel, map<Point,int> &pointOriginalLabel,Point p){// cells doesn't have to be hashed
	list<Cell> incidentCells;
	dt.finite_incident_cells(v,back_inserter(incidentCells));
	Vertex v1,v2,v3;
	for(std::list<Cell>::iterator it = incidentCells.begin(); it != incidentCells.end(); it++){ 
		Cell c = *it;// if c is in cells continue
		//if ( std::find(keep.begin(), keep.end(), c) != keep.end() )	continue; // c is in cells
		int idx = c->index(v);
		v1 = c->vertex( dt.vertex_triple_index(idx,0) );
		v2 = c->vertex( dt.vertex_triple_index(idx,1) );
		v3 = c->vertex( dt.vertex_triple_index(idx,2) );
		Point cen = CGAL::centroid(v1->point(),v2->point(),v3->point(),p);
		decrease_decimals(cen);
		pointLabel[cen] = c->info().m_label;
		pointOriginalLabel[cen] = c->info().m_original_label;
	}
}

//hash 
void hash_cells_in_conflict(Delaunay &dt, Vertex &v,Point &p,map<Point,int> &pointLabel,vector<Cell> &path){
	list<Cell> incidentCells;
	dt.finite_incident_cells(v,back_inserter(incidentCells));
	Point v1,v2,v3;
	for(std::list<Cell>::iterator it = incidentCells.begin(); it != incidentCells.end(); it++){ 
		Cell c = *it;
		if ( std::find(path.begin(), path.end(), c) != path.end() )	continue; // c is in path
		if( dt.side_of_sphere(c,p,true) == CGAL::ON_BOUNDED_SIDE  ){			
			int vindex = c->index(v);
			for( unsigned int i = 0; i < 3;++i){					
				int index  = dt.vertex_triple_index(vindex,i);
				v1 = c->vertex( dt.vertex_triple_index(index,0) )->point();
				v2 = c->vertex( dt.vertex_triple_index(index,1) )->point();
				v3 = c->vertex( dt.vertex_triple_index(index,2) )->point();
				Point cen = CGAL::centroid( v1, v2, v3, p );
				decrease_decimals(cen);
				pointLabel[cen] = c->info().m_label;
			}
		}
	}
}

void hash_cells_out_star(Delaunay &dt, Vertex &v,Point &p,map<Point,int> &pointLabel,vector<Cell> &cells_out_star){
	Point v1,v2,v3;
	for(unsigned int i = 0; i < cells_out_star.size(); ++i){
		Cell c = cells_out_star[i];
		for( unsigned int j = 0; j < 4;++j){	
			int index = j;
			v1 = c->vertex( dt.vertex_triple_index(index,0) )->point();
			v2 = c->vertex( dt.vertex_triple_index(index,1) )->point();
			v3 = c->vertex( dt.vertex_triple_index(index,2) )->point();
			Point cen = CGAL::centroid( v1, v2, v3, p );
			decrease_decimals(cen);
			pointLabel[cen] = c->info().m_label;
		}
	}
}

void relabel_path_with_hash(Delaunay &dt,Vertex &v,Vertex &w, int label, map<Point,int> &pointLabel, map<Point,int> &pointOriginalLabel, deque<Vertex> &singularVertices, deque<Cell> &holes){
	deque<Cell> withoutLabel;
	list<Cell> incidentCellsW;
	dt.finite_incident_cells(w,back_inserter(incidentCellsW));
	for(std::list<Cell>::iterator it = incidentCellsW.begin(); it != incidentCellsW.end(); it++){
		Cell c = *it;
		Point cenc = CGAL::centroid(dt.tetrahedron(c));
		decrease_decimals(cenc);
		int mapLabel = new_map(pointLabel,cenc);
		if( mapLabel != -1){
			if( mapLabel == HOLE)	hole_cell(c,holes);
			else	c->info().m_label = mapLabel;
			c->info().m_original_label = pointOriginalLabel[cenc];
		}
		else
			withoutLabel.push_back(c);
	} 
	if(withoutLabel.size()!=0){
		//relabel_with_neighbor(dt,withoutLabel,w,v,label,holes);
		for( unsigned int i = 0; i < withoutLabel.size(); ++i){
			Cell c = withoutLabel[i];
			c->info().m_label = label;
		}
	}
}

bool insert_to_dilate(Delaunay &dt, deque<Vertex> &singularVertices, vector<vector<Cell>> &cellSets,Vertex &v,int label,deque<Cell> &holes,bool draw ){

	Point p;
	vector<Cell> path;
	if( !path_to_join(dt,cellSets,v,path) ){
		cout<<"not path found"<<endl;
		return false;
	}
	int pathLabel = path[0]->info().m_label;

	/*cout<<"size path "<<path.size()<<endl;
	getchar();
	vector<Cell> reset;
	drawCells(dt,reset);*/

	vector<Cell> complement = complement_cells(dt,path,v);// the cells of the same label in the star: complement

	vector<Cell> cells_out_star;
	if( safe_point_to_dilate(dt,path,complement,v,p,cells_out_star,draw) ){ // this is ok
		map<Point,int> pointLabel;
		map<Point,int> pointOriginaLabel;
		hash_centroids_path(dt,v,complement,pointLabel,pointOriginaLabel,p);	
		hash_cells_in_conflict(dt,v,p,pointLabel,path);
		
		// hash cell out of star
		hash_cells_out_star(dt,v,p,pointLabel,cells_out_star);

		Vertex w = dt.insert(p, path[0]);	
		relabel_with_hash_2(dt,v,w,pathLabel,pointLabel,pointOriginaLabel,singularVertices,holes);
		//relabel_path_with_hash(dt, v, w, label, pointLabel,pointOriginaLabel, singularVertices, holes); /// HERE we need just one label ot the keep and erode set
		detect_new_singularities_2(dt,w,singularVertices);
		update_vertex_info(dt, v );

		//cout<<"size sing_regions: "<<v->info().m_singular_regions_and_subgraphs.size()<<endl;
		if(draw){
			//cout<<"we are here? "<<v->point()<<" * "<<w->point()<<endl;
			drawMesh(dt);
			aroundVertex(dt,v,w,label);
			getchar();
		}
		return true;
	}	
	else{
		//cout<<"not possible inserted to dilate"<<endl;
		return false;
	}
}

bool insert_point_2_region(Delaunay &dt, deque<Vertex> &singularVertices, vector<vector<Cell>> &cellSets,Vertex &v,int label,int min,int max,deque<Cell> &holes,bool draw ){
	Point p;
	vector<Cell> cells_out_star;
	cout<<" begin this"<<endl;
	if( safe_point2(dt,cellSets[max],cellSets[min],v,p,cells_out_star,draw) ){
		map<Point,int> pointLabel;
		map<Point,int> pointOriginaLabel;
		hash_centroids_to_labels_2(dt,v,cellSets[min],pointLabel,pointOriginaLabel,p);	

		hash_cells_out_star(dt,v,p,pointLabel,cells_out_star);

		Vertex w = dt.insert(p, cellSets[max][0]);	
		relabel_with_hash_2(dt,v,w,cellSets[min][0]->info().m_label,pointLabel,pointOriginaLabel,singularVertices,holes);
		detect_new_singularities_2(dt,w,singularVertices);
		update_vertex_info(dt, v );
		if(draw){
			aroundVertex(dt,v,w,label);
			getchar();
		}
		return true;
	}	
	else if( safe_point2(dt,cellSets[min],cellSets[max],v,p,cells_out_star,draw) ){
		map<Point,int> pointLabel;
		map<Point,int> pointOriginaLabel;
		hash_centroids_to_labels_2(dt,v,cellSets[max],pointLabel,pointOriginaLabel,p);	

		hash_cells_out_star(dt,v,p,pointLabel,cells_out_star);

		Vertex w = dt.insert(p, cellSets[min][0]);	
		relabel_with_hash_2(dt,v,w,cellSets[max][0]->info().m_label,pointLabel,pointOriginaLabel,singularVertices,holes);
		detect_new_singularities_2(dt,w,singularVertices);
		update_vertex_info( dt, v );
		if(draw){
			aroundVertex(dt,v,w,label);
			getchar();
		}
		return true;
	}
	else{
		//cout<<"not possible inserted"<<endl;
		return false;
	}
}

//flag to decide repairing 1 or 2
bool insert_to_erode(Delaunay &dt, deque<Vertex> &singularVertices, vector<vector<Cell>> &cellSets,Vertex &v,int label,int min,int max,deque<Cell> &holes,bool draw,bool flag = false){
	if( v->info().has_two_regions() ){
		//cout<<"insert 2 regions"<<endl;
		return insert_point_2(dt,singularVertices,cellSets,v,label,min,max,holes,draw,flag);
	}
	else{
		//cout<<"insert more than 2 regions"<<endl;
		return insert_point_2(dt,singularVertices,cellSets,v,label,min,max,holes,draw,flag);
	}
}