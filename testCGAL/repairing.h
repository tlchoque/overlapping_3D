
void mark_boundary(Delaunay &dt){
	Vertex v = dt.infinite_vertex();
	list<Vertex> incident;
	dt.finite_incident_vertices(v,back_inserter(incident));
	for(std::list<Vertex>::iterator it = incident.begin(); it != incident.end(); it++)
		(*it)->info().m_type = 0;
}

void boundary_vertices(Delaunay &dt,vector<Vertex> & boundaryVertices){
	boundaryVertices.clear();
	Vertex v1;
	for(Finite_facets_iterator fi=dt.finite_facets_begin();fi!=dt.finite_facets_end(); fi++){
		if( fi->first->info().m_label == fi->first->neighbor(fi->second)->info().m_label )		continue;// include facets adjacents to infinite cells
		if( dt.is_infinite(fi->first) || dt.is_infinite(fi->first->neighbor(fi->second)) )		continue;
		for(unsigned int i = 0; i < 3; ++i){
			v1 = fi->first->vertex( dt.vertex_triple_index(fi->second,i) );
			if(v1->info().m_state == 0){
				//update_regions_around_vertex( dt, v1 );
				update_vertex_info(dt,v1);
				v1->info().m_state = 1;
				boundaryVertices.push_back(v1);
			}
		}
	}
}

void singular_vertices(Delaunay &dt,vector<Vertex> & boundaryVertices, deque<Vertex> & singularVertices){
	singularVertices.clear();
	for(unsigned int i = 0; i < boundaryVertices.size(); ++i){
		//if( boundaryVertices[i]->info().is_singular() ){
		if( boundaryVertices[i]->info().is_singular_2() ){
			if( boundaryVertices[i]->info().has_two_regions() )
				singularVertices.push_back(boundaryVertices[i]);
			else
				singularVertices.push_front(boundaryVertices[i]);
			boundaryVertices[i]->info().m_state = 2;// singular vertex
		}
	}
}

//for new relabeling

bool is_label_around(Delaunay &dt, Vertex &v, int label){
	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		if( (*it)->info().m_label == label )
			return true;
	}
	return false;
}

bool add_tetrahedron(Delaunay &dt, Cell &c, int label){
	if( c->info().m_label == label) return false;
	Cell neighborCell;
	vector<int> facets;
	for(unsigned int i = 0; i < 4; ++i){
		neighborCell = c->neighbor(i);
		if(neighborCell->info().m_label == label)
			facets.push_back(i);
	}

	if( facets.size() == 0 ){
		for(unsigned int i = 0; i < 4; ++i){
			if( is_label_around(dt,c->vertex(i),label ) )
				return false;
		}
		return true;
	}		
	else if( facets.size() == 1 )
		return !is_label_around(dt,c->vertex(facets[0]),label);
	
	else if( facets.size() == 2 ){// i am checking the self tetra
		Cell_circulator cc_start = dt.incident_cells( c,facets[0],facets[1] );
		Cell_circulator cc = cc_start;	
		//cout<<"f = 2"<<endl;
		do{
			if(cc->info().m_label == label)
				return false;
			++cc;
		} while (cc!= cc_start);
		return true;
	}	
	else	return true;
}

bool is_singular_after_relabeling(Delaunay &dt, Vertex &v){
	vector<int> labels;
	list<Cell> incident;
	bool cen = false;
	dt.finite_incident_cells(v,back_inserter(incident));	
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		if( (*it)->info().m_visited ) continue;
		if( (*it)->info().m_label == -1 ) continue;
		int cellLabel =  (*it)->info().m_label;
		if( std::find(labels.begin(), labels.end(), cellLabel ) != labels.end() ){//found
			cen = true;
			break;
		}
		else{
			labels.push_back( cellLabel );
			visite_cells(dt, (*it) ,v , cellLabel );
		}
	}
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++)
		(*it)->info().m_visited = false;
	return cen;
}

void neighbor_vertices(Delaunay &dt, vector<Cell> &cells,vector<Vertex> &vertices,Vertex &v) {
	Vertex v1;
	for(unsigned int i = 0; i < cells.size(); ++i){
		int idx = cells[i]->index(v);
		for(unsigned int j = 0; j < 3; ++j){
			v1 = cells[i]->vertex( dt.vertex_triple_index(idx,j) );
			if( v1->info().m_auxiliar == false ){
				v1->info().m_auxiliar = true;
				vertices.push_back(v1);
			}
		}
	}
	for(unsigned int i = 0; i < vertices.size(); ++i)
		vertices[i]->info().m_auxiliar = false;
}

bool check_neighbor_vertices( Delaunay &dt, vector<Vertex> &vertices ){
	for(unsigned int i = 0; i < vertices.size(); ++i ){
		if( vertices[i]->info().has_more_than_two_regions() || vertices[i]->info().is_corner() ){
			if( is_singular_after_relabeling(dt, vertices[i] ) )
				return false;
		}
	}
	return true;
}

bool relabel_if_possible(Delaunay & dt, vector<Cell> &cells, int label,Vertex v){
	deque<Cell> cpyCells(cells.begin(),cells.end());
	size_t cen = cpyCells.size();
	Cell c;
	int old = cells[0]->info().m_label;

	/*cout<<endl<<"size cpycells: "<<cpyCells.size()<<endl;
	if(v->point() == Point(5,347,94)){
		vector<Point> points;
		points.push_back(v->point() );
		points.push_back(Point(3,432,5) );
		drawMesh(dt);
		vertexRegion(dt,v,points,v->info().m_regions_around[ v->info().m_regions_around[0] == label  ]);
		getchar();
	}*/

	while(cpyCells.size() != 0 && cen !=0){		
		c = cpyCells.back();
		cpyCells.pop_back();
		if( add_tetrahedron(dt,c,label) ){
			c->info().m_label = label;
			c->info().m_state = 1;
			cen = cpyCells.size();	
		}
		else{
			cpyCells.push_front(c);
			--cen;
		}
	}

	vector<Vertex> neighbors;
	neighbor_vertices(dt,cells,neighbors,v);

	if( cpyCells.size() != 0  ){
		for(unsigned int i = 0; i < cells.size(); ++i){
			cells[i]->info().m_label = old;
			cells[i]->info().m_state = 0;
		}
		return false;
	}
	else{		
		if( !check_neighbor_vertices(dt,neighbors) ){
			for(unsigned int i = 0; i < cells.size(); ++i){
				cells[i]->info().m_label = old;
				cells[i]->info().m_state = 0;
			}
			return false;
		}
		else{
			for(unsigned int i = 0; i < neighbors.size(); ++i ){
				update_regions_around_vertex( dt,neighbors[i] );
			}
			update_regions_around_vertex( dt, v );
			return true;
		}		
	}	
}

//for fix singularities

void get_neighbor_cells(Delaunay &dt,Cell &c, Vertex &v,int &label,vector<Cell> &cells){
	vector<Cell> visitedCells;
	visitedCells.push_back(c);
	c->info().m_visited = true;
	cells.push_back(c);

	Cell vc,c0;
	while(visitedCells.size() != 0){
		vc = visitedCells.back();
		visitedCells.pop_back();
		int index = vc->index(v);
		for(unsigned int i = 0; i < 3; ++i){
			c0 = vc->neighbor( dt.vertex_triple_index( index, i ) );
			if(!c0->info().m_visited && c0->info().m_label == label ){
				visitedCells.push_back(c0);
				c0->info().m_visited = true;
				cells.push_back(c0);
			}
		}
	}
}

void set_cells_around_singular_vertex(Delaunay &dt, Vertex &v, int label, vector<vector<Cell>> &cellSets){
	cellSets.clear();
	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		if( (*it)->info().m_visited ) {
			(*it)->info().m_visited = false;
			continue;
		}
		if( (*it)->info().m_label == label ){
			vector<Cell> cells;
			get_neighbor_cells(dt, (*it) , v, label, cells);
			cellSets.push_back(cells);
		}
		(*it)->info().m_visited = false;
	}
}

bool check_cells_state(vector<Cell> &cells){
	for(unsigned int i = 0; i < cells.size(); ++i){
		if( cells[i]->info().m_state )
			return false;
	}
	return true;
}

void relabel_cells_1(Delaunay &dt, deque<Vertex> &singularVertices, vector<Cell> &cells,Vertex &v, int label){
	int newlabel = v->info().m_regions_around[0];
	if( newlabel == label )	newlabel = v->info().m_regions_around[1];
	if( newlabel == -1 ) return;

	for(unsigned int i = 0; i < cells.size();++i){
		Cell c = cells[i];
		c->info().m_label = newlabel;
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
	update_regions_around_vertex( dt, v );
}

void relabel_cells_2(Delaunay &dt, deque<Vertex> &singularVertices, vector<Cell> &cells,Vertex &v, int label){
	for(unsigned int i = 0; i < cells.size();++i){
		Cell c = cells[i];
		double a2 = 0;int lbl = 0;int idx=0; // there will always be a facet 
		for(unsigned int j = 0; j < 3;++j){
			idx = dt.vertex_triple_index( c->index(v),j );
			Facet f = Facet(c, idx);
			double squaredAreaTriangle = dt.triangle( f ).squared_area();
			if( squaredAreaTriangle > a2 && c->neighbor( idx )->info().m_label != label && c->neighbor( idx )->info().m_label != -1) {
				a2 = squaredAreaTriangle;
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
	update_regions_around_vertex( dt, v );
}

void relabeling_if_checked( Delaunay &dt, deque<Vertex> &singularVertices,vector<vector<Cell>> &cellSets,Vertex &v,int label ){
	vector<double> ratio(cellSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxVolume=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		data.push_back( make_pair( larger_edge(dt,cellSets[l]) , volume_cell_set(cellSets[l]) )  );
		if(data[l].first > maxEdge)		maxEdge = data[l].first ;
		if(data[l].second > maxVolume)		maxVolume = data[l].second;
	}
	int max=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		ratio[l] = 0.7*data[l].first/maxEdge + 0.3*data[l].second/maxVolume; 
		if( ratio[max] < ratio[l])  max = l;
	}

	if(v->info().has_more_than_two_regions() ){
		for(unsigned int l = 0; l < cellSets.size(); ++l){
			if( l != max && check_cells_state(cellSets[l]) ){
				relabel_cells_2(dt,singularVertices,cellSets[l],v,label);//update all affected vertices
				cellSets.erase(cellSets.begin()+l);
				--l;
				--max;
			}
		}
	}
	else{
		for(unsigned int l = 0; l < cellSets.size(); ++l){
			if( l != max && check_cells_state(cellSets[l]) ){
				relabel_cells_1(dt,singularVertices,cellSets[l],v,label);//update all affected vertices
				cellSets.erase(cellSets.begin()+l);
				--l;
				--max;
			}
		}
	}
}

void relabeling( Delaunay &dt, deque<Vertex> &singularVertices,vector<vector<Cell>> &cellSets,Vertex &v,int label ){
	vector<double> ratio(cellSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxVolume=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		data.push_back( make_pair( larger_edge(dt,cellSets[l]) , volume_cell_set(cellSets[l]) )  );
		if(data[l].first > maxEdge)		maxEdge = data[l].first ;
		if(data[l].second > maxVolume)		maxVolume = data[l].second;
	}
	int max=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		ratio[l] = 0.7*data[l].first/maxEdge + 0.3*data[l].second/maxVolume; 
		if( ratio[max] < ratio[l])  max = l;
	}

	for(unsigned int l = 0; l < cellSets.size(); ++l){
		//if( l != max ){
			relabel_cells_2(dt,singularVertices,cellSets[l],v,label);//update all affected vertices
			//cellSets.erase(cellSets.begin()+l);
		//}
	}
}

void relabeling_decrease_set_to_two( Delaunay &dt, deque<Vertex> &singularVertices,vector<vector<Cell>> &cellSets,Vertex &v,int label ){
	vector<double> ratio(cellSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxVolume=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		data.push_back( make_pair( larger_edge(dt,cellSets[l]) , volume_cell_set(cellSets[l]) )  );
		if(data[l].first > maxEdge)		maxEdge = data[l].first ;
		if(data[l].second > maxVolume)		maxVolume = data[l].second;
	}
	for(unsigned int l = 0; l < ratio.size(); ++l){
		ratio[l] = 0.7*data[l].first/maxEdge + 0.3*data[l].second/maxVolume; 
	}
	int v1=0,v2=1;
	if(ratio[v2] > ratio[v1]){
		v1 = 1;
		v2 = 0;
	}
	for(unsigned int l = 2; l < ratio.size(); ++l){// choose the 2 greater volumes v1 > v2
		if(ratio[l] > ratio[v1]){
			int aux = v1;
			v1 = l;
			v2 = aux;
		}
		else if(ratio[l] > ratio[v2])   v2 = l;
	}					
	for(unsigned int l = 0; l < cellSets.size(); ++l){
		if( l != v1 && l!=v2){
			relabel_cells_2(dt,singularVertices,cellSets[l],v,label);
			cellSets.erase(cellSets.begin()+l);
		}
	}
}

bool is_delaunay_after_displacement(Delaunay &dt, Vertex v, const Point &p){
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
	// are incident bi-cells Delaunay?
	std::vector<Facet> facets;
	facets.reserve(128);
	dt.incident_facets(v, std::back_inserter(facets));
	size = facets.size();
	for(std::size_t i=0; i<size; i++){
		const Facet &f = facets[i];
		Cell c = f.first;
		int j = f.second;
		Cell cj = c->neighbor(j);
		int mj = dt.mirror_index(c, j);
		Vertex h1 = c->vertex(j);
		if(dt.is_infinite(h1)) {
			if( dt.side_of_sphere(c, cj->vertex(mj)->point(), true) != CGAL::ON_UNBOUNDED_SIDE) {
				v->set_point(ant);
				return false;
			}
		} else {
			if(dt.side_of_sphere(cj, h1->point(), true) != CGAL::ON_UNBOUNDED_SIDE) {
				v->set_point(ant);
				return false;
			}
		}
	}
	v->set_point(ant);
	return true;
}

bool insert_point(Delaunay &dt, deque<Vertex> &singularVertices, vector<vector<Cell>> &cellSets,Vertex &v,int label,bool draw ){
	deque<Cell> holes;

	//cout<<"inserting"<<endl;
	vector<double> ratio(cellSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxVolume=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		data.push_back( make_pair( larger_edge(dt,cellSets[l]) , volume_cell_set(cellSets[l]) )  );
		if(data[l].first > maxEdge)		maxEdge = data[l].first ;
		if(data[l].second > maxVolume)		maxVolume = data[l].second;
	}
	for(unsigned int l = 0; l < ratio.size(); ++l)
		ratio[l] = 0.7*data[l].first/maxEdge + 0.3*data[l].second/maxVolume; 
	
	int max=0,min=1;
	if(ratio[min] > ratio[max]){
		max = 1;
		min = 0;
	}

	if(!v->info().is_corner() ){
		Point pv1 = average_neighbor_boundary_vertices(dt,cellSets[max],v);		
		int step1 = getDisplacedPoint(dt,v,pv1);
		if(step1){//exist space to move
			if( is_delaunay_after_displacement(dt,v,pv1 ) ) { 
				//cout<<"booty moved"<<endl;
				//v->set_point(pv1);
			}											
		}
	}		
	Point pv2 = average_neighbor_boundary_vertices(dt,cellSets[min],v);
	int step2 = getDisplacedPoint(dt,v,pv2);
	if(step2){
		map<Point,int> pointLabel;
		hash_centroids_to_labels(dt,v,cellSets[max],pointLabel,pv2);	
		Vertex w = dt.insert(pv2);	
		check_if_cells_set_is_out(dt,cellSets[min],v,label,singularVertices);
		
		relabel_with_hash(dt,v,w,label,pointLabel,singularVertices,holes);
		detect_new_singularities(dt,w,singularVertices);
		update_regions_around_vertex(dt, v );

		if(draw){
			aroundVertex(dt,v,w,label);
			getchar();
		}
		return true;
	}
	else{
		cout<<"not possible inserted"<<endl;
		return false;
	}
}

void restart_vertex_state(Delaunay &dt){
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) { 
		vi->info().m_state = 0;
		//vi->info().m_labels.clear();
	} 
}

//for curvature optimization

void get_boundary_facets(Delaunay & dt,vector<Cell> &cells,vector<Facet> &facets,Vertex &v,int label){
	Cell vc,c;
	for(unsigned int i = 0; i < cells.size();++i){
		vc = cells[i];
		int index = vc->index(v);
		for(unsigned int j = 0; j < 3;++j){
			int index2 = dt.vertex_triple_index(index,j); 
			c = vc->neighbor( index2 );
			//if( dt.is_infinite(c) ) continue;
			if(c->info().m_label != label ){
				//cout<<"opposed label "<<c->info().m_label<<endl;
				Facet f = Facet(vc,index2);
				facets.push_back(f);
			}
		}
	}
}

vector<Vertex> ring_vertices(Delaunay &dt, vector<Cell> &cells,Vertex &v,int label){
	map<Vertex,Vertex> vv;
	vector<Facet> facets;
	//cout<<"stop 0"<<endl;

	get_boundary_facets(dt,cells,facets,v,label);
	//cout<<"stop 1"<<endl;
	//cout<<"facets: "<<facets.size()<<endl;
	//if( v->info().is_corner()  ){
	//	//drawMesh(dt);
	//	drawCells(dt,cells);
	//	//getchar();
	//}

	/*if( v->point() == Point(5,776,479) ){
		drawMesh(dt);
		drawCells(dt,cells);
		getchar();
	}*/

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
	//cout<<"stop 2"<<endl;
	return orderedVertices;
}

void vertices_to_points(vector<Vertex> &vertices, vector<Point> &points){
	for(unsigned int i = 0; i< vertices.size(); ++i)
		points.push_back(vertices[i]->point());
}

//it is still an open problem
void is_joinable(Delaunay & dt,vector<vector<Cell>> &cellSets,Vertex &v,int label,bool draw){
	vector<Vertex> vertices;
	vector<Point> points;
	vertices = ring_vertices(dt,cellSets[0],v,label);
	vertices_to_points(vertices,points);
				
	drawMesh(dt);
	/*aroundVertex(dt,v,v,label);
	getchar();*/

	Vector normtangent = norm_tangent_plane(vertices,v);

	Point test = CGAL::ORIGIN + ( ( v->point() - CGAL::ORIGIN ) + 20*normtangent);
	points.push_back( test );
	vertexRegion(dt,v,points,label);
	getchar();
}

int new_label_2(Delaunay &dt, vector<Cell> &cells,Vertex &v,int label){
	Cell vc,c;
	vector<int> neighbor;
	vector<double> areas;
	for(unsigned int i = 0; i < cells.size(); ++i ){
		vc = cells[i];
		int index = vc->index(v);
		for(unsigned int j = 0; j < 3;++j){
			int index2 = dt.vertex_triple_index(index,j); 
			c = vc->neighbor( index2 );
			if( c->info().m_label != label && c->info().m_label != -1 ){
				Facet f = Facet(vc,index2);
				double squaredAreaTriangle = dt.triangle( f ).squared_area();
				vector<int>::iterator ite = std::find(neighbor.begin(), neighbor.end(), c->info().m_label );
				if( ite != neighbor.end() ){// found
					int position = std::distance( neighbor.begin(), ite );										
					areas[position] += squaredAreaTriangle;
				}				
				else{// not found
					neighbor.push_back( c->info().m_label );
					areas.push_back(squaredAreaTriangle);
				}
			}
		}
	}
	//if( areas.size() == 0) return HOLE;
	int poslabel = std::distance( areas.begin() , max_element( areas.begin(),areas.end() ) );
	return neighbor[poslabel];
}

int new_label_for_relabeling(Delaunay &dt, vector<Cell> &cells,Vertex &v,int label){
	if(v->info().has_more_than_two_regions() )
		return new_label_2(dt,cells,v,label);
	else
		return v->info().m_regions_around[ v->info().m_regions_around[0] == label ];
}

bool relabeling_without_making_singularities( Delaunay &dt, deque<Vertex> &singularVertices,vector<vector<Cell>> &cellSets,Vertex &v,int label ){
	vector<double> ratio(cellSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxVolume=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		data.push_back( make_pair( larger_edge(dt,cellSets[l]) , volume_cell_set(cellSets[l]) )  );
		if(data[l].first > maxEdge)		maxEdge = data[l].first ;
		if(data[l].second > maxVolume)		maxVolume = data[l].second;
	}
	int max=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		ratio[l] = 0.7*data[l].first/maxEdge + 0.3*data[l].second/maxVolume; 
		if( ratio[max] < ratio[l])  max = l;
	}

	bool cen = true;
	for(unsigned int l = 0; l < cellSets.size(); ++l){
		//if( l != max && check_cells_state(cellSets[l]) ){
		if( l != max ){
			int newlabel = new_label_for_relabeling(dt,cellSets[l],v,label);
			if( !relabel_if_possible(dt,cellSets[l],newlabel,v ) ){
				cen = false;
			}
			else{
				cellSets.erase(cellSets.begin()+l);
				--l;
				--max;
			}
		}
	}
	return cen;
}

void fix_singularities(Delaunay &dt){
	mark_boundary(dt);
	vector<Vertex> boundaryVertices;
	deque<Vertex> singularVertices;		
	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);
	cout<<"number of non manifold vertices before "<<singularVertices.size()<<endl;
	vector<Point> points;

	bool draw = 0;
	while(singularVertices.size() != 0 ){
		Vertex v = singularVertices.front();
		singularVertices.pop_front();

		//cout<<"size singularities: "<<singularVertices.size()<<endl;
		//if(singularVertices.size() == 2592)	draw = 1;

		while(v->info().m_singular_regions_around.size() != 0){
			int label = v->info().m_singular_regions_around.back();
			v->info().m_singular_regions_around.pop_back();

			vector<vector<Cell>> cellSets;
			set_cells_around_singular_vertex(dt,v,label,cellSets);

			//is_joinable(dt,cellSets,v,label,draw);

			if(draw){				
				drawMesh(dt);
				aroundVertex(dt,v,v,label);
				getchar();
			}
			
			relabeling_if_checked(dt,singularVertices,cellSets,v,label );
			if(cellSets.size() == 1) {
				if(draw){
					/*vector<Vertex> vertices;
					vector<Point> points;
					vertices = ring_vertices(dt,cellSets[0],v,label);
					vertices_to_points(vertices,points);
					vertexRegion(dt,v,points,label);
					getchar();*/

					aroundVertex(dt,v,v,label);
					getchar();
				}
				continue;
			}

			if(cellSets.size() > 2)
				relabeling_decrease_set_to_two(dt,singularVertices,cellSets,v,label );

			if( cellSets.size() == 2){
				if( !insert_point(dt,singularVertices,cellSets,v,label,draw ) ){
					//cout<<"relabeling without check"<<endl;
					relabeling(dt,singularVertices,cellSets,v,label);
					if(draw){
						aroundVertex(dt,v,v,label);
						getchar();
					}
					//break;
				}					
			}
			//update_regions_around_vertex(dt,v);
		}
		v->info().m_state=1;
	}

	restart_vertex_state(dt);
	boundaryVertices.clear();
	singularVertices.clear();
	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);
	cout<<"number of non manifold vertices after "<<singularVertices.size()<<endl;
	
	//for smoothing
	restart_vertex_state(dt);
}

void fix_singularities2(Delaunay &dt){
	mark_boundary(dt);
	vector<Vertex> boundaryVertices;
	deque<Vertex> singularVertices;		
	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);
	cout<<"number of non manifold vertices before "<<singularVertices.size()<<endl;
	vector<Point> points;


	vector<Point> vtkPoints;
	for(unsigned int i = 0 ; i < singularVertices.size();++i)
		vtkPoints.push_back(singularVertices[i]->point());

	/*savePoints2(vtkPoints,"carp/original/new_singularities.vtk");
	getchar();*/

	bool draw = 1;
	while(singularVertices.size() != 0 ){
		Vertex v = singularVertices.front();
		singularVertices.pop_front();
		//cout<<"size singularities: "<<singularVertices.size()<<endl;

		while(v->info().m_singular_regions_around.size() != 0){
			int label = v->info().m_singular_regions_around.back();
			v->info().m_singular_regions_around.pop_back();
			
			vector<vector<Cell>> cellSets;
			set_cells_around_singular_vertex(dt,v,label,cellSets);

			if(draw){				
				drawMesh(dt);
				aroundVertex(dt,v,v,label);
				getchar();
			}

			if(cellSets.size() > 2)
				relabeling_decrease_set_to_two(dt,singularVertices,cellSets,v,label );

			if( !relabeling_without_making_singularities(dt,singularVertices,cellSets,v,label) ){
				if( !insert_point(dt,singularVertices,cellSets,v,label,draw ) ){
					break;
				}	
			}
			if(draw){
				aroundVertex(dt,v,v,label);
				getchar();
			}
		}
		v->info().m_state=1;
	}

	restart_vertex_state(dt);
	boundaryVertices.clear();
	singularVertices.clear();
	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);
	cout<<"number of non manifold vertices after "<<singularVertices.size()<<endl;

	//for smoothing
	restart_vertex_state(dt);
}