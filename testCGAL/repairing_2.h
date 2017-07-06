
bool is_erodible2(Delaunay &dt, vector<Cell> &cells,Vertex v, int label){
	vector<Vertex> vertices;
	vertices = ring_vertices(dt,cells,v,label);
	Vector vcur = normal_curvature(vertices,v);
	Vector vec = Vector(0,0,0);
	for(unsigned int i = 0; i < vertices.size(); ++i){
		Vertex w = vertices[i];
		if( vertices[i]->info().is_singular_2() ) continue;
		vector<vector<Cell>> cellSets;
		set_cells_around_singular_vertex(dt,w,label,cellSets);
		vector<Vertex> ring;
		ring = ring_vertices(dt,cellSets[0],w,label);
		vec = vec + normal_curvature(ring,w);
	}
	/*cout<<"mean curvature "<<mean_curvature(vertices,v )<<endl;
	cout<<"sum normal: "<<vec<<"\n"<<"v curature: "<<vcur<<endl;


	vector<Point> draw;
	draw.push_back(v->point() );
	draw.push_back( CGAL::ORIGIN + ( v->point() - CGAL::ORIGIN + 1000*vcur ) ); 
	draw.push_back( CGAL::ORIGIN + ( v->point() - CGAL::ORIGIN + 1000*vec ) ); 
	drawMesh(dt,draw);
	getchar();*/


	if(vec.z() > 0)// the vector points up
		return false;
	else 
		return true;
}

bool is_erodible(Delaunay &dt, vector<Cell> &cells,Vertex v, int label){
	vector<Vertex> vertices;
	vertices = ring_vertices(dt,cells,v,label);
	/*cout<<"stop 0 is_erodible"<<endl;
	cout<<"size vertices "<<vertices.size()<<endl;
	if( v->info().is_corner() ) cout<<"is corner"<<endl;*/
	/*vector<Point> pts = vtop(vertices);
	pts.push_back(v->point());
	drawMesh(dt,pts);
	getchar();*/

	//print(vertices);
	double mcv = mean_curvature2(vertices,v);
	//cout<<"stop 1 is_erodible"<<endl;
	double sum = 0;
	for(unsigned int i = 0; i < vertices.size(); ++i){
		Vertex w = vertices[i];
		if( vertices[i]->info().is_singular_2() ) continue;
		vector<vector<Cell>> cellSets;
		set_cells_around_singular_vertex(dt,w,label,cellSets);
		vector<Vertex> ring;
		ring = ring_vertices(dt,cellSets[0],w,label);
		sum = sum + mean_curvature2(ring,w);
	}
	//cout<<"vertex and sum: "<<mcv<<" * "<<sum<<endl;

	/*vector<Point> draw;
	draw.push_back(v->point() );
	draw.push_back( CGAL::ORIGIN + ( v->point() - CGAL::ORIGIN + 1000*vcur ) ); 
	draw.push_back( CGAL::ORIGIN + ( v->point() - CGAL::ORIGIN + 1000*vec ) ); 
	drawMesh(dt,draw);
	getchar();*/

	//cout<<"stop 2 is_erodible"<<endl;
	if( mcv*sum > 0)// the vector points up
		return false;
	return true;
}

bool produce_singularities( Delaunay &dt, vector<Vertex> &vertices ){
	for(unsigned int i = 0; i < vertices.size(); ++i ){
		//cout<<"verification  "<<vertices[i]->point()<<endl;
		//if( vertices[i]->info().has_more_than_two_regions() || vertices[i]->info().is_corner() ){
			if( vertices[i]->info().m_state != 2   )
				if( is_singular_after_relabeling_2(dt, vertices[i]) )
					return true;
		//}
	}
	return false;
}

void label_hole_cells_without_making_singularities(Delaunay &dt, deque<Cell> &holes){
	for( unsigned int i = 0; i < holes.size(); ++i ){
		Cell c = holes[i];
		holes[i]->info().m_hole_visited = false;
		if( c->info().m_label != HOLE){ // if the cell is not a hole
			holes.erase(holes.begin() + i);
			--i;
			continue;
		}
		vector<Vertex> neighbors;
		for( unsigned int j = 0; j < 4; ++j){
			neighbors.push_back(c->vertex(j) );			
		}

		c->info().m_label = c->info().m_original_label;
		if( c->info().m_original_label == HOLE)		continue;// if the new label is hole

		if( produce_singularities(dt,neighbors) )
			c->info().m_label = HOLE;
		else{
			holes.erase(holes.begin() + i);
			--i;
			for( unsigned int j = 0; j < 4; ++j)
				update_vertex_info( dt,neighbors[j] );
		}
	}
}

void label_boundary_cells(Delaunay &dt, vector<Cell> &holes, int other_label){
	for( unsigned int i = 0; i < holes.size(); ++i ){
		Cell c = holes[i];
		int original = c->info().m_label;
		c->info().m_label = other_label;
		vector<Vertex> neighbors;
		for( unsigned int j = 0; j < 4; ++j)
			neighbors.push_back(c->vertex(j) );		
		if( produce_singularities(dt,neighbors) )
			c->info().m_label = original;
		else{
			for( unsigned int j = 0; j < 4; ++j)
				update_vertex_info( dt,neighbors[j] );
		}
	}
}

bool produce_singularities_2( Delaunay &dt, vector<Vertex> &vertices ){
	for(unsigned int i = 0; i < vertices.size(); ++i ){
		if( is_singular_after_relabeling_2(dt, vertices[i]) )
			return true;
	}
	return false;
}

void boundary_cells(Delaunay & dt,vector<Cell> &cells,Vertex v,vector<Cell> &boundary){
	int original = cells[0]->info().m_label;
	for(unsigned int i = 0; i < cells.size(); ++i ){
		Cell c = cells[i];
		int vindex = c->index(v);
		for(unsigned int j = 0; j < 3; ++j ){
			int idx = dt.vertex_triple_index(vindex,j);
			int nlabel = c->neighbor(idx)->info().m_label;
			if( nlabel != original ){
				boundary.push_back(c);
				break;
			}
		}		
	}
}

bool relabel_if_possible_2(Delaunay & dt, vector<Cell> &cells, int other_label,Vertex v){
	vector<Vertex> neighbors;
	neighbor_vertices(dt,cells,neighbors,v);

	int old_label = cells[0]->info().m_label; 
	for(unsigned int i = 0; i < cells.size(); ++i )
		cells[i]->info().m_label = other_label;

	if( produce_singularities(dt,neighbors) ){
		for(unsigned int i = 0; i < cells.size(); ++i )
			cells[i]->info().m_label = old_label;
		/*vector<Cell> boundaries;
		boundary_cells(dt,cells,v,boundaries);
		label_boundary_cells(dt,boundaries,other_label);*/		
		return false;
	}
	else{
		for(unsigned int i = 0; i < neighbors.size(); ++i ){
			update_vertex_info( dt,neighbors[i] );
		}
		update_vertex_info( dt, v );
		return true;
	}
}

bool relabel_if_possible_3(Delaunay & dt, vector<Cell> &cells, int other_label,Vertex v, deque<Vertex> &singularities){
	if( !check_cells_state(cells) ){
		//cout<<"inserting point "<<endl;
		return false;
		//return relabel_if_possible_2(dt,cells,other_label,v);		
	}
	//cout<<"other label "<<other_label<<endl;

	int old_label = cells[0]->info().m_label; 
	for(unsigned int i = 0; i < cells.size(); ++i ){
		cells[i]->info().m_label = other_label;
		cells[i]->info().m_state = 1;	
	}

	vector<Vertex> neighbors;
	neighbor_vertices(dt,cells,neighbors,v);
	for(unsigned int i = 0; i < neighbors.size(); ++i ){
		Vertex w = neighbors[i];
		update_vertex_info( dt, w );
		if(  w->info().m_state!=2 && w->info().is_singular_2()  ){
			w->info().m_state = 2;
			if( w->info().has_two_regions() )
				singularities.push_back(w);
			else
				singularities.push_front(w);
		}
	}
	update_vertex_info( dt, v );
	return true;	
}

bool relabel_without_checking(Delaunay & dt, vector<Cell> &cells, int other_label,Vertex v, deque<Vertex> &singularities){
	//if( !check_cells_state(cells) ){
	//	//cout<<"inserting point "<<endl;
	//	return false;
	//}

	int old_label = cells[0]->info().m_label; 
	for(unsigned int i = 0; i < cells.size(); ++i ){
		cells[i]->info().m_label = other_label;
		cells[i]->info().m_state = 1;	
	}

	vector<Vertex> neighbors;
	neighbor_vertices(dt,cells,neighbors,v);
	for(unsigned int i = 0; i < neighbors.size(); ++i ){
		update_vertex_info( dt,neighbors[i] );
		if(  neighbors[i]->info().m_state!=2 && neighbors[i]->info().is_singular_2()  ){
			neighbors[i]->info().m_state = 2;
			singularities.push_back(neighbors[i] );
		}
	}
	update_vertex_info( dt, v );
	return true;	
}

//for relabeling path in one component

int adyacent_cell_label(Delaunay &dt, Cell &c, int alabel, int blabel){
	vector<int> neighbor;
	vector<double> areas;
	Cell vc;
	for(unsigned int j = 0; j < 4;++j){
		vc = c->neighbor( j );
		if( dt.is_infinite(vc) ) continue;
		if( vc->info().m_label == alabel || vc->info().m_label == blabel){// if neighbor label is equl to alabel or blabel
			Facet f = Facet(vc,j);
			double squaredAreaTriangle = dt.triangle( f ).squared_area();
			vector<int>::iterator ite = std::find(neighbor.begin(), neighbor.end(), vc->info().m_label ); // locate iterator in neighbor
			if( ite != neighbor.end() ){// found
				int position = std::distance( neighbor.begin(), ite );										
				areas[position] += squaredAreaTriangle;
			}				
			else{// not found
				neighbor.push_back( vc->info().m_label );
				areas.push_back(squaredAreaTriangle);
			}
		}
	}
	if( areas.size() != 0){
		int poslabel = std::distance( areas.begin() , max_element( areas.begin(),areas.end() ) );
		return neighbor[poslabel];
	}
	return -1;
}

bool relabel_path(Delaunay & dt, vector<Cell> &cells, int alabel,int blabel,Vertex v, deque<Vertex> &singularities){
	if( !check_cells_state(cells) ){
		//cout<<"fails checked "<<endl;
		return false;
	}

	deque<Cell> queue;
	for( unsigned int i = 0; i< cells.size(); ++i){
		//cout<<"label: "<<cells[i]->info().m_label<<endl;
		queue.push_back(cells[i]);	
	}

	while( queue.size() != 0 ){
		Cell c = queue.front();
		queue.pop_front();
		int lbl = adyacent_cell_label(dt,c,alabel,blabel);
		if( lbl != -1){
			c->info().m_label = lbl;
			//state of the cell
			c->info().m_state = true;
		}
		else{
			queue.push_back(c);
		}
	}
	vector<Vertex> restartIndex;
	int k = 1;
	for(unsigned int i = 0; i < cells.size(); ++i){
		Cell c = cells[i];
		int vindex = c->index(v);//changed
		for( unsigned int j = 0; j < 3; ++j){			
			int index = dt.vertex_triple_index(vindex,j);
			Vertex w = c->vertex(index);
			if(w->info().m_index == -1){
				w->info().m_index = k;
				restartIndex.push_back(w);

				update_vertex_info(dt,w);
				if( w->info().is_singular_2() && w->info().m_state!= 2)
					singularities.push_back(w);
			}
		}
	}
	update_vertex_info(dt,v);
	if( v->info().is_singular_2() && v->info().m_state!= 2)
		singularities.push_back(v);

	for(unsigned int i = 0 ; i < restartIndex.size(); ++i)
		restartIndex[i]->info().m_index=-1;

	//cout<<"exito relabeled cells: "<<cells.size()<<endl;
	return true;
}

void make_hole(Delaunay &dt, vector<Cell> &cells, deque<Vertex> &singularities,Vertex &v, deque<Cell> &holes){
	for(unsigned int i = 0; i < cells.size(); ++i ){
		Cell c = cells[i];
		hole_cell(c,holes);
	}
	vector<Vertex> neighbors;
	neighbor_vertices(dt,cells,neighbors,v);
	for(unsigned int i = 0; i < neighbors.size(); ++i ){
		update_vertex_info( dt,neighbors[i] );
		if( neighbors[i]->info().is_singular_2() && neighbors[i]->info().m_state != 2 ){
			singularities.push_back( neighbors[i] );
			neighbors[i]->info().m_state = 2;
		}
	}
	update_vertex_info( dt, v );
}

bool erode_without_making_singularities( Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw, deque<Cell> &holes ){	
	int min,max;
	bool change = false;
	sort_by_criteria(dt,cellSets,v,min,max,change);
	int newlabel = new_label_for_relabeling(dt,cellSets[max],v,label); // change because it should not be HOLE
	if( newlabel == HOLE )	{
		make_hole(dt,cellSets[max], singularities,v,holes);
		return false;
	}
	if( !relabel_if_possible_2(dt,cellSets[max],newlabel,v ) ){ // test if we canerode
		//if( !relabel_if_possible_2(dt,cellSets[min],newlabel,v ) ){
			if( !insert_point_2(dt,singularities,cellSets,v,label,min,max,holes,draw) ){
				//cout<<"make hole"<<endl;
				make_hole(dt,cellSets[max], singularities,v,holes);
				return false;
			}
		//}
	}
	return true;
}

bool erode_test( Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw, deque<Cell> &holes ){	
	int min,max;
	bool change = false;
	sort_by_criteria(dt,cellSets,v,min,max,change);
	make_hole(dt,cellSets[max], singularities,v,holes);
	return false;
}

bool erode_with_checking( Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw, deque<Cell> &holes ){	
	if(cellSets.size() != 2 ) return false;
	int min,max;
	//cout<<"stop erode 0"<<endl;
	bool change = false;
	sort_by_criteria(dt,cellSets,v,min,max,change);
	
	int newlabel = new_label_for_relabeling(dt,cellSets[max],v,label); // change because it should not be HOLE
	if( !relabel_if_possible_3(dt,cellSets[max],newlabel,v,singularities ) ){ // test if we canerode
		//if( !relabel_if_possible_2(dt,cellSets[min],newlabel,v ) ){
		if(change){// try to relabel the other subgraph
			if( relabel_if_possible_3(dt,cellSets[min],newlabel,v,singularities) )
				return true;
		}
		if( !insert_point_2(dt,singularities,cellSets,v,label,min,max,holes,draw) ){
			return false;
		}
	}
	return true;
}

bool dilate_with_checking( Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw){	
	if(cellSets.size() != 2 ) return false;
	vector<Cell> path;
	path_to_join(dt,cellSets,v,path);

	/*drawMesh(dt);
	drawCells(dt,path);
	cout<<"drawing cells "<<path.size()<<endl;
	getchar();*/

	double s1 = criteria_1(dt,cellSets[0],v);
	double s2 = criteria_1(dt,cellSets[1],v);
	double s3 = criteria_1(dt,path,v);

	double s = s1;
	if( s2 > s1) s = s2;

	//cout<<"s3 y s: "<<s3<<" * "<<s<<" * "<<s/s3<<endl;
	if( s/s3 > 2 ){ 
		return false;
	}
	//cout<<"have passed "<<s3<<" * "<<s<<endl;
	/*if(draw){
		drawCells(dt,path);
		getchar();
	}*/
	return relabel_path(dt,path,cellSets[0][0]->info().m_label,cellSets[1][0]->info().m_label,v,singularities) ; 
}

void set_cells_by_region(Delaunay &dt, Vertex v,vector<vector<Cell>> &cells_by_label,int label){
	cells_by_label.clear();
	map<int,int> label_to_index;
	int count=0;
	for(unsigned int i = 0; i < v->info().m_regions_around.size(); ++i){
		if( v->info().m_regions_around[i] == HOLE || v->info().m_regions_around[i] == label) continue;
		label_to_index[ v->info().m_regions_around[i] ] = count;
		++count;
	}
	vector<vector<Cell>> cellsSet( count );	
	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));	
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){// there exist -1 label
		if( (*it)->info().m_label == HOLE || (*it)->info().m_label == label ) continue;
		cellsSet[ label_to_index[ (*it)->info().m_label ]  ].push_back(*it);
	}
	cells_by_label = cellsSet;
}

//for 2 subgraphs
bool repair_1(Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw, deque<Cell> &holes){
	if( v->info().has_two_regions() ){
		if( !is_erodible(dt,cellSets[0],v,label) ){
			Vertex w;
			int opposite_label = v->info().m_regions_around[ v->info().m_regions_around[0] == label ];
			if( extreme_vertex_of_singular_edge(dt,v,w,opposite_label ) ){
				vector<vector<Cell>> cellSets_2;
				set_cells_around_singular_vertex(dt,w,opposite_label,cellSets_2);// get the other singular regions
				erode_without_making_singularities(dt, cellSets_2,singularities,w,opposite_label,draw,holes);
				//erode_test(dt, cellSets_2,singularities,w,opposite_label,draw,holes);
			}
		}
	}
	set_cells_around_singular_vertex(dt,v,label,cellSets);
	if( cellSets.size() == 1 )
		return true;
	else{
		return erode_without_making_singularities(dt, cellSets,singularities,v,label,draw,holes);	
		//return erode_test(dt, cellSets,singularities,v,label,draw,holes);	
	}
}

bool repair_1_2(Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw, deque<Cell> &holes){
	if( v->info().has_two_regions() ){
		if( !is_erodible(dt,cellSets[0],v,label) ){
			Vertex w;
			int opposite_label = v->info().m_regions_around[ v->info().m_regions_around[0] == label ];
			if( extreme_vertex_of_singular_edge(dt,v,w,opposite_label ) ){
				//cout<<"points found"<<endl;
				vector<	Point> points;
				points.push_back(w->point() );
				if(draw){				
					drawMesh(dt,points);
					aroundVertex(dt,w,w,opposite_label);
					getchar();
				}

				vector<vector<Cell>> cellSets_2;
				set_cells_around_singular_vertex(dt,w,opposite_label,cellSets_2);// get the other singular regions
				if( erode_with_checking(dt, cellSets_2,singularities,w,opposite_label,draw,holes) ){
					if(draw){				
						//drawMesh(dt,points);
						aroundVertex(dt,w,w,opposite_label);
						getchar();
					}
				}				
			}
		}
	}
	set_cells_around_singular_vertex(dt,v,label,cellSets);
	if( cellSets.size() == 1 )
		return true;
	else{
		if( erode_with_checking(dt, cellSets,singularities,v,label,draw,holes) )
			return true;
		else
			 return dilate_with_checking(dt,cellSets,singularities,v,label,draw);
		//return erode_with_checking(dt, cellSets,singularities,v,label,draw,holes);	
	}
}

//for 1 subgraph
bool repair_2(Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw, deque<Cell> &holes){
	vector<vector<Cell>> cells_by_label;
	set_cells_by_region(dt,v,cells_by_label,label);

	vector<pair<vector<Cell> , double>> cells_criteria;
	for(unsigned int i = 0; i < cells_by_label.size(); ++i){
		double criteria = criteria_1(dt,cells_by_label[i],v);
		cells_criteria.push_back(make_pair(cells_by_label[i],criteria ) );
	}
	std::sort(cells_criteria.begin(), cells_criteria.end(), 
		boost::bind(&std::pair<vector<Cell> , double>::second, _1) > boost::bind(&std::pair<vector<Cell> , double>::second, _2));

	for(unsigned int i = 0; i < cells_criteria.size(); ++i){// we can try inserting points
	//for(unsigned int i = 0; i < 1; ++i){// we can try inserting points
		if(cells_criteria[i].first[0]->info().m_label == 0)	continue;

		int newlabel = new_label_for_relabeling(dt,cells_criteria[i].first,v,cells_criteria[i].first[0]->info().m_label);	
		if( newlabel != HOLE && relabel_if_possible_2(dt,cells_criteria[i].first,newlabel,v ) ){
		//if( relabel_without_checking(dt,cells_criteria[i].first,label,v,singularities ) ){
			if( !v->info().is_singular_2() )
				return true;
		}
		else{// try to insert a point
			if( cells_criteria.size() == 2){				
				vector<vector<Cell>> separate_set;
				separate_set.push_back(cells_criteria[i].first);
				int idx = (i == 0);
				separate_set.push_back(cells_criteria[idx].first);
				if( insert_point_2(dt,singularities,separate_set,v,label,0,1,holes,draw) )
					if( !v->info().is_singular_2() )
						return true;
			}
		}
	}
	make_hole(dt,cellSets[0],singularities,v,holes);
	return false;
}

bool repair_2_2(Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw, deque<Cell> &holes){
	if(draw){		
		drawMesh(dt);
		aroundVertex(dt,v,v,label);
		getchar();
	}

	vector<vector<Cell>> cells_by_label;
	set_cells_by_region(dt,v,cells_by_label,label);

	vector<pair<vector<Cell> , double>> cells_criteria;
	for(unsigned int i = 0; i < cells_by_label.size(); ++i){
		double criteria = criteria_1(dt,cells_by_label[i],v);
		cells_criteria.push_back(make_pair(cells_by_label[i],criteria ) );
	}

	std::sort(cells_criteria.begin(), cells_criteria.end(), 
		boost::bind(&std::pair<vector<Cell> , double>::second, _1) > boost::bind(&std::pair<vector<Cell> , double>::second, _2));

	//for(unsigned int i = 0; i < cells_criteria.size(); ++i){// for all the multimaterial meshes
	for(unsigned int i = 0; i < 1; ++i){// we can try inserting points
		
		/*if( cells_criteria.size() == 2){	
			vector<Cell> path;
			vector<vector<Cell>> separate_set;
			separate_set.push_back(cells_criteria[0].first);
			separate_set.push_back(cells_criteria[1].first);
			path_to_join(dt,separate_set,v,path);

			drawMesh(dt);
			drawCells(dt,path);
			getchar();

			if( relabel_path(dt,path,separate_set[0][0]->info().m_label,separate_set[1][0]->info().m_label,v,singularities) ){
				drawMesh(dt);
				drawCells(dt,path);
				cout<<"size path"<<path.size()<<endl;
				getchar();
				break;
			}			
		}*/

		//if(cells_criteria[i].first[0]->info().m_label == 0)	continue;
		int newlabel = new_label_for_relabeling(dt,cells_criteria[i].first,v,cells_criteria[i].first[0]->info().m_label);		
		if( relabel_if_possible_3(dt,cells_criteria[i].first,newlabel,v,singularities ) ){
			if( !v->info().is_singular_2() ){
				if(draw){		
					drawMesh(dt);
					aroundVertex(dt,v,v,label);
					getchar();
				}
				return true;
			}
		}
		else{
			if( cells_criteria.size() == 2){	
				vector<vector<Cell>> separate_set;
				separate_set.push_back(cells_criteria[i].first);
				separate_set.push_back(cells_criteria[i == 0].first);

				if( insert_point_2(dt,singularities,separate_set,v,label,0,1,holes,draw) ){
					if(draw){		
						drawMesh(dt);
						aroundVertex(dt,v,v,label);
						getchar();
					}
					//cout<<"no pos guau"<<endl;
					return true;
				}
				else{// find the path 
					return dilate_with_checking(dt,separate_set,singularities,v,label,draw);
				}
			}
		}
	}
	if(draw){		
		drawMesh(dt);
		aroundVertex(dt,v,v,label);
		getchar();
		draw = 0;
	}
	cout<<"one problem"<<endl;
	return false;
}

//for more than 2 graphs
void repair_3(Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw, deque<Cell> &holes){
	vector<double> ratio(cellSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxVolume=0;

	//cout<<"<top 0 repair_3"<<endl;
	for(unsigned int i = 0; i < ratio.size(); ++i){
		ratio[i] = criteria_1(dt,cellSets[i],v);
	}
	int v1=0,v2=1;
	if(ratio[v2] < ratio[v1]){
		v1 = 1;
		v2 = 0;
	}
	for(unsigned int l = 2; l < ratio.size(); ++l){// choose the 2 greater volumes v1 < v2
		if(ratio[l] < ratio[v1]){
			int aux = v1;
			v1 = l;
			v2 = aux;
		}
		else if(ratio[l] < ratio[v2])   v2 = l;
	}
	for(unsigned int l = 0; l < cellSets.size(); ++l){
		if( l != v1 && l!=v2){
			int newlabel = new_label_for_relabeling(dt,cellSets[l],v,label);

			if( newlabel == HOLE) {
				make_hole(dt,cellSets[l],singularities,v,holes);
			}
			else if(! relabel_if_possible_2(dt,cellSets[l],newlabel,v ) ){
				make_hole(dt,cellSets[l],singularities,v,holes);
			}
			cellSets.erase(cellSets.begin()+l);
			--l;
			--v1;--v2;
		}
	}
	repair_1(dt,cellSets,singularities,v,label,draw,holes);
}

bool repair_3_2(Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw, deque<Cell> &holes){
	if(draw){		
		cout<<"size: "<<cellSets.size()<<endl;
		drawMesh(dt);
		aroundVertex(dt,v,v,label);
		getchar();
	}

	vector<double> ratio(cellSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxVolume=0;
	for(unsigned int i = 0; i < ratio.size(); ++i){
		ratio[i] = criteria_1(dt,cellSets[i],v);
	}
	int v1=0,v2=1;
	if(ratio[v2] < ratio[v1]){
		v1 = 1;
		v2 = 0;
	}
	for(unsigned int l = 2; l < ratio.size(); ++l){// choose the 2 greater volumes v1 < v2
		if(ratio[l] < ratio[v1]){
			int aux = v1;
			v1 = l;
			v2 = aux;
		}
		else if(ratio[l] < ratio[v2])   v2 = l;
	}
	for(unsigned int l = 0; l < cellSets.size(); ++l){		
		if( l != v1 && l!=v2){
			int newlabel = new_label_for_relabeling(dt,cellSets[l],v,label);
			if(! relabel_if_possible_3(dt,cellSets[l],newlabel,v,singularities ) ){
				continue;
			}
			cellSets.erase(cellSets.begin()+l);
			--l;
			--v1;--v2;
		}
	}

	set_cells_around_singular_vertex(dt,v,label,cellSets);
	if(draw){		
		cout<<"size: "<<cellSets.size()<<endl;
		drawMesh(dt);
		aroundVertex(dt,v,v,label);
		getchar();
		// = 0;
	}
	return repair_1_2(dt,cellSets,singularities,v,label,draw,holes);
}

// for simulated annealing

void copy_label(Delaunay &dt){
	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		Cell c = it;
		c->info().m_old_solution = c->info().m_label;
	}
}

void select_hole_cells(Delaunay &dt,deque<Cell> &holes){
	for(unsigned int i = 0 ; i < holes.size(); ++i){
		if( holes[i]->info().m_hole_visited ||  !dt.is_cell(holes[i]) ){// for repeated cells and visited for this array
			holes.erase(holes.begin() + i);
			--i;
			continue;
		}
		holes[i]->info().m_hole_visited = true;
		if( holes[i]->info().m_label != HOLE  ){
			holes[i]->info().m_hole_visited = false;
			holes.erase(holes.begin() + i);
			--i;
		}		
	}
}

deque<Cell> get_holes(Delaunay &dt){
	deque<Cell> holes;
	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		if( it->info().m_label == HOLE){
			holes.push_back(it);
		}
	}
	return holes;
}

double testing(Delaunay &dt, Vertex &v){
	list<Cell> incidentCellsW;
	dt.finite_incident_cells(v,back_inserter(incidentCellsW));
	double count=0;
	for(std::list<Cell>::iterator it = incidentCellsW.begin(); it != incidentCellsW.end(); it++){
		if( (*it)->info().m_label == HOLE)
			++count;
	}
	return count;
}

void prepare_for_sa(Delaunay &dt,deque<Vertex> &singularVertices,deque<Cell> &holes){
	bool draw = 1;
	int count = 1;
	Vertex testing;
	while(singularVertices.size() != 0 ){
		Vertex v = singularVertices.front();
		//if( v->point() == Point(882,886,65) ) {draw = 1;testing = v;}
		singularVertices.pop_front();
		vector<pair<int,int>> &sing_regions = v->info().m_singular_regions_and_subgraphs;
		while( sing_regions.size() != 0){
			int label= sing_regions.front().first;
			vector<vector<Cell>> cellSets;
			set_cells_around_singular_vertex(dt,v,label,cellSets);
			if(draw){				
				drawMesh(dt);
				aroundVertex(dt,v,v,label);
				getchar();
			}
			if( cellSets.size() == 2 ){
				//cout<<"stop 2"<<endl;
				repair_1(dt,cellSets,singularVertices,v,label,draw,holes);
			}
			else if( cellSets.size() == 1 ){
				//cout<<"stop 3 "<<endl;
				repair_2(dt,cellSets,singularVertices,v,label,draw,holes);
			}
			else{
				//draw=1;
				//cout<<"stop 4 "<<endl;
				repair_3(dt,cellSets,singularVertices,v,label,draw,holes);
			}
			if(draw){
				aroundVertex(dt,v,v,label);
				getchar();
			}
		}
		v->info().m_state=1;
	}

	//cout<<"<top 5"<<endl;
	select_hole_cells(dt,holes);

	/*sort( holes.begin(), holes.end() );
	holes.erase( unique( holes.begin(), holes.end() ), holes.end() );*/

	cout<<"selected holes from deque: "<<holes.size();

	//deque<Cell> new_holes = get_holes(dt);
	//cout<<"real holes: "<<new_holes.size();
	label_hole_cells_without_making_singularities(dt,holes);
}

void set_label_for_sa(Delaunay &dt){// set an original label for all the cells
	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		Cell c = it;
		c->info().m_old_solution = c->info().m_label;
		c->info().m_new_solution = c->info().m_label;
	}
}

void repair_singularities(Delaunay &dt,deque<Cell> &holes){
	vector<Vertex> boundaryVertices;
	deque<Vertex> singularVertices;		

	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);
	cout<<"number of non manifold vertices before "<<singularVertices.size()<<endl;
	cout<<"number of vertices "<<dt.number_of_vertices()<<endl;
	bool draw = 0;	
	int count = 0;

	//vector<Point> vtkPoints;
	//for(unsigned int i = 0 ; i < singularVertices.size();++i){
	//	vtkPoints.push_back( singularVertices[i]->point() );
	//	//cout<<singularVertices[i]->point()<<endl;
	//}
	//savePoints(vtkPoints,"test1/original/old_singularities.vtk");

	while(singularVertices.size() != 0 ){		
		++count;

		//if( count == 1617) draw = 1;
		//if( v->point() == Point(882,886,65) ) draw = 1;
		//cout<<"count: "<<count<<endl;

		Vertex v = singularVertices.front();
		singularVertices.pop_front();
		v->info().m_state = 1;				
		vector<pair<int,int>> &sing_regions = v->info().m_singular_regions_and_subgraphs;
		while( sing_regions.size() != 0){
			int label = sing_regions.front().first;
			//cout<<"label "<<label<<endl;
			vector<vector<Cell>> cellSets;
			set_cells_around_singular_vertex(dt,v,label,cellSets);
			if(draw){				
				drawMesh(dt);
				aroundVertex(dt,v,v,label);
				getchar();
			}

			/*if( cellSets.size() == 2 ){
				if( insert_to_dilate(dt,singularVertices,cellSets,v,label,holes,draw) )
					continue;
			}*/
			
			if( cellSets.size() == 2 ){
				//cout<<"stop 0"<<endl;
				if( !repair_1_2(dt,cellSets,singularVertices,v,label,draw,holes) ){
					cout<<"stop 1"<<endl;
					v->info().m_state=2;
					break;
				}
				//break;
			}
			else if( cellSets.size() == 1 ){
				//draw = 1;
				//cout<<"stop 1"<<endl;
				if( !repair_2_2(dt,cellSets,singularVertices,v,label,draw,holes) ){
					cout<<"stop 2"<<endl;
					v->info().m_state=2;
					break;
				}
			}
			else{
				//cout<<"stop 2"<<endl;
				if( !repair_3_2(dt,cellSets,singularVertices,v,label,draw,holes) ){
					cout<<"stop 3"<<endl;
					v->info().m_state=2;
					break;
				}
			}
			if(draw){
				aroundVertex(dt,v,v,label);
				getchar();
			}
		}
		//v->info().m_state=1;
	}

	/*myfile1.close();
	myfile2.close();*/

	cout<<"holes "<<holes.size()<<endl;
	restart_vertex_state(dt);
	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);
	cout<<"number of non manifold vertices after "<<singularVertices.size()<<endl;
	cout<<"number of vertices "<<dt.number_of_vertices()<<endl;

	//prepare_for_sa(dt,singularVertices,holes);// make whole without singularities
	//cout<<"holes "<<holes.size()<<endl;
	/*restart_vertex_state(dt);
	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);
	cout<<"number of non manifold vertices after label holes "<<singularVertices.size()<<endl;
	cout<<"number of vertices "<<dt.number_of_vertices()<<endl;*/

	/*vector<Point> vtkPoints;
	for(unsigned int i = 0 ; i < singularVertices.size();++i){
		vtkPoints.push_back( singularVertices[i]->point() );
		cout<<singularVertices[i]->point()<<endl;
	}
	savePoints2(vtkPoints,"test1/original/singularities.vtk");*/

	/*drawMesh(dt,vtkPoints);
	getchar();*/
	
	//for smoothing
	restart_vertex_state(dt);
}


