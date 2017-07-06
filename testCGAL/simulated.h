
double volume_cell_set_2(deque<Cell> &cells){
	double vol=0;
	for(unsigned int i = 0; i < cells.size(); ++i)
		vol+= CGAL::volume( cells[i]->vertex(0)->point(),cells[i]->vertex(1)->point(),cells[i]->vertex(2)->point(),cells[i]->vertex(3)->point() ) ;
	return vol;	
}

bool erode_without_making_singularities_sa( Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool draw, deque<Cell> &holes ){	
	int min,max;
	bool change = false;
	sort_by_criteria(dt,cellSets,v,min,max,change);
	int newlabel = new_label_for_relabeling(dt,cellSets[min],v,label); // change because it should not be HOLE

	int index;
	double r = ((double) rand() / (RAND_MAX));
	if( r < 0.5) index = min;
	else index = max;
	

	if( newlabel == HOLE )	{
		make_hole(dt,cellSets[index], singularities,v,holes);
		return false;
	}
	if( !relabel_if_possible_2(dt,cellSets[index],newlabel,v ) ){ // test if we canerode change to MIN
		//if( !relabel_if_possible_2(dt,cellSets[min],newlabel,v ) ){
		if( !insert_point_2(dt,singularities,cellSets,v,label,min,max,holes,draw) ){
			//cout<<"make hole"<<endl;
			make_hole(dt,cellSets[index], singularities,v,holes); //change to MIN
			return false;
		}
	}
	return true;
}

bool repair_sa(Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool draw, deque<Cell> &holes){
	if( v->info().has_two_regions() ){
		if( !is_erodible(dt,cellSets[0],v,label) ){
			Vertex w;
			int opposite_label = v->info().m_regions_around[ v->info().m_regions_around[0] == label ];
			if( extreme_vertex_of_singular_edge(dt,v,w,opposite_label ) ){
				vector<vector<Cell>> cellSets_2;
				set_cells_around_singular_vertex(dt,w,opposite_label,cellSets_2);// get the other singular regions
				erode_without_making_singularities_sa(dt, cellSets_2,singularities,w,opposite_label,draw,holes);
			}
		}
	}
	set_cells_around_singular_vertex(dt,v,label,cellSets);
	if( cellSets.size() == 1 )
		return true;
	else
		return erode_without_making_singularities_sa(dt, cellSets,singularities,v,label,draw,holes);	
}

void solve_singularities_with_holes(Delaunay &dt,deque<Vertex> &singularVertices,deque<Cell> &holes){
	bool draw = 0;
	int count = 1;
	while(singularVertices.size() != 0 ){
		Vertex v = singularVertices.front();
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
				//cout<<"stop 1"<<endl;
				repair_sa(dt,cellSets,singularVertices,v,label,draw,holes);
			}
			else if( cellSets.size() == 1 ){
				//cout<<"stop 2"<<endl;
				repair_2(dt,cellSets,singularVertices,v,label,draw,holes);
			}
			else{
				//cout<<"stop 3"<<endl;
				repair_3(dt,cellSets,singularVertices,v,label,draw,holes);
			}
			if(draw){
				aroundVertex(dt,v,v,label);
				getchar();
			}
		}
		v->info().m_state=1;
	}
}

double vol_cell(Cell &c){
	return CGAL::volume( c->vertex(0)->point(),c->vertex(1)->point(),c->vertex(2)->point(),c->vertex(3)->point() );
}

int label_for_cell(Delaunay &dt, Cell &c){
	vector<int> neighbor;
	vector<double> areas;
	Cell vc;
	for(unsigned int j = 0; j < 4;++j){
		vc = c->neighbor( j );
		if( dt.is_infinite(vc) ) continue;
		if( vc->info().m_label != -1 && vc->info().m_label != HOLE){// if neighbor label is differen of the c original label
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

void set_stars_vertex(Delaunay &dt, Vertex &v){
	vector<vector<Vertex>> &stars = v->info().m_sa_stars;
	stars.clear();
	map<int,int> label_to_index;
	for(unsigned int i = 0; i < v->info().m_regions_around.size();++i){
		vector<Vertex> star;
		stars.push_back(star);
		label_to_index[ v->info().m_regions_around[i] ] = i;
	}

	//cout<<"stop 2"<<endl;
	vector<vector<Cell>> cellsSet( stars.size() );	
	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));	
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		//cout<<"label "<<(*it)->info().m_label<<endl;
		cellsSet[ label_to_index[ (*it)->info().m_label ]  ].push_back(*it);
	}
	//cout<<"stop 3"<<endl;

	if( v->info().m_regions_around.size() == 2 &&  !v->info().is_corner() ){
		//cout<<"option 0"<<endl;
		vector<Vertex> r1 = ring_vertices( dt,cellsSet[0],v,v->info().m_regions_around[0] );
		stars[0] = r1;
		reverse(r1.begin(),r1.end());
		stars[1] = r1;
	}
	else{
		//cout<<"option 1"<<endl;
		for(unsigned int i = 0; i < cellsSet.size(); ++i){
			if( v->info().m_regions_around[i] == HOLE ) continue;
			stars[i] = ring_vertices( dt,cellsSet[i],v,v->info().m_regions_around[i] );
			/*print(stars[i]);
			vector<Point> ptos = vtop(stars[i]);

			drawMesh(dt,ptos);
			getchar();*/
		}
	}
}

double mesh_curvature(Delaunay &dt){
	double total = 0;
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) { 
		Vertex v = vi;
		if( !v->info().is_boundary() ) continue;
		set_stars_vertex(dt,v);
		vector<double> curvatures;
		v->info().m_curvature_sa = vertex_curvature(v,curvatures);
		total+= v->info().m_curvature_sa;
		/*cout<<v->info().m_curvature_sa <<endl;

		vector<Point> lala;
		lala.push_back(v->point());
		drawMesh(dt,lala);
		getchar();*/
	}
	return total;
}

double cost(Delaunay &dt,deque<Cell> &holes){
	double c1 = volume_cell_set_2(holes);
	//double c2 = mesh_curvature(dt);
	//return c1 + c2;
	//return c1;
	return holes.size();
}

void neighbor(Delaunay &dt,deque<Cell> &new_holes){
	//we relabel the hole cells until it wouldn´t be possible
	new_holes = get_holes(dt);
	/*double volume = vol_cell( new_holes[0] );
	unsigned int index = 0;
	index = rand()%new_holes.size();
	cout<<"index "<<index<<endl;*/
	Cell c;
	int new_label;
	//new_label = label_for_cell(dt,c);
	//double r = ((double) rand() / (RAND_MAX));
	//for( unsigned int i = 0; i < new_holes.size(); ++i){
	while(true){
		unsigned int index = rand()%new_holes.size();
		c = new_holes[index];
		if( c->info().m_original_label != HOLE){
			new_label = c->info().m_original_label;
			break;
		}
		else{
			int lbl = label_for_cell(dt,c);
			if(lbl == -1) {continue;}
			new_label = lbl;
			break;
		}
	}
	cout<<"new_label "<<new_label<<endl;
	c->info().m_label = new_label;

	deque<Vertex> singular_vertices;
	for( unsigned int i = 0; i < 4; ++i){
		Vertex v = c->vertex(i);
		update_vertex_info(dt,v);
		if( v->info().is_singular_2() ){
			singular_vertices.push_back(v);
			v->info().m_state = 2;
		}
	}
	cout<<"before solving "<<new_label<<endl;
	solve_singularities_with_holes(dt,singular_vertices,new_holes);
	select_hole_cells(dt,new_holes);
	label_hole_cells_without_making_singularities(dt,new_holes);
}

void neighbor_2(Delaunay &dt,deque<Cell> &new_holes){
	new_holes = get_holes(dt);
	deque<Vertex> singularities;
	deque<Cell> hls = new_holes;

	//cout<<"start while"<<endl;
	//for( unsigned int i = 0; i < new_holes.size(); ++i){
	while(hls.size() != 0){

		//cout<<"size "<<hls.size()<<endl;
		Cell c = hls.front();
		hls.pop_front();
		int new_label = label_for_cell(dt,c);;
		if( new_label == -1){
			//cout<<"new_label "<<new_label<<endl;
			hls.push_back(c);
			continue;
		}
		c->info().m_label = new_label;

		//c->info().m_label = c->info().m_original_label;
		for( unsigned int j = 0; j < 4; ++j){
			Vertex v = c->vertex(j);
			update_vertex_info(dt,v);
			if( v->info().is_singular_2() && v->info().m_state != 2 ){
				singularities.push_back(v);
				v->info().m_state = 2;
			}
		}
	}
	solve_singularities_with_holes(dt,singularities,new_holes);
	select_hole_cells(dt,new_holes);
	label_hole_cells_without_making_singularities(dt,new_holes);

	/*vector<Vertex> boundaryVertices;
	deque<Vertex> singularVertices;	
	restart_vertex_state(dt);
	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);
	cout<<"number of non manifold vertices simulated "<<singularVertices.size()<<endl;
	cout<<"number of vertices simulated "<<dt.number_of_vertices()<<endl;*/
}

double acceptance_probability(double old_cost,double new_cost, double T ){
	return exp(old_cost - new_cost)/T;
}

double random_sa(){
	return (double) (rand()/ (double) RAND_MAX);
}

void simulated_annealing(Delaunay &dt,deque<Cell> holes){
	if(holes.size() == 0) return;
	cout<<"begin"<<endl;
	double old_cost = cost(dt,holes);
	double T = 1.0;
	double T_min = 0.00001;
    double alpha = 0.9;
	while ( T > T_min ){
        int i = 1;
        while (i <= 100){
			deque<Cell> new_holes;
			Delaunay new_dt = dt;
			neighbor_2(new_dt,new_holes);
			if(new_holes.size() == 0) {dt = new_dt;T_min =2; break;}
            double new_cost = cost(new_dt,new_holes);
            double ap = acceptance_probability(old_cost, new_cost, T);
			cout<<"old and new cost "<< volume_cell_set_2(holes)<<"  "<<volume_cell_set_2(new_holes)<<endl;
			cout<<"old and new size holes "<<holes.size()<<"  "<<new_holes.size()<<endl;
			drawMesh(dt);
			getchar();
            if( ap > random_sa() ){
				cout<<"accepted"<<endl;
                dt = new_dt;			
				holes.swap(new_holes);
                old_cost = new_cost;				
			}
			drawMesh(dt);
			getchar();
            i += 1;			
		}
        T = T*alpha;
	}

	vector<Vertex> boundaryVertices;
	deque<Vertex> singularVertices;	
	restart_vertex_state(dt);
	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);
	cout<<"number of non manifold vertices simulated "<<singularVertices.size()<<endl;
	cout<<"number of vertices simulated "<<dt.number_of_vertices()<<endl;
}