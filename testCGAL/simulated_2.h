
bool relabel_for_simulated(Delaunay & dt, vector<Cell> &cells, int other_label,Vertex v, vector<Cell> &modified){
	//cout<<"relabel for sa stop 0"<<endl;
	vector<Vertex> neighbors;
	neighbor_vertices(dt,cells,neighbors,v);
	//cout<<"relabel for sa stop 0 1/4: "<<cells.size()<<endl;
	int old_label = cells[0]->info().m_label; 
	//cout<<"relabel for sa stop 0 1/2"<<endl;
	for(unsigned int i = 0; i < cells.size(); ++i )
		cells[i]->info().m_label = other_label;

	if( produce_singularities(dt,neighbors) ){
		//cout<<"relabel for sa stop 1"<<endl;
		for(unsigned int i = 0; i < cells.size(); ++i )
			cells[i]->info().m_label = old_label;
		return false;
	}
	else{
		//cout<<"relabel for sa stop 2"<<endl;
		for(unsigned int i = 0; i < neighbors.size(); ++i ){
			update_vertex_info( dt,neighbors[i] );
		}
		modified.insert(modified.end(), cells.begin(), cells.end());//append the c cells to modified cells
		update_vertex_info( dt, v );
		return true;
	}
}

bool dilate_sa( Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw,deque<Cell> &holes, vector<Cell> &modified){	
	if(cellSets.size() != 2 ) return false;
	//cout<<"dilate stop 0"<<endl;
	vector<Cell> path;
	if( !path_to_join(dt,cellSets,v,path) ) return false;
	//cout<<"dilate stop 1"<<endl;

	double s1 = criteria_1(dt,cellSets[0],v);
	double s2 = criteria_1(dt,cellSets[1],v);
	double s3 = criteria_1(dt,path,v);
	//cout<<"dilate stop 2"<<endl;

	double s = s1;
	if( s2 > s1) s = s2;
	//cout<<"s3 y s: "<<s3<<" * "<<s<<" * "<<s/s3<<endl;
	if( s/s3 > 2 ) return false;
	//return relabel_path(dt,path,cellSets[0][0]->info().m_label,cellSets[1][0]->info().m_label,v,singularities) ; 

	//cout<<"dilate stop 3"<<endl;
	if( !relabel_for_simulated(dt,path,cellSets[0][0]->info().m_label,v,modified ) ){ // test if we canerode_sa change to MIN
		return false;
	}
	return true;
}

bool erode_sa( Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool draw, deque<Cell> &holes, vector<Cell> &modified){	
	int min,max;
	bool change = false;
	sort_by_criteria(dt,cellSets,v,min,max,change);
	int newlabel = new_label_for_relabeling(dt,cellSets[min],v,label); // change because it should not be HOLE
	int index;

	double r = ((double) rand() / (RAND_MAX));
	if( r < 0.1) index = min;///3
	else index = max;

	if( newlabel == HOLE )	{
		make_hole(dt,cellSets[index], singularities,v,holes);
		return false;
	}

	double s = ((double) rand() / (RAND_MAX));
	if( s < 0.5) {
		if( !relabel_for_simulated(dt,cellSets[index],newlabel,v,modified ) ){ // test if we canerode_sa change to MIN
			//if( !insert_point_2(dt,singularities,cellSets,v,label,min,max,holes,draw) ){
				//cout<<"make hole"<<endl;
			if( !dilate_sa(dt,cellSets,singularities,v,label,draw,holes,modified) ){
				make_hole(dt,cellSets[index], singularities,v,holes); //change to MIN
				return false;
			}
		}
	}
	else{
		if( !relabel_for_simulated(dt,cellSets[index],newlabel,v,modified ) ){ // test if we canerode_sa change to MIN
			make_hole(dt,cellSets[index], singularities,v,holes); //change to MIN
			return false;
		}
	}

	//if( !relabel_for_simulated(dt,cellSets[index],newlabel,v,modified ) ){ // test if we canerode_sa change to MIN
	//	//if( !insert_point_2(dt,singularities,cellSets,v,label,min,max,holes,draw) ){
	//		//cout<<"make hole"<<endl;
	//	if( !dilate_sa(dt,cellSets,singularities,v,label,draw,holes,modified) ){
	//		make_hole(dt,cellSets[index], singularities,v,holes); //change to MIN
	//		return false;
	//	}
	//}
	return true;
}

bool solve_1(Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool draw, deque<Cell> &holes, vector<Cell> &modified){
	if( v->info().has_two_regions() ){// it has the HOLE label
		if( !is_erodible(dt,cellSets[0],v,label) ){
			Vertex w;
			int opposite_label = v->info().m_regions_around[ v->info().m_regions_around[0] == label ];
			if( extreme_vertex_of_singular_edge(dt,v,w,opposite_label ) ){
				vector<vector<Cell>> cellSets_2;
				set_cells_around_singular_vertex(dt,w,opposite_label,cellSets_2);// get the other singular regions
				erode_sa(dt, cellSets_2,singularities,w,opposite_label,draw,holes,modified);
			}
			else{
				dilate_sa(dt,cellSets,singularities,v,label,draw,holes,modified);
			}
		}
	}
	set_cells_around_singular_vertex(dt,v,label,cellSets);
	if( cellSets.size() == 1 )
		return true;
	else
		return erode_sa(dt, cellSets,singularities,v,label,draw,holes,modified);	
}

bool solve_2(Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool draw, deque<Cell> &holes, vector<Cell> &modified){
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
		if(cells_criteria[i].first[0]->info().m_label == 0)	continue;
		int newlabel = new_label_for_relabeling(dt,cells_criteria[i].first,v,cells_criteria[i].first[0]->info().m_label);	
		if( newlabel != HOLE && relabel_for_simulated(dt,cells_criteria[i].first,newlabel,v,modified ) ){
		//if( relabel_without_checking(dt,cells_criteria[i].first,label,v,singularities ) ){
			if( !v->info().is_singular_2() )
				return true;
		}
		else{// try to insert a point
			/*if( cells_criteria.size() == 2){				
				vector<vector<Cell>> separate_set;
				separate_set.push_back(cells_criteria[i].first);
				int idx = (i == 0);
				separate_set.push_back(cells_criteria[idx].first);
				if( insert_point_2(dt,singularities,separate_set,v,label,0,1,holes,draw) )
					if( !v->info().is_singular_2() )
						return true;
			}*/
		}
	}

	cout<<"not possible solve_2"<<endl;
	make_hole(dt,cellSets[0],singularities,v,holes);
	return false;
}

void walk_graph_exclude(Delaunay &dt,Cell &c, Vertex &v,int &label,vector<Cell> &set){	
	vector<Cell> visitedCells;
	visitedCells.push_back(c);
	set.push_back(c);
	c->info().m_visited = true;
	Cell vc,c0;
	while(visitedCells.size() != 0){
		vc = visitedCells.back();
		visitedCells.pop_back();
		int index = vc->index(v);
		for(unsigned int i = 0; i < 3;++i){
			c0 =  vc->neighbor( dt.vertex_triple_index(index,i) );
			if( !c0->info().m_visited && c0->info().m_label != label ){
				//cout<<"inserted: "<<visitedCells.size()<<endl;
				visitedCells.push_back(c0);
				set.push_back(c0);
				c0->info().m_visited = true;
			}
		}
	}
	//cout<<"visitedCells size: "<<visitedCells.size()<<endl;
}

vector<vector<Cell>> get_two_sets(Delaunay &dt, Vertex v, int label){
	vector<vector<Cell>> sets;
	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));	
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		Cell c = *it;
		if( c->info().m_label == label || c->info().m_visited ) continue;
		vector<Cell> set;
		walk_graph_exclude(dt,c,v,label,set);
		//cout<<" set size "<<set.size()<<endl;
		sets.push_back(set);
		
	}
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++)
		(*it)->info().m_visited = false;
	return sets;
}

vector<pair<vector<Cell> , double>> order_sets_by_criteria(Delaunay &dt, Vertex v, int label){
	vector<vector<Cell>> cells_by_label = get_two_sets(dt,v,label);
//	cout<<"size cells by label: "<<cells_by_label.size()<<endl;
	vector<pair<vector<Cell> , double>> cells_criteria;
	for(unsigned int i = 0; i < cells_by_label.size(); ++i){
		double criteria = criteria_1(dt,cells_by_label[i],v);
		cells_criteria.push_back(make_pair(cells_by_label[i],criteria ) );
	}

	std::sort(cells_criteria.begin(), cells_criteria.end(), 
		boost::bind(&std::pair<vector<Cell> , double>::second, _1) > boost::bind(&std::pair<vector<Cell> , double>::second, _2));
	return cells_criteria;
}

bool solve_2_2(Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool draw, deque<Cell> &holes, vector<Cell> &modified){
	vector<pair<vector<Cell> , double>> cells_criteria = order_sets_by_criteria(dt,v,label);
	//cout<<"size cells_criteria: "<<cells_criteria.size()<<endl;
	if( cells_criteria.size() == 2){
		//cout<<"case a"<<endl;
		vector<vector<Cell>> separate_set;
		separate_set.push_back(cells_criteria[1].first);
		separate_set.push_back(cells_criteria[0].first);
		//int newlabel = new_label_for_relabeling(dt,cells_criteria[0].first,v,cells_criteria[0].first[0]->info().m_label);
		double r = ((double) rand() / (RAND_MAX));
		if( r < 0.85){
			//cout<<"1st option "<<endl;
			if( relabel_for_simulated(dt,cells_criteria[0].first,label,v,modified ) ) return true;
		}
		else{
			//cout<<"2nd option "<<endl;
			if( dilate_sa(dt,separate_set,singularities,v,label,draw,holes,modified) ) return true;
		}

		/*if( !relabel_for_simulated(dt,cells_criteria[0].first,label,v,modified ) ){
			cout<<"stop 1"<<endl;
			if( dilate_sa(dt,separate_set,singularities,v,label,draw,holes,modified) ){
				cout<<"stop 2"<<endl;
				return true;
			}
		}
		else return true;*/
	}

	int newlabel = new_label_for_relabeling(dt,cellSets[0],v,label);
	if( relabel_for_simulated(dt,cellSets[0],newlabel,v,modified) ) return true;	

	//else if( cells_criteria.size() == 1 ){
	//	cout<<"case b"<<endl;
	//	if( relabel_for_simulated(dt,cells_criteria[0].first,label,v,modified ) ){
	//		return true;
	//	}
	//}
	//else if( cells_criteria.size() > 2 ){
	//	cout<<"case c: "<<cells_criteria.size() <<endl;
	//	vector< vector<Cell> > sets;
	//	sets.push_back( cells_criteria[1].first );
	//	sets.push_back( cells_criteria.back().first );
	//	if( dilate_sa(dt,sets,singularities,v,label,draw,holes,modified) ){
	//		cells_criteria = order_set_by_criteria(dt,v,label);
	//		if( !v->info().is_singular_2() ){
	//			return true;
	//		}
	//		else cells_criteria = order_set_by_criteria(dt,v,label);
	//	}
	
	//	for(unsigned int i = 0; i < cells_criteria.size(); ++i){// for all the multimaterial meshes
	//		int newlabel = new_label_for_relabeling(dt,cells_criteria[i].first,v,cells_criteria[i].first[0]->info().m_label);		
	//		if( relabel_for_simulated(dt,cells_criteria[i].first,newlabel,v,modified ) ){
	//			if( !v->info().is_singular_2() ){
	//				return true;
	//			}
	//		}
	//	}
	//}

	cout<<"not possible solve_2"<<endl;
	make_hole(dt,cellSets[0],singularities,v,holes);
	return false;
}

bool solve_3(Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool draw, deque<Cell> &holes, vector<Cell> &modified){
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

			if( newlabel == HOLE) {
				make_hole(dt,cellSets[l],singularities,v,holes);
			}
			else if(!relabel_for_simulated(dt,cellSets[l],newlabel,v,modified ) ){
				make_hole(dt,cellSets[l],singularities,v,holes);
			}
			cellSets.erase(cellSets.begin()+l);
			--l;
			--v1;--v2;
		}
	}
	return solve_1(dt,cellSets,singularities,v,label,draw,holes,modified);
}

void choose_holes(Delaunay &dt,deque<Cell> &holes){
	deque<Cell> aux = holes;
	for(unsigned int i = 0 ; i < holes.size(); ++i){
		Cell c = holes[i];
		if( c->info().m_label != HOLE || !dt.is_cell(c) || c->info().m_hole_visited ){//hole visited for repeated cells
			c->info().m_hole_state = false;//3 times repetead
			if( c->info().m_hole_visited ){
				//cout<<"repeated cell visited "<<endl;
				c->info().m_hole_state = true;
			}
			holes.erase(holes.begin() + i);
			--i;
		}		
		else{
			c->info().m_hole_visited = true;
			c->info().m_hole_state = true;
		}
	}
}

void separate_hole(Delaunay &dt, deque<Cell> &holes,vector<Cell> &modified){
	for( unsigned int i = 0; i < holes.size(); ++i ){
		Cell c = holes[i];
		c->info().m_hole_visited = false; // to return original state

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
		if( c->info().m_original_label == HOLE )		continue;// if the new label is hole

		if( produce_singularities(dt,neighbors) )
			c->info().m_label = HOLE;
		else{
			modified.push_back(c);
			c->info().m_hole_state = false;

			holes.erase( holes.begin() + i );
			--i;
			for( unsigned int j = 0; j < 4; ++j)
				update_vertex_info( dt,neighbors[j] );
		}
	}
}

deque<Cell> get_holes_solution(Delaunay &dt){
	deque<Cell> holes;
	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		if( it->info().m_solution == HOLE){
			holes.push_back(it);
		}
	}
	return holes;
}

deque<Cell> get_holes_solution_2(Delaunay &dt){
	deque<Cell> holes;
	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		if( it->info().m_label == HOLE){
			holes.push_back(it);
		}
	}
	return holes;
}

void mark_hole_visited(deque<Cell> &holes){
	for( unsigned int i = 0; i < holes.size(); ++i ){
		Cell c = holes[i];
		c->info().m_hole_visited = false;
	}
}

void solve_holes(Delaunay &dt,deque<Vertex> &singularVertices,deque<Cell> &holes,vector<Cell> &modified){
	bool draw = 0;
	int count = 1;
	while(singularVertices.size() != 0 ){
		Vertex v = singularVertices.front();
		singularVertices.pop_front();
		//if( v->point() == Point(358,	411,281) ) draw = 1;
		v->info().m_state=1;

		vector<pair<int,int>> &sing_regions = v->info().m_singular_regions_and_subgraphs;
		while( sing_regions.size() != 0){
			int label= sing_regions.front().first;
			vector<vector<Cell>> cellSets;
			set_cells_around_singular_vertex(dt,v,label,cellSets);
			/*if(draw){				
				drawMesh(dt);
				aroundVertex(dt,v,v,label);
				getchar();
			}*/
			if( cellSets.size() == 2 ){
				//cout<<"stop 1"<<endl;	
				solve_1(dt,cellSets,singularVertices,v,label,draw,holes,modified);
			}
			else if( cellSets.size() == 1 ){
				//cout<<"stop 2"<<endl;
				solve_2_2(dt,cellSets,singularVertices,v,label,draw,holes,modified);
			}
			else{
				//cout<<"stop 3"<<endl;
				solve_3(dt,cellSets,singularVertices,v,label,draw,holes,modified);
			}
			/*if(draw){
				aroundVertex(dt,v,v,label);
				getchar();
			}*/
		}
		//v->info().m_state=1;
		/*deque<Cell> ch = get_holes_solution_2(dt);
		cout<<"iteration, original and real: "<<holes.size()<<" "<<ch.size()<<endl;*/
	}

	//deque<Cell> ch = get_holes_solution_2(dt);
	//cout<<"original and real holes: "<<holes.size()<<" "<<ch.size()<<endl;

	//deque<Cell> aux_holes = holes;
	//choose_holes(dt,holes);// avoid redundance

	holes = get_holes_solution_2(dt);
	//mark_hole_visited(aux_holes);

	//cout<<"holes after chose it  "<<holes.size()<<endl;
	separate_hole(dt,holes,modified);// maybe the holes not generate singularities
	//cout<<"holes after separate it  "<<holes.size()<<endl;
}

void neighbor_state(Delaunay &dt,deque<Cell> &holes,vector<Cell> &modified){
	Cell c;
	int new_label;
	while(true){
		unsigned int index = rand()%holes.size();
		c = holes[index];
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
	c->info().m_label = new_label;
	modified.push_back(c);

	deque<Vertex> singularities;
	for( unsigned int i = 0; i < 4; ++i){
		Vertex v = c->vertex(i);
		update_vertex_info(dt,v);
		if( v->info().is_singular_2() ){
			singularities.push_back(v);
			v->info().m_state = 2;
		}
	}
	solve_holes(dt,singularities,holes,modified);	
}

void neighbor_state_2(Delaunay &dt,deque<Cell> &holes,vector<Cell> &modified){
	/*deque<Vertex> singularities;
	for(unsigned int i = 0; i < holes.size(); ++i){
		Cell c = holes[i];
		c->info().m_label = c->info().m_original_label;
		modified.push_back(c);
		for( unsigned int j = 0; j < 4; ++j){
			Vertex v = c->vertex(j);
			update_vertex_info(dt,v);
			if( v->info().is_singular_2() && v->info().m_state != 2 ){
				singularities.push_back(v);
				v->info().m_state = 2;
			}
		}
	}
	solve_holes(dt,singularities,holes,modified);*/

	deque<Vertex> singularities;
	for(unsigned int i = 0; i < holes.size(); ++i){
		Cell c = holes[i];
		//cout<<" m_label before: "<<c->info().m_label<<endl;
		c->info().m_label = c->info().m_original_label;
		double r = ((double) rand() / (RAND_MAX));
		if( r < 0.3) {//see
			int lbl = label_for_cell(dt,c);
			if(lbl != -1) {
				c->info().m_label =lbl;
			}
		}
		//cout<<" m_label before: "<<c->info().m_label<<endl;
		modified.push_back(c);
		for( unsigned int j = 0; j < 4; ++j){
			Vertex v = c->vertex(j);
			update_vertex_info(dt,v);
			if( v->info().is_singular_2() && v->info().m_state != 2 ){
				singularities.push_back(v);
				v->info().m_state = 2;
			}
		}
	}
	solve_holes(dt,singularities,holes,modified);
}

void initiate_simulated(Delaunay &dt){
	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		Cell c = it;
		// here need to change
		c->info().m_original_label = c->info().m_label;

		c->info().m_solution = c->info().m_label;
	}
}

void accept_state(deque<Cell> &holes, vector<Cell> &modified){	
	for( unsigned int i = 0; i < modified.size(); ++i){
		modified[i]->info().m_solution = modified[i]->info().m_label;
	}
	for( unsigned int i = 0; i < holes.size(); ++i){
		holes[i]->info().m_solution = holes[i]->info().m_label;
	}
}

void not_accept_state(deque<Cell> &holes, vector<Cell> &modified){
	for( unsigned int i = 0; i < modified.size(); ++i){
		modified[i]->info().m_label = modified[i]->info().m_solution;
	}
	for( unsigned int i = 0; i < holes.size(); ++i){
		holes[i]->info().m_label = holes[i]->info().m_solution;
	}
}

int check_process(Delaunay &dt){
	vector<Cell> incorrect;
	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		Cell c = it;
		if( c->info().m_label != c->info().m_solution || c->info().m_label == HOLE ){
			//cout<<"label and solution "<<c->info().m_label<<" "<<c->info().m_solution<<" "<<c->info().m_hole_state <<endl;
			incorrect.push_back(c);
		}
	}
	return (int)incorrect.size();
}

void check_lalbel(Delaunay &dt){
	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		if(it->info().m_label == -1){
			cout<<"error lalbel -1"<<endl;
			getchar();
		}
	}
}

void simulated_annealing_whithout_points(Delaunay &dt){
	double count = 0;
	std::stringstream manifold_simulated;
	cout<<"begin"<<endl;
	vector<Vertex> boundaryVertices;
	deque<Vertex> singularVertices;
	deque<Cell> holes;
	vector<Cell> modifiedCells;

	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);

	cout<<"sa number of non manifold vertices before "<<singularVertices.size()<<endl;
	cout<<"sa number of vertices "<<dt.number_of_vertices()<<endl;

	clock_t begin = clock();

	initiate_simulated(dt);// equal soluton to m_label
	solve_holes(dt,singularVertices,holes,modifiedCells);
	accept_state(holes,modifiedCells);
	
	//cout<<"stop 1 "<<endl;
	double old_cost = cost(dt,holes);
	double T = 1.0;
	double T_min = 0.00001;
	double alpha = 0.9;
	while ( T > T_min && holes.size() != 0){
		int i = 1;
		while (i <= 100 && holes.size() != 0){		
			++count;
			manifold_simulated.str("");
			manifold_simulated<<count<<"manifold_simulated";

			modifiedCells.clear();
			deque<Cell> new_holes = holes;
			neighbor_state_2(dt,new_holes,modifiedCells);// it return less holes 
			//if(new_holes.size() == 0) { T_min = 2; break;}
			double new_cost = cost(dt,new_holes);
			double ap = acceptance_probability(old_cost, new_cost, T);

			/*deque<Cell> verify = get_holes(dt);
			cout<<"verified holes before "<<verify.size()<<endl;*/
			//cout<<"old and new size costs "<<holes.size()<<"  "<<new_holes.size()<<endl;
			cout<<"old and new size holes "<<holes.size()<<"  "<<new_holes.size()<<endl;
			
			/*drawMesh(dt,1);
			getchar();*/

			if( ap > random_sa() ){
				accept_state(new_holes,modifiedCells);
				holes = new_holes;
				old_cost = new_cost;

				//deque<Cell> ch = get_holes_solution(dt);
				////cout<<"accepted"<<endl;
				//if( ch.size() != new_holes.size() ){
				//	cout<<"real and new_holes "<<ch.size()<<" * "<<new_holes.size()<<endl;
				//	cout<<"problem with holes yes"<<endl;
				//	drawMesh(dt,1);
				//	getchar();
				//}
				
				//saveMesh(dt,"chest/simulating", ( manifold_simulated.str() ).c_str()  );				
			}
			else{
				not_accept_state(new_holes,modifiedCells);		
				//deque<Cell> ch = get_holes_solution_2(dt);
				////cout<<"accepted"<<endl;
				//if( ch.size() != holes.size() ){
				//	cout<<"real and new_holes "<<ch.size()<<" * "<<holes.size()<<endl;
				//	cout<<"problem with holes not"<<endl;
				//	drawMesh(dt,1);
				//	getchar();
				//}
			}
			/*drawMesh(dt,1);
			getchar();*/

			//int check = check_process(dt);


			check_lalbel(dt);

			/*verify = get_holes(dt);
			cout<<"verified holes "<<verify.size()<<endl;*/
			i += 1;			
		}
		T = T*alpha;
	}

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<<"time simulated: "<<elapsed_secs<<endl;

	restart_vertex_state(dt);
	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);
	cout<<"number of non manifold vertices simulated "<<singularVertices.size()<<endl;
	cout<<"number of vertices simulated "<<dt.number_of_vertices()<<endl;

	vector<Point> vtkPoints;
	for(unsigned int i = 0 ; i < singularVertices.size();++i){
		vtkPoints.push_back(singularVertices[i]->point());
		cout<<singularVertices[i]->point()<<endl;
	}

	/*drawMesh(dt,vtkPoints);
	getchar();*/

	restart_vertex_state(dt);
	/*drawMesh(dt,vtkPoints);
	getchar();*/
}