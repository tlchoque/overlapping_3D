
void visite_cells(Delaunay &dt,Cell &c, Vertex &v,int &label){
	vector<Cell> visitedCells;
	visitedCells.push_back(c);
	Cell vc,c0;
	while(visitedCells.size() != 0){
		vc = visitedCells.back();
		visitedCells.pop_back();
		int index = vc->index(v);
		for(unsigned int i = 0; i < 3;++i){
			c0 =  vc->neighbor( dt.vertex_triple_index(index,i) );
			if( !c0->info().m_visited && c0->info().m_label == label ){
				visitedCells.push_back(c0);
				c0->info().m_visited = true;
			}
		}
	}
}

void update_regions_around_vertex(Delaunay &dt, Vertex &v){
	vector<int> &regions = v->info().m_singular_regions_around;
	vector<int> &labels = v->info().m_regions_around;
	regions.clear();
	labels.clear();

	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));	
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		if( (*it)->info().m_visited ) continue;
		if( (*it)->info().m_label == -1 ) continue;
		int cellLabel =  (*it)->info().m_label;
		if( std::find(labels.begin(), labels.end(), cellLabel ) != labels.end() ){//found
			if( !(std::find(regions.begin(), regions.end(), cellLabel ) != regions.end()) )
				regions.push_back( cellLabel );
		}
		else{
			labels.push_back( cellLabel );
			visite_cells(dt, (*it) ,v , cellLabel );
		}
	}
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++)
		(*it)->info().m_visited = false;
}

void walk_graph(Delaunay &dt,Cell &c, Vertex &v,int &label){
	vector<Cell> visitedCells;
	visitedCells.push_back(c);
	c->info().m_visited = true;
	Cell vc,c0;
	while(visitedCells.size() != 0){
		vc = visitedCells.back();
		visitedCells.pop_back();
		int index = vc->index(v);
		for(unsigned int i = 0; i < 3;++i){
			c0 =  vc->neighbor( dt.vertex_triple_index(index,i) );
			if( !c0->info().m_visited && c0->info().m_label == label ){
				visitedCells.push_back(c0);
				c0->info().m_visited = true;
			}
		}
	}
}

void walk_graph_2(Delaunay &dt,Cell &c, Vertex &v,int &label){
	vector<Cell> visitedCells;
	visitedCells.push_back(c);
	c->info().m_visited = true;

	Cell vc,c0;
	while(visitedCells.size() != 0){
		vc = visitedCells.back();
		visitedCells.pop_back();
		int index = vc->index(v);
		for(unsigned int i = 0; i < 3;++i){
			c0 =  vc->neighbor( dt.vertex_triple_index(index,i) );
			if( !c0->info().m_visited && c0->info().m_label != label ){
				visitedCells.push_back(c0);
				c0->info().m_visited = true;
			}
		}
	}
}

bool test_graph(Delaunay &dt, Vertex &v, int label,int &nro_subgraphs){
	nro_subgraphs = 0;
	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));	
	int count=0;//number of subgraphs
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		if( (*it)->info().m_visited ) continue;
		if( (*it)->info().m_label == -1 ) continue;
		int cellLabel =  (*it)->info().m_label;
		if( cellLabel == label){
			walk_graph(dt,(*it),v,cellLabel);
			++nro_subgraphs;
			++count;
		}
		else{
			walk_graph_2(dt,(*it),v,label);
			++count;
		}
	}
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++)
		(*it)->info().m_visited = false;

	if( count > 2) return true;
	else return false;
}

void update_vertex_info(Delaunay &dt, Vertex &v){
	vector<pair<int,int>> &regions = v->info().m_singular_regions_and_subgraphs;
	vector<int> &labels = v->info().m_regions_around;
	regions.clear();
	labels.clear();

	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));	
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		if( (*it)->info().m_label == -1 ) continue;
		int cellLabel =  (*it)->info().m_label;

		if( std::find(labels.begin(), labels.end(), cellLabel ) == labels.end() ){//not found
			labels.push_back( cellLabel );
			if(cellLabel == HOLE) continue;			
			int nro_subgraphs;			
			if(test_graph( dt,v,cellLabel,nro_subgraphs ) )
				regions.push_back( make_pair( cellLabel, nro_subgraphs ) );
		}
	}
	std::sort(regions.begin(), regions.end(), boost::bind(&std::pair<int, int>::second, _1) > boost::bind(&std::pair<int, int>::second, _2));

	/*for(unsigned int i = 0; i < regions.size(); ++i){
		cout<<"label and nro sub: "<<v->point()<<" "<<regions[i].first<<" "<<regions[i].second<<endl;
	}*/	
}

bool is_singular_after_relabeling_2(Delaunay &dt, Vertex &v){
	vector<int> labels ;
	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));	
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		if( (*it)->info().m_label == -1 ) continue;
		int cellLabel =  (*it)->info().m_label;

		if( std::find(labels.begin(), labels.end(), cellLabel ) == labels.end() ){//not found
			labels.push_back( cellLabel );
			if(cellLabel == HOLE) continue;			
			int nro_subgraphs;			
			if(test_graph( dt,v,cellLabel,nro_subgraphs ) )
				return true;
		}
	}
	return false;
}

void update_vertex_info_sa(Delaunay &dt, Vertex &v){
	vector<pair<int,int>> &regions = v->info().m_singular_regions_and_subgraphs;
	vector<int> &labels = v->info().m_regions_around;
	regions.clear();
	labels.clear();

	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));	
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		if( (*it)->info().m_label == -1 ) continue;
		int cellLabel =  (*it)->info().m_label;

		if( std::find(labels.begin(), labels.end(), cellLabel ) == labels.end() ){//not found
			labels.push_back( cellLabel );
			if(cellLabel == HOLE) continue;			
			int nro_subgraphs;			
			if(test_graph( dt,v,cellLabel,nro_subgraphs ) )
				regions.push_back( make_pair( cellLabel, nro_subgraphs ) );
		}
	}
	std::sort(regions.begin(), regions.end(), boost::bind(&std::pair<int, int>::second, _1) > boost::bind(&std::pair<int, int>::second, _2));
}