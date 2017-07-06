double larger_edge(Delaunay &dt,vector<Cell> &cells){
	double max=0;
	double squaredLength;
	for(unsigned int i = 0; i < cells.size(); ++i){
		squaredLength = dt.segment(cells[i],0,1).squared_length();
		if( squaredLength > max)
			max = squaredLength;

		squaredLength = dt.segment(cells[i],0,2).squared_length();
		if( squaredLength > max)
			max = squaredLength;

		squaredLength = dt.segment(cells[i],0,3).squared_length();
		if( squaredLength > max)
			max = squaredLength;

		squaredLength = dt.segment(cells[i],1,2).squared_length();
		if( squaredLength > max)
			max = squaredLength;

		squaredLength = dt.segment(cells[i],1,3).squared_length();
		if( squaredLength > max)
			max = squaredLength;

		squaredLength = dt.segment(cells[i],2,3).squared_length();
		if( squaredLength > max)
			max = squaredLength;
	}
	return max;
}

double volume_cell_set(vector<Cell> &cells){
	double vol=0;
	for(unsigned int i = 0; i < cells.size(); ++i)
		vol+= CGAL::volume( cells[i]->vertex(0)->point(),cells[i]->vertex(1)->point(),cells[i]->vertex(2)->point(),cells[i]->vertex(3)->point() ) ;
	return vol;	
}

double criteria_1( Delaunay &dt,vector<Cell> &cells, Vertex v){
	return 1/larger_edge(dt,cells);
}

void sort_by_criteria( Delaunay &dt,vector<vector<Cell>> &cellSets, Vertex v,int &min, int &max, bool &changeable){
	vector<double> ratio(cellSets.size());
	//cout<<"sort size"<<ratio.size()<<endl;
	for(unsigned int i = 0; i < cellSets.size(); ++i)
		ratio[i] = criteria_1(dt,cellSets[i],v);
	max=0;min=1;
	if(ratio[min] > ratio[max]){
		max = 1;
		min = 0;
	}
	if( ratio[max]/ratio[min] < 2 )//for chest is different =2.5
		changeable = true;
	else
		changeable = false;
}

bool is_singular_edge(Delaunay &dt,Vertex v, Vertex w,int &opposite_label){
	vector<pair<int,int>> &regv = v->info().m_singular_regions_and_subgraphs;
	vector<pair<int,int>> &regw = w->info().m_singular_regions_and_subgraphs;

	for( unsigned int i = 0; i < regw.size(); ++i){
		if( regw[i].first == opposite_label){
			if( regw[i].second == 2 ){
				//cout<<"size w sg"<<regw[i].second<<endl;
				return true;
			}
		}
	}


	return false;
	/*if( regv[0].first == regw[1].first && regv[1].first == regw[0].first )
		return true;
	else
		return false;*/
}


bool extreme_vertex_of_singular_edge(Delaunay &dt,Vertex &v,Vertex &w,int &opposite_label){
	list<Vertex> vertices;
	dt.finite_incident_vertices(v,back_inserter(vertices));
	for(std::list<Vertex>::iterator it = vertices.begin(); it != vertices.end(); it++){ 
		Vertex vw =*it;
		if( vw->info().is_singular_2() ){
			if( is_singular_edge(dt,v,vw,opposite_label) ){
				w = vw;
				return true;
			}
		}
	}
	return false;
}

