bool testCell(Delaunay & dt,Cell &c, int label,Vertex &v){
	for(unsigned int i = 0; i < 4; ++i){
		update_regions_around_vertex(dt,c->vertex(i) );
		if( c->vertex(i)->info().is_singular() ){
			v = c->vertex(i);
			cout<<"singular region: "<<c->vertex(i)->info().m_singular_regions_around[0]<<endl;
			return false;
		}
	}
	return true;
}

void checking(Delaunay &dt,Vertex v){
	list<Vertex> incident;
	dt.finite_incident_vertices(v,back_inserter(incident));	
	for(std::list<Vertex>::iterator it = incident.begin(); it != incident.end(); it++){
		Vertex vv = *it;
		update_regions_around_vertex(dt,vv);
		if(vv->info().is_singular()){
			cout<<"epic fail"<<endl;
			cout<<v->point()<<endl;
			cout<<vv->point()<<endl;
			getchar();
		}
	}
	update_regions_around_vertex(dt,v);
	if(v->info().is_singular()){
		cout<<"epic fail"<<endl;
		cout<<v->point()<<endl;
		getchar();
	}
}

//bool is_label_around(Delaunay &dt, Vertex &v, int label){
//	list<Cell> incident;
//	dt.finite_incident_cells(v,back_inserter(incident));
//	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
//		if( (*it)->info().m_label == label )
//			return true;
//	}
//	return false;
//}
//
//bool add_tetrahedron(Delaunay &dt, Cell &c, int label){
//	if( c->info().m_label == label) return false;
//	Cell neighborCell;
//	vector<int> facets;
//	for(unsigned int i = 0; i < 4; ++i){
//		neighborCell = c->neighbor(i);
//		if(neighborCell->info().m_label == label)
//			facets.push_back(i);
//	}
//
//	if( facets.size() == 0 ){
//		for(unsigned int i = 0; i < 4; ++i){
//			if( is_label_around(dt,c->vertex(i),label ) )
//				return false;
//		}
//		return true;
//	}		
//	else if( facets.size() == 1 )
//		return !is_label_around(dt,c->vertex(facets[0]),label);
//	
//	else if( facets.size() == 2 ){// i am checking the self tetra
//		Cell_circulator cc_start = dt.incident_cells( c,facets[0],facets[1] );
//		Cell_circulator cc = cc_start;	
//		//cout<<"f = 2"<<endl;
//		do{
//			if(cc->info().m_label == label)
//				return false;
//			++cc;
//		} while (cc!= cc_start);
//		return true;
//	}	
//	else	return true;
//}
//
//bool is_singular_after_relabeling(Delaunay &dt, Vertex &v){
//	vector<int> labels;
//	list<Cell> incident;
//	bool cen = false;
//	dt.finite_incident_cells(v,back_inserter(incident));	
//	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
//		if( (*it)->info().m_visited ) continue;
//		if( (*it)->info().m_label == -1 ) continue;
//		int cellLabel =  (*it)->info().m_label;
//		if( std::find(labels.begin(), labels.end(), cellLabel ) != labels.end() ){//found
//			cen = true;
//			break;
//		}
//		else{
//			labels.push_back( cellLabel );
//			visite_cells(dt, (*it) ,v , cellLabel );
//		}
//	}
//	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++)
//		(*it)->info().m_visited = false;
//	return cen;
//}
//
//void neighbor_vertices(Delaunay &dt, vector<Cell> &cells,vector<Vertex> &vertices,Vertex &v) {
//	Vertex v1;
//	for(unsigned int i = 0; i < cells.size(); ++i){
//		int idx = cells[i]->index(v);
//		for(unsigned int j = 0; j < 3; ++j){
//			v1 = cells[i]->vertex( dt.vertex_triple_index(idx,j) );
//			if( v1->info().m_auxiliar == false ){
//				v1->info().m_auxiliar = true;
//				vertices.push_back(v1);
//			}
//		}
//	}
//	for(unsigned int i = 0; i < vertices.size(); ++i)
//		vertices[i]->info().m_auxiliar = false;
//}
//
//bool check_neighbor_vertices( Delaunay &dt, vector<Vertex> &vertices ){
//	for(unsigned int i = 0; i < vertices.size(); ++i ){
//		if( vertices[i]->info().has_more_than_two_regions() || vertices[i]->info().m_boundary){
//			if( is_singular_after_relabeling(dt, vertices[i] ) )
//				return false;
//		}
//	}
//	return true;
//}
//
//bool relabel_if_possible(Delaunay & dt, vector<Cell> &cells, int label,Vertex v){
//	deque<Cell> cpyCells(cells.begin(),cells.end());
//	int cen = cpyCells.size();
//	Cell c;
//	int old = cells[0]->info().m_label;
//
//	/*cout<<endl<<"size cpycells: "<<cpyCells.size()<<endl;
//	if(v->point() == Point(5,347,94)){
//		vector<Point> points;
//		points.push_back(v->point() );
//		points.push_back(Point(3,432,5) );
//		drawMesh(dt);
//		vertexRegion(dt,v,points,v->info().m_regions_around[ v->info().m_regions_around[0] == label  ]);
//		getchar();
//	}*/
//
//	while(cpyCells.size() != 0 && cen !=0){		
//		c = cpyCells.back();
//		cpyCells.pop_back();
//		if( add_tetrahedron(dt,c,label) ){
//			c->info().m_label = label;
//			cen = cpyCells.size();	
//
//			/*Vertex ver;
//			if(testCell(dt,c,label,ver)) cout<<"can add"<<endl;
//			else cout<<"can not add"<<endl;
//
//			if(v->point() == Point(5,347,94)){
//				if(!testCell(dt,c,label,ver)){
//					vector<Point> points;
//					points.push_back(ver->point() );
//					points.push_back( v->point() );
//					drawMesh(dt);
//					vertexRegion(dt,ver,points,ver->info().m_regions_around[ ver->info().m_regions_around[0] == label  ]);
//					getchar();
//				}
//			}*/
//		}
//		else{
//			cpyCells.push_front(c);
//			--cen;
//		}
//	}
//
//	vector<Vertex> neighbors;
//	neighbor_vertices(dt,cells,neighbors,v);
//
//	if( cpyCells.size() != 0  ){
//		for(unsigned int i = 0; i < cells.size(); ++i)
//			cells[i]->info().m_label = old;
//		return false;
//	}
//	else{		
//		if( !check_neighbor_vertices(dt,neighbors) ){
//			for(unsigned int i = 0; i < cells.size(); ++i)
//				cells[i]->info().m_label = old;
//			return false;
//		}
//		else{
//			for(unsigned int i = 0; i < neighbors.size(); ++i ){
//				update_regions_around_vertex( dt,neighbors[i] );
//			}
//			update_regions_around_vertex( dt, v );
//			return true;
//		}		
//	}	
//}

bool remove_pick(Delaunay & dt, Vertex &v, double minAngle){
	int label0 = v->info().m_regions_around[0];
	int label1 = v->info().m_regions_around[1];
	vector<Cell> cells0;
	vector<Cell> cells1;
	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));
	Vertex v1,v2,v3;
	double total = 0;
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		if( (*it)->info().m_label == label0 ){
			cells0.push_back(*it);
			Cell cell = *it;
			int indx = cell->index(v);	
			v1 = cell->vertex( dt.vertex_triple_index(indx,0) );
			v2 = cell->vertex( dt.vertex_triple_index(indx,1) );
			v3 = cell->vertex( dt.vertex_triple_index(indx,2) );
			total += solid_angle( cell->vertex(indx)->point(),v1->point(),v2->point(),v3->point() );
		}else{
			cells1.push_back(*it);
		}
	}

	if( total > 360){
		total = 720 - total;
		if( total < minAngle){
			return relabel_if_possible(dt,cells1,label0,v);
		}
	}
	else{
		if( total < minAngle){
			return relabel_if_possible(dt,cells0,label1,v) ;
		}
	}
	return false;
}

void remove_pick2(Delaunay & dt, Vertex &v, double minAngle){

	map<int,int> label_to_index;
	for(unsigned int i = 0; i < v->info().m_regions_around.size(); ++i){
		label_to_index[ v->info().m_regions_around[i] ] = i;
	}

	vector<vector<Cell>> cellsSet( v->info().m_regions_around.size() );	
	vector<double> regionAngle( v->info().m_regions_around.size() ,0 );	

	Vertex v1,v2,v3;
	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));	
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		Cell c =*it;
		cellsSet[ label_to_index[ c->info().m_label ]  ].push_back(c);
		int indx = c->index(v);	
		v1 = c->vertex( dt.vertex_triple_index(indx,0) );
		v2 = c->vertex( dt.vertex_triple_index(indx,1) );
		v3 = c->vertex( dt.vertex_triple_index(indx,2) );
		regionAngle[ label_to_index[ c->info().m_label ] ] += solid_angle( c->vertex(indx)->point(),v1->point(),v2->point(),v3->point() );
	}

	for(unsigned int i = 0; i < regionAngle.size(); ++i){
		//cout<<"angle: "<<i<<" "<<regionAngle[i]<<endl;
		if( regionAngle[i] < minAngle ){
			int newlabel = new_label_for_relabeling(dt,cellsSet[i],v,v->info().m_regions_around[i]);
			relabel_if_possible(dt,cellsSet[i],v->info().m_regions_around[i],v);
		}
	}
}

void remove_boundary_picks(Delaunay &dt, double minAngle){
	vector<Vertex> boundaryVertices;
	boundary_vertices(dt,boundaryVertices);

	//vector<Point> points;
	//for(unsigned int i = 0 ; i< boundaryVertices.size(); ++i){
	//	if( boundaryVertices[i]->info().has_more_than_two_regions() ){
	//		points.push_back(boundaryVertices[i]->point());
	//		//cout<<"labels: "<<boundaryVertices[i]->info().m_regions_around.size()<<endl;
	//	}
	//}

	//drawMesh(dt,points);
	//getchar();

	for(unsigned int i = 0 ; i< boundaryVertices.size(); ++i){
		Vertex v = boundaryVertices[i];
		if( v->info().is_singular() || v->info().is_corner() ) continue;
		if( v->info().has_two_regions() )
			remove_pick( dt,v,minAngle );
		else{
			if( v->info().has_more_than_two_regions() ){
				//cout<<"has_more_than_two_regions"<<endl;
				remove_pick2( dt,v,minAngle );
			}
		}
	}
	restart_vertex_state(dt);
	cout<<"ends remove boundary picks"<<endl;
}