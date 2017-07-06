
bool erode( Delaunay &dt,vector<vector<Cell>> &cellSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw, deque<Cell> &holes ){	
	if(cellSets.size() != 2 ) return false;
	int min,max;
	bool change = false;
	//sort_by_criteria(dt,cellSets,v,min,max,change);	
	sort_by_criteria_w(dt,cellSets,v,min,max,change);	

	//if( !insert_to_erode(dt,singularities,cellSets,v,label,min,max,holes,draw) ){
	//	int newlabel = new_label_for_relabeling(dt,cellSets[max],v,label); 
	//	if( !relabel_if_possible_3(dt,cellSets[max],newlabel,v,singularities ) ){
	//		if( change){
	//			//cout<<"eroded changed"<<endl;
	//			newlabel = new_label_for_relabeling(dt,cellSets[min],v,label); 
	//			if( !relabel_if_possible_3(dt,cellSets[min],newlabel,v,singularities) ){
	//				return false;
	//			}
	//		}
	//		else return false;
	//	}
	//}

	int newlabel = new_label_for_relabeling(dt,cellSets[max],v,label); 
	if( !relabel_if_possible_3(dt,cellSets[max],newlabel,v,singularities) ){	
		if( change){
			newlabel = new_label_for_relabeling(dt,cellSets[min],v,label); 
			//cout<<"new label 2 "<<newlabel<<endl;
			if( !relabel_if_possible_3(dt,cellSets[min],newlabel,v,singularities ) ){						
				if( !insert_to_erode(dt,singularities,cellSets,v,label,min,max,holes,draw) ){
					return false;
				}
			}
		}
		else{
			if( !insert_to_erode(dt,singularities,cellSets,v,label,min,max,holes,draw) ){
				return false;
			}
		}
	}	

	// always uncommented
	return true;
}

vector<pair<vector<Cell> , double>> order_set_by_criteria(Delaunay &dt, Vertex v, int label){
	vector<vector<Cell>> cells_by_label;
	set_cells_by_region(dt,v,cells_by_label,label);

	vector<pair<vector<Cell> , double>> cells_criteria;
	for(unsigned int i = 0; i < cells_by_label.size(); ++i){
		double criteria = criteria_1(dt,cells_by_label[i],v);
		cells_criteria.push_back(make_pair(cells_by_label[i],criteria ) );
	}

	std::sort(cells_criteria.begin(), cells_criteria.end(), 
		boost::bind(&std::pair<vector<Cell> , double>::second, _1) > boost::bind(&std::pair<vector<Cell> , double>::second, _2));
	return cells_criteria;
}

void repairing_walk_graph(Delaunay &dt,Cell &c,vector<Cell> &new_group, Vertex &v,int &label){	
	vector<Cell> visitedCells;
	visitedCells.push_back(c);
	new_group.push_back(c);
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
				new_group.push_back(c0);
				c0->info().m_visited = true;
			}
		}
	}
}

void separate_cells_for_repairing_1(Delaunay &dt,vector<vector<Cell>> &groups, Vertex v, int label){
	int nro_subgraphs = 0;
	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));	

	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		if( (*it)->info().m_visited ) continue;
		if( (*it)->info().m_label == -1 || (*it)->info().m_label == label) continue;
		int cellLabel =  (*it)->info().m_label;
		if( cellLabel != label){
			vector<Cell> new_group;
			repairing_walk_graph(dt,(*it),new_group,v,label);
			//cout<<" new_group size "<<new_group.size()<<endl;
			groups.push_back(new_group);
		}
	}
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++)
		(*it)->info().m_visited = false;
}

bool repairing_1(Delaunay &dt, vector< vector<Cell> > &cellSets, deque<Vertex> &singularities, Vertex v, int label, bool draw, deque<Cell> &holes ) {
	vector<vector<Cell>> groups;
	separate_cells_for_repairing_1(dt,groups,v,label);
	vector<pair<vector<Cell> , double>> cells_criteria;
	for(unsigned int i = 0; i < groups.size(); ++i){
		double criteria = criteria_1(dt,groups[i],v);
		cells_criteria.push_back(make_pair(groups[i],criteria ) );
		//cout<<"criteria and size: "<<criteria<<" "<<groups[i].size()<<endl;
	}
	std::sort(cells_criteria.begin(), cells_criteria.end(), 
		boost::bind(&std::pair<vector<Cell> , double>::second, _1) > boost::bind(&std::pair<vector<Cell> , double>::second, _2));

	if( cells_criteria.size() == 2){
		vector<vector<Cell>> separate_set;
		separate_set.push_back(cells_criteria[1].first);
		separate_set.push_back(cells_criteria[0].first);
		//if( insert_to_erode(dt,singularities,separate_set,v,label,0,1,holes,draw,1) ){
		if(1==2){
			return true;
		}
		else{// find the path 
			//return insert_to_dilate(dt,singularities,separate_set,v,separate_set[0][0]->info().m_label,holes,draw);
			//return dilate_with_checking(dt,separate_set,singularities,v,label,draw) ;
			int newlabel = new_label_for_relabeling(dt,cells_criteria[0].first,v,cells_criteria[0].first[0]->info().m_label);
			if( !relabel_if_possible_3(dt,cells_criteria[0].first,newlabel,v,singularities ) ){//original ==> 0
				if( insert_to_erode(dt,singularities,separate_set,v,label,0,1,holes,draw,1) )	return true;
				else	return dilate_with_checking(dt,separate_set,singularities,v,label,draw) ;

				//return dilate_with_checking(dt,separate_set,singularities,v,label,draw) ;
			}
			else	return true;
		}
	}
	else if( cells_criteria.size() == 1 ){
		if( relabel_if_possible_3(dt,cells_criteria[0].first,label,v,singularities ) ){
			return true;
		}
	}
	else{
		//if cells_criteria = 1
		//cout<<"cells size +2 "<<cells_criteria.size()<<endl;
		/*for(unsigned int i = 0; i < cells_criteria.size(); ++i){
			cout<<i<<" label: "<<cells_criteria[i].first[0]->info().m_label<<endl;
		}*/
		vector< vector<Cell> > sets;
		sets.push_back( cells_criteria[1].first );
		sets.push_back( cells_criteria.back().first );
		if( dilate_with_checking(dt,sets,singularities,v,label,draw) ){
			cells_criteria = order_set_by_criteria(dt,v,label);
			if( !v->info().is_singular_2() ){
				//cout<<"cells size +2, dilated"<<endl;
				return true;
			}
			else cells_criteria = order_set_by_criteria(dt,v,label);
		}

		for(unsigned int i = 0; i < cells_criteria.size(); ++i){// for all the multimaterial meshes
			int newlabel = new_label_for_relabeling(dt,cells_criteria[i].first,v,cells_criteria[i].first[0]->info().m_label);		
			if( relabel_if_possible_3(dt,cells_criteria[i].first,newlabel,v,singularities ) ){
				if( !v->info().is_singular_2() ){
					/*if(draw){		
						drawMesh(dt);
						aroundVertex(dt,v,v,label);
						getchar();
					}*/
					return true;
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

	return false;
}

bool repairing_2(Delaunay &dt, vector< vector<Cell> > &cellSets, deque<Vertex> &singularities, Vertex v, int label, bool draw, deque<Cell> &holes )
{
	
	//if(v->info().has_more_than_two_regions() ) cout<<"+ 2 regions"<<endl;
	//if( ( !is_erodible(dt,cellSets[0],v,label) && !is_erodible(dt,cellSets[1],v,label) ) || v->info().has_more_than_two_regions() ){		

	//if( ( !is_erodible(dt,cellSets[0],v,label) ||  !is_erodible(dt,cellSets[1],v,label) ) ){		
	if(1 == 2){
		if( v->info().has_two_regions() ){
			//cout<<"dilate stop 0"<<endl;
			Vertex w;
			int opposite_label = v->info().m_regions_around[ v->info().m_regions_around[0] == label ];
			if( extreme_vertex_of_singular_edge(dt,v,w,opposite_label ) ){
				//cout<<"dilate stop 1"<<endl;
				vector<	Point> points;
				points.push_back( w->point() );
				/*if(draw){				
					drawMesh(dt,points);
					aroundVertex(dt,w,w,opposite_label);
					getchar();
				}*/
				vector<vector<Cell>> cellSets_2;
				set_cells_around_singular_vertex(dt,w,opposite_label,cellSets_2);// get the other singular regions
				//cout<<"dilate stop 2 - "<<cellSets_2.size()<<endl;
				if( erode(dt,cellSets_2,singularities,w,opposite_label,draw,holes) ){
					//if(draw){				
					//	//drawMesh(dt,points);
					//	aroundVertex(dt,w,w,opposite_label);
					//	getchar();
					//}
				}				
			}
			else{// it is probably a singular vertex
				//cout<<"dilate stop 2"<<endl;
				/*if( !dilate_with_checking(dt,cellSets,singularities,v,label,draw) )
					insert_to_dilate(dt,singularities,cellSets,v,label,holes,draw);*/
				//cout<<"dilate stop 2"<<endl;
				//if( !insert_to_dilate(dt,singularities,cellSets,v,label,holes,draw) )
					dilate_with_checking(dt,cellSets,singularities,v,label,draw);
				//else return true;//recently
			}
		}
		else{
			//cout<<"dilate stop 3"<<endl;
			//if( !insert_to_dilate(dt,singularities,cellSets,v,label,holes,draw) ) // is not the final code
				dilate_with_checking(dt,cellSets,singularities,v,label,draw);

			//else return true;//recently
		}
	}	
	set_cells_around_singular_vertex(dt,v,label,cellSets);
	
	if( cellSets.size() == 1 ){
		return true;
	}
	else{
		//cout<<"erode, then dilate"<<endl;
		if( erode(dt, cellSets,singularities,v,label,draw,holes) ){
			//cout<<" then eroded"<<endl;
			return true;
		}
		else
			 return dilate_with_checking(dt,cellSets,singularities,v,label,draw);
		//return erode_with_checking(dt, cellSets,singularities,v,label,draw,holes);	
	}
}

bool repairing_2_1(Delaunay &dt, vector< vector<Cell> > &cellSets, deque<Vertex> &singularities, Vertex v, int label, bool draw, deque<Cell> &holes ) {
	/*if(draw){
		drawMesh(dt);
		vector<Point> pts;
		two_ring(dt,v,pts,label);
		getchar();
	}*/
	//cout<<"is erodible the other side: "<<is_erodible(dt,cellSets[1],v,label)<<endl;
	if( !is_erodible(dt,cellSets[0],v,label) || !is_erodible(dt,cellSets[1],v,label) ){		
		//cout<<"dilate"<<endl;
		if( v->info().has_two_regions() ){
			Vertex w;
			int opposite_label = v->info().m_regions_around[ v->info().m_regions_around[0] == label ];
			if( extreme_vertex_of_singular_edge(dt,v,w,opposite_label ) ){
				vector<	Point> points;
				points.push_back( w->point() );
				if(draw){				
					drawMesh(dt,points);
					aroundVertex(dt,w,w,opposite_label);
					getchar();
				}
				vector<vector<Cell>> cellSets_2;
				set_cells_around_singular_vertex(dt,w,opposite_label,cellSets_2);// get the other singular regions
				if( erode(dt, cellSets_2,singularities,w,opposite_label,draw,holes) ){
					if(draw){				
						//drawMesh(dt,points);
						aroundVertex(dt,w,w,opposite_label);
						getchar();
					}
				}				
			}
			else{// it is probably a singular vertex
				if( !dilate_with_checking(dt,cellSets,singularities,v,label,draw) )
					insert_to_dilate(dt,singularities,cellSets,v,label,holes,draw);
			}
		}
		else{
			dilate_with_checking(dt,cellSets,singularities,v,label,draw);
		}
	}

	set_cells_around_singular_vertex(dt,v,label,cellSets);
	if( cellSets.size() == 1 )
		return true;
	else{
		if( erode(dt, cellSets,singularities,v,label,draw,holes) ){
			return true;
		}
		else
			 return dilate_with_checking(dt,cellSets,singularities,v,label,draw);
		//return erode_with_checking(dt, cellSets,singularities,v,label,draw,holes);	
	}
}

bool dilate_reparing_3(Delaunay &dt, vector< vector<Cell> > &cellSets, deque<Vertex> &singularities, Vertex v, int label, bool draw, deque<Cell> &holes ){
	int cen = true ;
	while( cellSets.size() > 1 ){
		cen = true;
		for( unsigned int i = 0; i < cellSets.size() - 1; ++i){
			vector< vector<Cell> > sets;
			sets.push_back( cellSets[i] );
			sets.push_back( cellSets[i + 1] );
			if( dilate_with_checking(dt,sets,singularities,v,label,draw) ){// dilated
				cen = false;
				set_cells_around_singular_vertex(dt,v,label,cellSets);
				break;
			}
		}
		if (cen) return false;
	}
	return true;
}

bool set_is_erodible(Delaunay &dt, vector< vector<Cell> > &cellSets,Vertex v,int label){
	double count = 0;
	for(unsigned int i = 0; i < cellSets.size(); ++i){
		if( !is_erodible(dt,cellSets[i],v,label) )
			++count;
	}
	if( count > 0) return false;
	return true;
}

bool repairing_3(Delaunay &dt, vector< vector<Cell> > &cellSets, deque<Vertex> &singularities, Vertex v, int label, bool draw, deque<Cell> &holes ) {

	//if(draw){		
	//	//cout<<"size: "<<cellSets.size()<<endl;
	//	drawMesh(dt);
	//	aroundVertex(dt,v,v,label);
	//	getchar();
	//}

	/*if(draw){
		cout<<"2 ring : "<<endl;
		drawMesh(dt);
		vector<Point> pts;
		two_ring(dt,v,pts,label);
		getchar();
	}*/

	bool result = false;
	if( !set_is_erodible(dt,cellSets,v,label) || v->info().has_more_than_two_regions() )
		result = dilate_reparing_3(dt,cellSets,singularities,v,label,draw,holes);

	if ( result ) return true;

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
		//if( l != v1 ){
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
		//cout<<"size: "<<cellSets.size()<<endl;
		drawMesh(dt);
		aroundVertex(dt,v,v,label);
		getchar();
		// = 0;
	}
	return repairing_2(dt,cellSets,singularities,v,label,draw,holes);
}

void repairing(Delaunay &dt){
	mark_boundary(dt);
	deque<Cell> holes;
	vector<Vertex> boundaries;
	deque<Vertex> singularities;			

	clock_t begin0 = clock();

	boundary_vertices(dt,boundaries);
	singular_vertices(dt,boundaries,singularities);

	clock_t end0 = clock();
	double elapsed_secs0 = double(end0 - begin0) / (CLOCKS_PER_SEC/1000);
	cout<<"time get singularities: "<<elapsed_secs0<<endl;

	cout<<"number of non manifold vertices before "<<singularities.size()<<endl;
	cout<<"number of vertices "<<dt.number_of_vertices()<<endl;
	
	/*vector<Point> vtkPoints;
	for(unsigned int i = 0 ; i < singularities.size();++i){
		vtkPoints.push_back( singularities[i]->point() );
	}
	off_file_points(vtkPoints,"hyena/manifold","off_points" );*/


	/*savePoints2(vtkPoints,"hyena/original/singularities7.vtk",7);
	savePoints2(vtkPoints,"hyena/original/singularities6.vtk",6);
	savePoints2(vtkPoints,"hyena/original/singularities5.vtk",5); 
	savePoints2(vtkPoints,"hyena/original/singularities3.vtk",3); 
	savePoints2(vtkPoints,"hyena/original/singularities2.vtk",2); 
	savePoints2(vtkPoints,"hyena/original/singularities15.vtk",1.5);
	savePoints2(vtkPoints,"hyena/original/singularities1.vtk",1); 
	cout<<"finish off singularities"<<endl;
	getchar();
	*/


	clock_t begin = clock();

	bool draw = 1;	
	int count = 0;
	while(singularities.size() != 0 ){	
		++count;		
		//if( count == 296 ) draw = 1;// 3570 - 3590    ------> 3587 bad solution for singularity 590 for head

		/*draw = 0;
		if( count == 250 ) draw = 1;
		if( count == 255 ) draw = 1;
		if( count == 260 ) draw = 1;
		if( count == 265 ) draw = 1;
		if( count == 270 ) draw = 1;
		if( count == 275 ) draw = 1;
		if( count == 280 ) draw = 1;
		if( count == 285 ) draw = 1;
		if( count == 290 ) draw = 1;
		if( count == 295 ) draw = 1;
		if( count == 300 ) draw = 1;
		if( count == 400 ) draw = 1;
		if( count == 450 ) draw = 1;
		if( count == 500 ) draw = 1;*/

		//if( count == 3200 ) getchar();
		//cout<<"count: "<<count<<endl;

		//cout<<"count: "<<count<<endl;
		Vertex v = singularities.front();
		//Vertex v = singularities[1];

		singularities.pop_front();
		v->info().m_state = 1;	
		
		vector<pair<int,int>> &sing_regions = v->info().m_singular_regions_and_subgraphs;
		
		while( sing_regions.size() != 0){
			int label = sing_regions.front().first;
			//cout<<"label: "<<label<<endl;
			//int label = sing_regions[1].first;
			vector<vector<Cell>> cellSets;
			set_cells_around_singular_vertex(dt,v,label,cellSets);
			if(draw){				
				//save_vertex_star(dt,"identification/relabeling","3comp_original");
				drawMesh(dt);
				aroundVertex(dt,v,v,label);
				getchar();
			}
			
			if( cellSets.size() == 1 ){
				//cout<<"case 1"<<endl;
				if( !repairing_1(dt,cellSets,singularities,v,label,draw,holes) ){
					//cout<<"stop 1"<<endl;
					v->info().m_state=2;
					break;
				}
				//break;
			}
			else if( cellSets.size() == 2 ){
				//cout<<"case 2"<<endl;
				if( !repairing_2(dt,cellSets,singularities,v,label,draw,holes) ){
					//cout<<"stop 2"<<endl;
					v->info().m_state=2;
					break;
				}
			}
			else{
				//cout<<"case 3"<<endl;
				if( !repairing_3(dt,cellSets,singularities,v,label,draw,holes) ){
					//cout<<"stop 3"<<endl;
					v->info().m_state=2;
					break;
				}
			}
			//cout<<"size sing_regions: "<<sing_regions.size()<<endl;
			if(draw){				
				aroundVertex(dt,v,v,label);
				getchar();
			}
		}
	}

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<<"time repairing: "<<elapsed_secs<<endl;

	restart_vertex_state(dt);
	boundary_vertices(dt,boundaries);
	singular_vertices(dt,boundaries,singularities);
	cout<<"number of non manifold vertices after "<<singularities.size()<<endl;
	cout<<"number of vertices "<<dt.number_of_vertices()<<endl;

	/*vtkPoints.clear();
	for(unsigned int i = 0 ; i < singularities.size();++i){
		vtkPoints.push_back( singularities[i]->point() );
	}*/

	//vector<Point> vtkPoints;

	/*vtkPoints.clear();
	for(unsigned int i = 0 ; i < singularities.size();++i){
		vtkPoints.push_back( singularities[i]->point() );
	}

	savePoints2(vtkPoints,"head/manifold/remaining7.vtk",7);
	savePoints2(vtkPoints,"head/manifold/remaining5.vtk",5); 
	savePoints2(vtkPoints,"head/manifold/remaining3.vtk",3); 
	savePoints2(vtkPoints,"head/manifold/remaining1.vtk",1); */

	restart_vertex_state(dt);
}