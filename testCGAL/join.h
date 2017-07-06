
bool is_boundary_cell(Delaunay &dt,Cell c, Vertex v){
	int label = c->info().m_label;
	int vindex = c->index(v);
	for( unsigned int i = 0; i < 3; ++i){
		int index = dt.vertex_triple_index(vindex,i);
		Cell nc = c->neighbor(index);
		if( nc->info().m_label != label)
			return true;
	} 
	return false;
}

int new_destiny(Delaunay &dt, vector<Point> &group, Point cen){
	double distance = 10000000000;
	Vector a = cen - CGAL::ORIGIN;
	int index = 0;
	for(unsigned int i = 0; i < group.size(); ++i){
		Vector b = group[i] - CGAL::ORIGIN;
		double current = a*b;
		if( current < distance){
			distance = current;
			index = i;
		}
	}
	return index;
}

bool path_to_join_1(Delaunay &dt, vector<vector<Cell>> &cellSets,Vertex &v,vector<Cell> &path){
	path.clear();
	if( cellSets.size() != 2 ) return false;	
	vector<Point> agroup,bgroup;
	vector<Cell> acells,bcells;
	for( unsigned int i = 0; i < cellSets[0].size(); ++i){// i can take jus the boundary cells
		Cell c = cellSets[0][i];
		c->info().m_path_state = true;
		if( !is_boundary_cell(dt,c,v) ) continue;		
		acells.push_back(c);
		Point cen = CGAL::centroid(c->vertex(0)->point(), c->vertex(1)->point(), c->vertex(2)->point(), c->vertex(3)->point() );
		agroup.push_back(cen);
	}
	for( unsigned int i = 0; i < cellSets[1].size(); ++i){
		Cell c = cellSets[1][i];
		if( !is_boundary_cell(dt,c,v) ) continue;
		bcells.push_back(c);
		Point cen = CGAL::centroid(c->vertex(0)->point(), c->vertex(1)->point(), c->vertex(2)->point(), c->vertex(3)->point() );
		bgroup.push_back(cen);
	}
	double distance = 10000000000;
	int origin = 0, destiny = 0;
	for( unsigned int i = 0; i < agroup.size(); ++i){
		Vector a = agroup[i] - CGAL::ORIGIN;
		for( unsigned int j = 0; j < bgroup.size(); ++j){
			Vector b = bgroup[j] - CGAL::ORIGIN;
			double current = a*b;
			if( current < distance ){
				origin = i;
				destiny = j;
			}
		}
	}	

	bool aux = true;
	Cell node = acells[origin];
	Cell end_node = bcells[destiny];
	Cell next = acells[origin];
	node->info().m_path_state = true;

	bool is_path = true;

	while(aux){
		node = next;
		int cindex = node->index(v);
		double dif = 1000000000;

		Vector a;
		Vector b = bgroup[destiny] - CGAL::ORIGIN;
		for( unsigned int i = 0; i < 3; ++i){
			int nindex = dt.vertex_triple_index(cindex,i);
			Cell n_node = node->neighbor(nindex);//compare distances HERE
			//if( n_node->info().m_path_state == true || n_node->info().m_label  ==  acells[origin]->info().m_label ) continue;
			if( n_node->info().m_path_state || dt.is_infinite(n_node) ) continue;
			if( std::find( bcells.begin(), bcells.end(), n_node) != bcells.end() ){
				aux = false;
				break;
			}
			Point cen = CGAL::centroid(n_node->vertex(0)->point(), n_node->vertex(1)->point(), n_node->vertex(2)->point(), n_node->vertex(3)->point() );
			a = cen - CGAL::ORIGIN;
			double current = a*b;
			if( current < dif ){
				dif = current;
				next = n_node;
			}
		} 

		if(aux){
			if( next== node) {
				is_path = false;
				break;
			}
			next->info().m_path_state = true;
			path.push_back(next);
			//HERE found the new destiny
			/*destiny = new_destiny(dt,bgroup, CGAL::ORIGIN + a );
			cout<<"new destiny: "<<destiny<<endl;*/

			/*cout<<"label"<<path.back()->info().m_label<<endl;
			drawMesh(dt);
			drawCells(dt,path);
			cout<<"drawing inside join"<<endl;
			getchar();*/
		}
	}

	for( unsigned int i = 0; i < cellSets[0].size(); ++i){// i can take jus the boundary cells
		Cell c = cellSets[0][i];
		c->info().m_path_state = false;
	}

	for( unsigned int i = 0; i < path.size(); ++i){
		path[i]->info().m_path_state = false;
	}

	if(is_path)
		return true;
	else
		return false;
}

bool path_to_join_2(Delaunay &dt, vector<vector<Cell>> &cellSets,Vertex &v,vector<Cell> &path){
	path.clear();
	if( cellSets.size() != 2 ) return false;		
	vector<Point> agroup,bgroup;
	vector<Facet> afacets,bfacets;
	for( unsigned int i = 0; i < cellSets[0].size(); ++i){// i can take jus the boundary cells
		Cell c = cellSets[0][i];
		c->info().m_path_state = true;
		int label = c->info().m_label;
		int vindex = c->index(v);
		for( unsigned int i = 0; i < 3; ++i){
			int index = dt.vertex_triple_index(vindex,i);
			Cell nc = c->neighbor(index);
			if( !dt.is_infinite(nc) && nc->info().m_label != label){
				Facet f = Facet(nc,index);
				afacets.push_back(f);
				Point cen = CGAL::centroid( dt.triangle(f) );
				agroup.push_back(cen);
			}
		} 
	}

	for( unsigned int i = 0; i < cellSets[1].size(); ++i){// i can take just the boundary cells
		Cell c = cellSets[1][i];
		int label = c->info().m_label;
		int vindex = c->index(v);
		for( unsigned int i = 0; i < 3; ++i){
			int index = dt.vertex_triple_index(vindex,i);
			Cell nc = c->neighbor(index);
			if( !dt.is_infinite(nc) && nc->info().m_label != label){
				Facet f = Facet(nc,index);
				bfacets.push_back(f);
				Point cen = CGAL::centroid( dt.triangle(f) );
				bgroup.push_back(cen);
			}
		} 
	}

	double distance = 10000000000;
	int origin = 0, destiny = 0;
	for( unsigned int i = 0; i < agroup.size(); ++i){
		Vector a = agroup[i] - CGAL::ORIGIN;
		for( unsigned int j = 0; j < bgroup.size(); ++j){
			Vector b = bgroup[j] - CGAL::ORIGIN;
			double current = a*b;
			if( current < distance ){
				origin = i;
				destiny = j;
			}
		}
	}	

	bool aux = true;
	Cell node = afacets[origin].first;
	Cell end_node = bfacets[destiny].first;

	node->info().m_path_state = true;
	path.push_back(node);
	Cell next_cell = node; 

	Vector b = bgroup[destiny] - CGAL::ORIGIN;

	bool is_path = true;
	while(aux){
		Cell c = next_cell;
		int vindex = c->index(v);
		Facet next_facet;
		double dif = 1000000000;
		for( unsigned int i = 0; i < 3; ++i){// iterate all the facets
			int index = dt.vertex_triple_index(vindex,i);
			Cell nc = c->neighbor(index);
			if( std::find( cellSets[1].begin(), cellSets[1].end(), nc) != cellSets[1].end() ){
				aux = false;
				break;
			}
			if( !dt.is_infinite(nc) && !nc->info().m_path_state){
				Facet f = Facet(nc,index);
				bfacets.push_back(f);
				Point cen = CGAL::centroid( dt.triangle(f) );
				Vector a = cen - CGAL::ORIGIN;				
				double current = a*b;
				if( current < dif ){
					dif = current;
					next_facet = f;
					next_cell = next_facet.first;
				}
			}
		} 

		if(aux){	
			if( next_cell == node) {
				is_path = false;
				break;
			}
			next_cell->info().m_path_state = true;
			path.push_back(next_cell);

			/*cout<<"label"<<path.back()->info().m_label<<endl;
			drawMesh(dt);
			drawCells(dt,path);
			cout<<"drawing inside join"<<endl;
			getchar();*/
		}
	}

	for( unsigned int i = 0; i < cellSets[0].size(); ++i){// i can take jus the boundary cells
		Cell c = cellSets[0][i];
		c->info().m_path_state = false;
	}
	for( unsigned int i = 0; i < path.size(); ++i){
		path[i]->info().m_path_state = false;
	}
	if(is_path)
		return true;
	else
		return false;
}

bool path_to_join_3(Delaunay &dt, vector<vector<Cell>> &cellSets,Vertex &v,vector<Cell> &path){
	path.clear();
	if( cellSets.size() != 2 ) return false;	
	vector<Point> agroup,bgroup;
	vector<Cell> acells,bcells;


	for( unsigned int i = 0; i < cellSets[0].size(); ++i){// i can take jus the boundary cells
		Cell c = cellSets[0][i];
		//c->info().m_path_state = true;
		if( !is_boundary_cell(dt,c,v) ) continue;		
		acells.push_back(c);
		Point cen = CGAL::centroid(c->vertex(0)->point(), c->vertex(1)->point(), c->vertex(2)->point(), c->vertex(3)->point() );
		agroup.push_back(cen);
	}
	for( unsigned int i = 0; i < cellSets[1].size(); ++i){
		Cell c = cellSets[1][i];
		if( !is_boundary_cell(dt,c,v) ) continue;
		bcells.push_back(c);
		Point cen = CGAL::centroid(c->vertex(0)->point(), c->vertex(1)->point(), c->vertex(2)->point(), c->vertex(3)->point() );
		bgroup.push_back(cen);
	}

	double distance = 10000000000;
	int origin = 0, destiny = 0;
	for( unsigned int i = 0; i < agroup.size(); ++i){
		Vector a = agroup[i] - CGAL::ORIGIN;
		for( unsigned int j = 0; j < bgroup.size(); ++j){
			Vector b = bgroup[j] - CGAL::ORIGIN;
			Vector c = a - b;	
			double current = c*c;
			if( current < distance ){
				distance = current;
				origin = i;
				destiny = j;
			}
		}
	}	

	
	Cell c = acells[origin];
	//Cell c = afacets[origin].first;
	Cell next = c;
	Point p = bgroup[destiny];


	bool is_path = false;
	bool try_next_cell = true;
	while(try_next_cell){
		try_next_cell = false;
		int cindex = c->index(v);
		const Point* pts[4] = { &(c->vertex(0)->point()),&(c->vertex(1)->point()),&(c->vertex(2)->point()),&(c->vertex(3)->point()) };
		for( unsigned int i = 0; !try_next_cell && i != 3; ++i){
			int index = dt.vertex_triple_index(cindex,i);
			Cell next = c->neighbor( index );
			if( next->info().m_path_state || dt.is_infinite(next) ) continue;// visited cell or infinite cell
			//if( std::find( bcells.begin(), bcells.end(), next) != bcells.end() )	{is_path = true; break;} //we arrive to destiny

			const Point* backup = pts[index];
			pts[index] = &p;
			if ( orientation(*pts[0], *pts[1], *pts[2], *pts[3]) != CGAL::NEGATIVE ){
				pts[index] = backup;
			}
			else{
				next->info().m_path_state = true;
				c = next;
				path.push_back(next);
				//try_next_cell = true;
				
				for( unsigned int j = 0 ; j < 4; ++j){
					Cell nc = next->neighbor(j);
					if( std::find( bcells.begin(), bcells.end(), nc) != bcells.end() )	{is_path = true; } //we arrive to destiny
				}

				if(is_path)
					break;
				else
					try_next_cell = true;;
			}
		}
	}

	for( unsigned int i = 0; i < cellSets[0].size(); ++i){
		cellSets[0][i]->info().m_path_state = false;
	}

	for( unsigned int i = 0; i < path.size(); ++i){
		if( !path[i]->info().m_path_state ){
			path.erase( path.begin() + i);
			--i;
		}
		else	path[i]->info().m_path_state = false;
	}

	if( is_path) return true;
	else	return false;
}

bool path_to_join(Delaunay &dt, vector<vector<Cell>> &cellSets,Vertex &v,vector<Cell> &path){
	path.clear();
	if( cellSets.size() != 2 ) return false;	
	vector<Point> agroup,bgroup;
	vector<Cell> acells,bcells;

	/*drawMesh(dt);
	drawCells(dt,cellSets[0]);
	getchar();*/

	vector<Facet> afacets;
	for( unsigned int i = 0; i < cellSets[0].size(); ++i){// i can take jus the boundary cells
		Cell c = cellSets[0][i];
		int label = c->info().m_label;
		int vindex = c->index(v);
		for( unsigned int j = 0; j < 3; ++j){
			int index = dt.vertex_triple_index(vindex,j);
			Cell nc = c->neighbor(index);
			if( !dt.is_infinite(nc) && nc->info().m_label != label){
				Facet f = Facet(c,index);
				afacets.push_back(f);
				Point cen = CGAL::centroid( dt.triangle(f) );
				agroup.push_back(cen);
			}
		} 
	}

	for( unsigned int i = 0; i < cellSets[1].size(); ++i){// i can take just the boundary cells
		Cell c = cellSets[1][i];
		bcells.push_back(c);

		int label = c->info().m_label;
		int vindex = c->index(v);
		for( unsigned int j = 0; j < 3; ++j){
			int index = dt.vertex_triple_index(vindex,j);//index of c
			Cell nc = c->neighbor(index);
			if( !dt.is_infinite(nc) && nc->info().m_label != label){
				Facet f = Facet(c,index);
				Point cen = CGAL::centroid( dt.triangle(f) );
				bgroup.push_back(cen);
			}
		} 
	}
	double distance = 10000000000;
	int origin = 0, destiny = 0;
	for( unsigned int i = 0; i < agroup.size(); ++i){
		Vector a = agroup[i] - CGAL::ORIGIN;
		for( unsigned int j = 0; j < bgroup.size(); ++j){
			Vector b = bgroup[j] - CGAL::ORIGIN;
			Vector c = a - b;			
			double current = c*c;	
			if( current < distance ){
				distance = current;
				origin = i;
				destiny = j;
			}
		}
	}	

	/*vector<Point> pts;
	pts.push_back(agroup[origin]);
	pts.push_back(bgroup[destiny]);
	drawMesh(dt,pts);
	getchar();*/


	Cell c = afacets[origin].first;
	Cell next = c;
	Point p = bgroup[destiny];

	bool is_path = false;
	bool try_next_cell = true;
	while(try_next_cell){
		try_next_cell = false;
		int cindex = c->index(v);
		const Point* pts[4] = { &(c->vertex(0)->point()),&(c->vertex(1)->point()),&(c->vertex(2)->point()),&(c->vertex(3)->point()) };
		for( unsigned int i = 0; !try_next_cell && i != 3; ++i){
			int index = dt.vertex_triple_index(cindex,i);
			Cell next = c->neighbor( index );
			if( next->info().m_path_state || dt.is_infinite(next) ) continue;// visited cell or infinite cell
			//if( std::find( bcells.begin(), bcells.end(), next) != bcells.end() )	{is_path = true; break;} //we arrive to destiny
			const Point* backup = pts[index];
			pts[index] = &p;
			if ( orientation(*pts[0], *pts[1], *pts[2], *pts[3]) != CGAL::NEGATIVE )
				pts[index] = backup;
			else{
				next->info().m_path_state = true;
				c = next;
				path.push_back(next);			
				
				for( unsigned int j = 0 ; j < 3; ++j){
					int index = dt.vertex_triple_index(next->index(v),j);
					Cell nc = next->neighbor(index);
					if( std::find( bcells.begin(), bcells.end(), nc) != bcells.end() )	{is_path = true; } //we arrive to destiny
				}

				if(is_path)
					break;
				else
					try_next_cell = true;;
			}
		}
	}

	for( unsigned int i = 0; i < cellSets[0].size(); ++i)
		cellSets[0][i]->info().m_path_state = false;

	for( unsigned int i = 0; i < path.size(); ++i){
		if( !path[i]->info().m_path_state ){
			path.erase( path.begin() + i);
			--i;
		}
		else	path[i]->info().m_path_state = false;
	}

	/*cout<<"size path: "<<path.size()<<endl;
	drawMesh(dt);
	drawCells(dt,path);
	getchar();*/


	if( is_path) return true;
	else	return false;
}