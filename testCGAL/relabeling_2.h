#include <queue>

struct Celement{
	Vertex v;
	double k;
	int index;
	int label;
	Celement( Vertex vi,int labeli,double ki,int idx ){
		v = vi;
		label = labeli;
		k = ki;
		index =idx;
	}
	bool operator<(const Celement& rhs) const {
        return k < rhs.k;
    }
};

void update_cuvatures(Delaunay &dt, Vertex v){
	// change for 2 label case
	vector<double>  &curvatures =v->info().m_curvatures;
	curvatures.clear();

	if( v->info().m_regions_around.size() == 2 ){
		double k = mean_curvature2(v->info().m_sa_stars[0],v);
		curvatures.push_back(k);
		curvatures.push_back(-k);
	}
	else{
		for(unsigned int i = 0 ; i < v->info().m_sa_stars.size(); ++i){
			double k = mean_curvature2(v->info().m_sa_stars[i],v);
			curvatures.push_back(k);
		}
	}	
}

void relabeling_initialization(Delaunay &dt, vector<Vertex> &boundary){
	for(unsigned int i = 0 ; i < boundary.size(); ++i){
		Vertex v = boundary[i];
		//cout<<"stop 0 relabeling_initialization"<<endl;
		if(v->info().is_corner() ) continue;
		set_stars_vertex(dt,v);
		update_cuvatures(dt,v);
		/*if(v->info().is_corner() ){
			cout<<"m_sa_stars size : "<< v->info().m_sa_stars.size()<<endl;
			for(unsigned int i = 0; i < v->info().m_sa_starss.size(); ++i){
				vector<Point> pts = vtop(v->info().m_sa_stars[i]);
				pts.push_back(v->point());
				cout<<"pts size : "<< pts.size()<<" "<< v->point()<<endl;
				drawMesh(dt,pts);
				getchar();
			}
		}*/

		//cout<<"stop 1 relabeling_initialization"<<endl;
	}
}

double curvature_by_label(Delaunay &dt, Vertex v, int label){
	int index = -1;
	for(unsigned int i = 0 ; i< v->info().m_regions_around.size();  ++i ){
		if( label == v->info().m_regions_around[i] ){
			index = i;
		}
	}
	if( index != -1 )
		return v->info().m_curvatures[index];
	else {
		cout<<"problem, size labels "<< v->info().m_regions_around.size()<<endl;

		if( v->info().is_corner() ) cout<<"and is corner"<<endl;

		/*vector<Point> pts;
		pts.push_back(v->point());
		drawMesh(dt,pts);
		getchar();*/

		return -1;
	}
}

bool is_feature_2(Delaunay &dt, Vertex v, int label,int index){
	return true;

	vector<Vertex> vertices = v->info().m_sa_stars[index];	
	double mcv = v->info().m_curvatures[index];
	double sum = 0;
	for(unsigned int i = 0; i < vertices.size(); ++i){
		Vertex w = vertices[i];
		if(w->info().is_corner() ) continue;
		//if( w->info().is_singular_2() || w->info().is_corner() ) continue; // how to include the vertices ni the corner
		if( w->info().is_singular_2() ) continue;
		double k = curvature_by_label(dt,w,label);
		/*if(k == -1){
			vector<Point> pts = vtop(vertices);
			pts.push_back(v->point());
			cout<<"label and size : "<<label<<" "<< pts.size()<<" "<< v->point()<<endl;
			drawMesh(dt,pts);
			getchar();
		}*/
			 
		sum = sum + k;
	}
	if( mcv*sum > 0)// the vector points up
		return true;
	return false;
}

void select_elements_1(Delaunay &dt, Vertex v, priority_queue<Celement> &peaks){//for 2 region vertices
	int index = 0;
	if( v->info().m_curvatures[1] > v->info().m_curvatures[0] )	index =1;
	if( is_feature_2(dt,v,v->info().m_regions_around[index],index) ){
		v->info().m_feature = true;
		peaks.push(Celement(v,v->info().m_regions_around[index],v->info().m_curvatures[index],index) );
	}
}

void select_elements_2(Delaunay &dt, Vertex v, priority_queue<Celement> &peaks){//for +2 region vertices
	vector<Celement> elements;
	bool add = true;
	for(unsigned int i = 0; i < v->info().m_regions_around.size(); ++i){	
		if( is_feature_2(dt,v,v->info().m_regions_around[i],i) ){
			add = false;
			v->info().m_feature = true;
			break;
		}
		if( v->info().m_curvatures[i] > 0 )
			elements.push_back( Celement(v,v->info().m_regions_around[i],v->info().m_curvatures[i],i) ); 
	}
	if( add ){
		for(unsigned int l = 0; l < elements.size(); ++l){
			Celement c = elements[l];
			peaks.push(c); 
		}
	}
}

void available_vertices(Delaunay &dt, vector<Vertex> &boundary, priority_queue<Celement> &peaks, double &total_curvature){	
	for(unsigned int i = 0 ; i < boundary.size(); ++i){
		Vertex v = boundary[i];
		if(v->info().is_corner() ) continue;
		if( v->info().m_regions_around.size() == 2)
			select_elements_1(dt,v,peaks);
		else
			select_elements_2(dt,v,peaks);
	}
}

double vertex_curvature_2(Vertex &v){
	double total=0;
	for( unsigned int i = 0 ; i < v->info().m_curvatures.size(); ++i){
		if( v->info().m_regions_around[i] == HOLE ) continue;
		total += abs(v->info().m_curvatures[i]);
	}
	return total;
}

double vertex_curvature_3(Vertex &v, int & count){
	double total=0;
	for( unsigned int i = 0 ; i < v->info().m_curvatures.size(); ++i){
		if( v->info().m_regions_around[i] == HOLE ) continue;
		++count;
		total += abs(v->info().m_curvatures[i]);
	}
	return total;
}

bool is_possible_to_relabel_2(Delaunay &dt, Celement c,vector<Vertex> &neighbors){
	if( !( std::find(c.v->info().m_regions_around.begin(), c.v->info().m_regions_around.end(), c.label ) != c.v->info().m_regions_around.end() ) ) return false;
	vector<vector<Cell>> cellSets;
	set_cells_around_singular_vertex(dt,c.v,c.label,cellSets);
	if( is_relabeled(cellSets[0]) ) return false;

	int new_label = new_label_for_relabeling(dt,cellSets[0],c.v,c.label); 
	if( !right_neigborhood(dt,c.v,cellSets[0],new_label) ) {
		return false;
	}

	vector<Vertex> vertices;// = c.v->info().m_sa_stars[c.index];
	neighbor_vertices(dt,cellSets[0],vertices,c.v);
	vertices.push_back(c.v);

	double sum_before=0;
	int count_before =0;
	for(unsigned int i = 0; i< vertices.size(); ++i){
		Vertex v = vertices[i];
		if( v->info().is_corner() ) continue;
		if( v->info().is_boundary() ){
			sum_before+=vertex_curvature_2(v);
			//sum_before+=vertex_curvature_3(v,count_before);
		}
	}

	for(unsigned int i = 0; i< cellSets[0].size(); ++i)//relabel
		cellSets[0][i]->info().m_label = new_label;

	double sum_after=0;
	int count_after =0;
	for(unsigned int i = 0; i< vertices.size(); ++i){
		Vertex v = vertices[i];
		if(v->info().is_corner() ) continue;
		update_vertex_info(dt,v);		
		if( v->info().is_boundary() && !v->info().is_singular_2() ){
			set_stars_vertex(dt,v);
			update_cuvatures(dt,v);
			sum_after+=vertex_curvature_2(v);
		}
	}	

	//cout<<sum_before<<" * "<<sum_after<<endl;

	//if( has_corner_vertex(vertices) || edge_restriction(dt,cellSets[0],c.v,c.label ) || produce_singularities(dt,vertices) || sum_before/sum_after < 1 ){//4 is for hyena  //3 for carp
	if( has_corner_vertex(vertices) || produce_singularities(dt,vertices) || sum_before/sum_after < 1.5 ){
		for(unsigned int i = 0; i< cellSets[0].size(); ++i)
			cellSets[0][i]->info().m_label = c.label;
				
		for(unsigned int i = 0; i< vertices.size(); ++i){
			Vertex u = vertices[i];
			if(u->info().is_corner() ) continue;
			update_vertex_info(dt,u);
			set_stars_vertex(dt,u);
			update_cuvatures(dt,u);
		}
		return false;
	}
	else{
		//cout<<"relabeled"<<endl;
		for(unsigned int i = 0; i< cellSets[0].size(); ++i)
			cellSets[0][i]->info().m_relabeled = true;

		//c.v->info().m_state=4;//why

		/*for(unsigned int i = 0; i< vertices.size(); ++i){
			Vertex v = vertices[i];
			if( v->info().m_state == 1 && v->info().is_boundary() )
				neighbors.push_back(v);
		}*/
		return true;
	}
}

void relabeling_to_smooth_2(Delaunay &dt){
	mark_boundary(dt);
	vector<Vertex> boundary;		
	boundary_vertices(dt,boundary);
	relabeling_initialization(dt,boundary);
	double total_curvature;
	std::priority_queue<Celement> peaks;

	//measure time
	clock_t begin = clock();

	cout<<"start mesh curvature"<<endl;
	available_vertices(dt,boundary,peaks,total_curvature);
	cout<<"start relabeling"<<endl;
	while( peaks.size() !=0 ){
		Celement c = peaks.top();
		peaks.pop();
		if( c.v->info().is_singular_2() || c.v->info().is_corner() ) continue;
		vector<Vertex> neighbors;	

		/*vector<Point> ptos;
		cout<<" before label and curvature"<<c.label<<" "<<c.k<<endl;
		ptos.push_back(c.v->point());
		drawMesh(dt,ptos),
		getchar();*/

		//cout<<" stop 0"<<endl;
		if( is_possible_to_relabel_2(dt,c,neighbors) ){
			/*cout<<" after"<<endl;
			drawMesh(dt,ptos),
			getchar();*/

			//for( unsigned int i = 0 ; i < neighbors.size(); ++i ){
			//	//cout<<" stop 1"<<endl;
			//	add_component(dt,neighbors[i],peaks);
			//}


			/*restart_vertex_state(dt);
			deque<Vertex> singularities;
			boundary_vertices(dt,boundary);
			singular_vertices(dt,boundary,singularities);
			cout<<"relabeling number of non manifold vertices after "<<singularities.size()<<endl;
			cout<<"relabeling number of vertices "<<dt.number_of_vertices()<<endl;
			if( singularities.size() > 0){
				cout<<singularities[0]->point()<<endl;
				vector<Point> ptos;
				ptos.push_back(c.v->point());
				ptos.push_back(singularities[0]->point());
				drawMesh(dt,ptos),
				getchar();
			}*/

		}
	}
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<<"time relabeling: "<<elapsed_secs<<endl;

	restart_vertex_state(dt);	
}