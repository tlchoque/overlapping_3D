#include <queue>

struct Ccomponent{
	Vertex v;
	double k;
	int label;
	Ccomponent( Vertex vi,int labeli,double ki ){
		v = vi;
		label = labeli;
		k = ki;
	}
	bool operator<(const Ccomponent& rhs) const {
        return k < rhs.k;
    }
};

bool comparator2(const Ccomponent &a, const Ccomponent &b){
    return (a.k < b.k);
}

bool is_feature(Delaunay &dt, vector<Cell> &cells,Vertex v, int label){
	vector<Vertex> vertices;
	vertices = ring_vertices(dt,cells,v,label);
	double mcv = mean_curvature2(vertices,v);
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

	if( mcv*sum > 0)// the vector points up
		return true;
	return false;
}

void mesh_curvature(Delaunay &dt, vector<Vertex> &boundary, priority_queue<Ccomponent> &peaks, double &total_curvature){	
	for(unsigned int i = 0 ; i < boundary.size(); ++i){
		vector<Ccomponent> components;
		Vertex v = boundary[i];
		bool add = true;
		//cout<<"stop 0 mesh crvature "<<v->point()<<endl;
		for(unsigned int j = 0; j < v->info().m_regions_around.size(); ++j){			
			int label = v->info().m_regions_around[j];
			//cout<<"stop 1 mesh crvature: "<<label<<endl;
			vector<vector<Cell>> cellSets;
			set_cells_around_singular_vertex(dt,v,label,cellSets);
			vector<Vertex> vertices = ring_vertices(dt,cellSets[0],v,label);
			double k = mean_curvature2(vertices,v);
			if( is_feature(dt,cellSets[0],v,label) ){
				add = false;
				v->info().m_feature = true;
				break;
			}
			if( k > 0 )
				components.push_back(Ccomponent(v,label,k ) ); 
		}
		if( add ){
			for(unsigned int l = 0; l < components.size(); ++l){
				Ccomponent c = components[l];
				peaks.push(c); 
			}
		}
	}
}

bool is_relabeled(vector<Cell> &cells){
	for(unsigned int i = 0; i < cells.size(); ++i){
		Cell f = cells[i];
		if( f->info().m_relabeled )
			return true;
	}
	return false;
}

void add_component(Delaunay &dt,Vertex v,priority_queue<Ccomponent> &peaks){
	if( v->info().m_feature ) return;	
	bool add = true;
	vector<Ccomponent> components;
	for(unsigned int i = 0; i < v->info().m_regions_around.size(); ++i){
		int label = v->info().m_regions_around[i];
		vector<vector<Cell>> cellSets;
		set_cells_around_singular_vertex(dt,v,label,cellSets);
		vector<Vertex> vertices = ring_vertices(dt,cellSets[0],v,label);
		double k = mean_curvature2(vertices,v);
		if( k > 0 && !is_relabeled(cellSets[0]) )
			components.push_back(Ccomponent(v,label,k ) ); 
	}
	if( add ){
		for(unsigned int l = 0; l < components.size(); ++l){
			Ccomponent c = components[l];
			peaks.push(c); 
		}
	}
}

bool right_neigborhood(Delaunay &dt,Vertex v,vector<Cell> &cells,int label){
	for(unsigned int i = 0; i< cells.size(); ++i){
		Cell f =cells[i];
		Cell n = f->neighbor( f->index(v) );
		if( n->info().m_label == label  )
			return false;
	}
	return true;
}

bool has_corner_vertex(vector<Vertex> &vertices){
	for(unsigned int i = 0; i< vertices.size(); ++i){
		Vertex v = vertices[i];
		if( v->info().is_corner() )
			return true;
	}
	return false;
}

vector<Cell> cells_complement(Delaunay &dt,Vertex v, int label){
	vector<Cell> comp;
	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		Cell c = *it;
		if( c->info().m_label != label )
			comp.push_back(c);
	}
	return comp;
}

bool edge_restriction(Delaunay &dt,vector<Cell> &cells, Vertex v, int label){
	vector<Cell> comp = cells_complement(dt,v,label);
	double c1 = larger_edge(dt,cells);
	double c2 = larger_edge(dt,comp);
	if( c1/c2 < 1.5)
		return false;
	else
		return true;
}

bool is_possible_to_relabel(Delaunay &dt, Ccomponent c,vector<Vertex> &neighbors){
	//is label info of vertex 
	if( !( std::find(c.v->info().m_regions_around.begin(), c.v->info().m_regions_around.end(), c.label ) != c.v->info().m_regions_around.end() ) ) return false;
	vector<vector<Cell>> cellSets;
	set_cells_around_singular_vertex(dt,c.v,c.label,cellSets);
	if( is_relabeled(cellSets[0]) ) return false;

	int new_label = new_label_for_relabeling(dt,cellSets[0],c.v,c.label); 
	if( !right_neigborhood(dt,c.v,cellSets[0],new_label) ) {
		return false;
	}
	vector<Vertex> vertices;
	neighbor_vertices(dt,cellSets[0],vertices,c.v);
	vertices.push_back(c.v);
	double sum_before=0;
	for(unsigned int i = 0; i< vertices.size(); ++i){
		Vertex v = vertices[i];
		if( v->info().is_boundary() ){
			set_stars_vertex(dt,v);
			vector<double> curvatures;
			sum_before+=vertex_curvature(v,curvatures);
		}
	}

	for(unsigned int i = 0; i< cellSets[0].size(); ++i)//relabel
		cellSets[0][i]->info().m_label = new_label;

	double sum_after=0;
	for(unsigned int i = 0; i< vertices.size(); ++i){
		Vertex v = vertices[i];
		update_vertex_info(dt,v);
		if( v->info().is_boundary() && !v->info().is_singular_2() ){
			set_stars_vertex(dt,v);
			vector<double> curvatures;
			sum_after+=vertex_curvature(v,curvatures);
		}
	}	

	//cout<<sum_before<<" * "<<sum_after<<endl;
	if( has_corner_vertex(vertices) || edge_restriction(dt,cellSets[0],c.v,c.label ) || produce_singularities(dt,vertices) || sum_before/sum_after < 3.5 ){//4 is for hyena  //3 for carp
		for(unsigned int i = 0; i< cellSets[0].size(); ++i)
			cellSets[0][i]->info().m_label = c.label;
				
		for(unsigned int i = 0; i< vertices.size(); ++i)
			update_vertex_info(dt,vertices[i]);
		return false;
	}
	else{
		for(unsigned int i = 0; i< cellSets[0].size(); ++i)
			cellSets[0][i]->info().m_relabeled = true;

		c.v->info().m_state=4;//why

		/*for(unsigned int i = 0; i< vertices.size(); ++i){
			Vertex v = vertices[i];
			if( v->info().m_state == 1 && v->info().is_boundary() )
				neighbors.push_back(v);
		}*/
		return true;
	}
}

void relabeling_to_smooth(Delaunay &dt){
	vector<Vertex> boundary;	
	boundary_vertices(dt,boundary);
	double total_curvature;
	std::priority_queue<Ccomponent> peaks;
	cout<<"start mesh crvature"<<endl;
	mesh_curvature(dt,boundary,peaks,total_curvature);
	cout<<"start relabeling"<<endl;
	while( peaks.size() !=0 ){
		Ccomponent c = peaks.top();
		peaks.pop();
		if( c.v->info().is_singular_2() || c.v->info().is_corner() ) continue;

		vector<Vertex> neighbors;		
		//vector<Point> ptos;
		//cout<<" before "<<c.label<<endl;
		//ptos.push_back(c.v->point());
		/*drawMesh(dt,ptos),
		getchar();*/

		//cout<<" stop 0"<<endl;
		if( is_possible_to_relabel(dt,c,neighbors) ){
			//cout<<" after"<<endl;
			/*drawMesh(dt,ptos),
			getchar();*/

			//for( unsigned int i = 0 ; i < neighbors.size(); ++i ){
			//	//cout<<" stop 1"<<endl;
			//	add_component(dt,neighbors[i],peaks);
			//}
		}
	}
	cout<<"end relabeling"<<endl;
	restart_vertex_state(dt);	
}