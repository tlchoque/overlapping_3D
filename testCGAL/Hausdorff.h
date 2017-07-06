

void off_file( Delaunay &dt, const char* filename, string file){
	vector<Point> points;
	vector<Facet> facets;
	int k = 0;
	for(Finite_facets_iterator fi=dt.finite_facets_begin();fi!=dt.finite_facets_end(); fi++){
		if( fi->first->info().m_label == fi->first->neighbor(fi->second)->info().m_label )		continue;
		if( dt.is_infinite(fi->first) || dt.is_infinite(fi->first->neighbor(fi->second)) )		continue;

		Cell c = fi->first;
		int index_facet = fi->second;
		if( fi->first->neighbor(fi->second)->info().m_label <  c->info().m_label ) {
			c = fi->first->neighbor(fi->second);
			index_facet = dt.mirror_index(fi->first, fi->second );
		}
		facets.push_back(Facet(c,index_facet));

		for(unsigned int i = 0; i < 3;++i){
			Vertex v = c->vertex ( dt.vertex_triple_index( index_facet,i) );
			if( v->info().m_index != -1	)  continue;
			v->info().m_index  = k;
			++k;
			points.push_back(v->point() );
		}
	}

	std::ofstream os;
	std::stringstream ss;
	ss<<filename<<"/"<<file<<".off";
	os.open( ss.str() );
	os <<"OFF"<<endl;
	os << points.size()<< ' ' << facets.size()<< ' ' <<0<< "\n";

	for(unsigned int i = 0; i < points.size(); ++i){
		os << points[i]<< "\n";
	}
	for(unsigned int i = 0; i < facets.size(); ++i){
		Facet f = facets[i];
		os <<"3";
		for(unsigned int i = 0; i < 3;++i){
		//for(unsigned int i = 2; i >-1;--i){
			Vertex v = f.first->vertex ( dt.vertex_triple_index( f.second,i) );
			os << ' ' << v->info().m_index;
		}
		os << '\n';
	}

	os.close();

	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) 
		vi->info().m_index = -1;

	cout<<"finsish OFF file"<<endl;
}


void off_file_points( vector<Point> &points, const char* filename, string file){
	std::ofstream os;
	std::stringstream ss;
	ss<<filename<<"/"<<file<<".off";
	os.open( ss.str() );
	os <<"OFF"<<endl;
	os << points.size()<< ' ' << 0<< ' ' <<0<< "\n";
	for(unsigned int i = 0; i < points.size(); ++i){
		os << points[i]<< "\n";
	}
	os.close();
}

void obj_file( Delaunay &dt, const char* filename, string file){
	vector<Point> points;
	vector<Facet> facets;
	int k = 1;
	for(Finite_facets_iterator fi=dt.finite_facets_begin();fi!=dt.finite_facets_end(); fi++){
		if( fi->first->info().m_label == fi->first->neighbor(fi->second)->info().m_label )		continue;
		if( dt.is_infinite(fi->first) || dt.is_infinite(fi->first->neighbor(fi->second)) )		continue;

		Cell c = fi->first;
		int index_facet = fi->second;
		if( fi->first->neighbor(fi->second)->info().m_label <  c->info().m_label ) {
			c = fi->first->neighbor(fi->second);
			index_facet = dt.mirror_index(fi->first, fi->second );
		}
		facets.push_back(Facet(c,index_facet));

		for(unsigned int i = 0; i < 3;++i){
			Vertex v = c->vertex ( dt.vertex_triple_index( index_facet,i) );
			if( v->info().m_index != -1	)  continue;
			v->info().m_index  = k;
			++k;
			points.push_back(v->point() );
		}
	}

	std::ofstream os;
	std::stringstream ss;
	ss<<filename<<"/"<<file<<".obj";
	os.open( ss.str() );
	//os <<"OFF"<<endl;
	//os << points.size()<< ' ' << facets.size()<< ' ' <<0<< "\n";

	for(unsigned int i = 0; i < points.size(); ++i){
		os <<"v "<< points[i]<<" "<<1<<" "<<1<<" "<<0<< "\n";
	}
	for(unsigned int i = 0; i < facets.size(); ++i){
		Facet f = facets[i];
		os <<"f";
		for(unsigned int i = 0; i < 3;++i){
		//for(unsigned int i = 2; i >-1;--i){
			Vertex v = f.first->vertex ( dt.vertex_triple_index( f.second,i) );
			os << ' ' << v->info().m_index<<"//"<<v->info().m_index;
		}
		os << '\n';
	}

	os.close();

	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) 
		vi->info().m_index = -1;

	cout<<"finsish obj file"<<endl;
}