vector<Facet> extract_facets(Delaunay &dt,vector<vector<Cell>> &cellSets, Vertex v) {
	vector<Facet> facets;
	for (unsigned int i = 0; i < cellSets.size(); ++i) {
		for (unsigned int j = 0; j < cellSets[i].size(); ++j) {
			Cell c = cellSets[i][j];
			int idx = c->index(v);
			for (unsigned int k = 0; k < 3; ++k) {
				int facet_idx = dt.vertex_triple_index(idx,k); 
				if ( c->neighbor(facet_idx)->info().m_label != c->info().m_label ) {
					Facet e = Facet( c, facet_idx);
					facets.push_back(e);
				}
			}
		}
	}
	return facets;
}

vector<Facet> comp_facets(Delaunay &dt, vector<Cell> &cells, Vertex v) {
	vector<Facet> facets;
	for (unsigned int j = 0; j < cells.size(); ++j) {
		Cell c = cells[j];
		int idx = c->index(v);
		for (unsigned int k = 0; k < 3; ++k) {
			int facet_idx = dt.vertex_triple_index(idx, k);
			if (c->neighbor(facet_idx)->info().m_label != c->info().m_label) {
				Facet e = Facet(c, facet_idx);
				facets.push_back(e);
			}
		}
	}
	return facets;
}

double solid_angle(Vector a, Vector b, Vector c){
	double a_n = norm(a);
	double b_n = norm(b);
	double c_n = norm(c);
	double det = CGAL::determinant(a,b,c);
	//return 2*atan2( det, a_n*b_n*c_n  + c_n*a*b + a_n*b*c + b_n*c*a );
	return -2 * atan2(det, a_n*b_n*c_n + c_n*a*b + a_n*b*c + b_n*c*a);
}

double point_w( Delaunay &dt, Point p, vector<Facet> &facets, Vertex v){
	double total = 0;
	for (unsigned int i = 0; i < facets.size(); ++i) {
		Facet f = facets[i];
		int v_idx = f.second;
		Point A = f.first->vertex( dt.vertex_triple_index(v_idx,0) )->point();
		Point B = f.first->vertex( dt.vertex_triple_index(v_idx,1) )->point();
		Point C = f.first->vertex( dt.vertex_triple_index(v_idx,2) )->point();
		Vector a = A - p;
		Vector b = B - p;
		Vector c = C - p;
		total += solid_angle(a,b,c);

		/*vector<Point> pts;
		pts.push_back(p);
		pts.push_back(A);
		pts.push_back(B);
		pts.push_back(C);
		cout<<"w: "<<solid_angle(a,b,c)<<endl;
		drawMesh(dt,facets,pts);
		getchar();*/
	}
	//return total/(4*PI);
	cout << "solid angle: " << total << endl;
	return total;
}

double criteria_w(Delaunay &dt, vector<Cell> &cells, vector<Facet> &facets, Vertex v ) {
	double w_average = 0;
	vector<Facet> component_facets = comp_facets(dt,cells, v);
	for (unsigned int i = 0; i < cells.size(); ++i) {
		Cell c = cells[i];
		double volume = CGAL::volume( c->vertex(0)->point(), c->vertex(1)->point(), c->vertex(2)->point(), c->vertex(3)->point() ) ;
		Point p = CGAL::centroid( c->vertex(0)->point(), c->vertex(1)->point(), c->vertex(2)->point(), c->vertex(3)->point() );
		double w = point_w(dt,p, facets,v);
		//double w = point_w(dt, p, facets, v);

		w_average += w * volume;
		//w_average += w;
	}
	cout << "w_avergae: " << w_average << endl;
	return 1/w_average;
	//return cells.size() / w_average;
}

void sort_by_criteria_w(Delaunay &dt, vector<vector<Cell>> &cellSets, Vertex v, int &min, int &max, bool &changeable) {
	vector<Facet> facets = extract_facets(dt,cellSets, v);
	
	vector<double> ratio( cellSets.size() );
	for (unsigned int i = 0; i < cellSets.size(); ++i) {
		ratio[i] = criteria_w(dt, cellSets[i], facets, v);
		cout << "ratio: " << ratio[i] << endl;
	}
	max = 0; min = 1;
	if (ratio[min] > ratio[max]) {
		max = 1;
		min = 0;
	}
	if (ratio[max] / ratio[min] < 1)//for chest is different 2.65 3
		changeable = true;
	else
		changeable = false;
}

