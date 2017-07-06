
void singularities(Delaunay &dt, vector<Point> &points,vector< pair<Point,int> > &centroids, Vertex &v	){
	load_points2(points);
	dt.insert(points.begin(),points.end());
	forManifoldVertex(dt,centroids);	
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) { 
		if( vi->point().x() == 25 && vi->point().y() == 21 && vi->point().z() == 30 )
			v=vi;
	} 
}

void test_curvature(Delaunay &dt){
	vector< pair<Point,int> > centroids;
	vector<pair<Point,Point>> edges;
	vector<Point> points;
	vector<Facet> facets;	
	Vertex v;
	singularities(dt,points,centroids,v);
	vector<vector<Cell>> cellSets;
	set_cells_around_singular_vertex(dt,v,0,cellSets);

	/*cout<<cellSets.size()<<endl;
	cout<<cellSets[0].size()<<endl;*/

	vector<Vertex> vertices;
	vertices = ring_vertices(dt,cellSets[0],v,0);
	vector<Point> ring_points;
	
	//vertices_to_points(vertices,ring_points);

	/*for(unsigned int i = 0; i< vertices.size(); ++i)
		cout<<ring_points[i]<<endl;*/

	//double k1,k2;

	double curvature = mean_curvature2(vertices,v);
	//double curvature = mean_curvature(vertices,v);
	double gaussian = gaussian_curvature(vertices,v);
	//principal_curvatures(vertices,v,k1,k2);
	Vector n1 = normal_curvature(vertices,v);

	cout<<"1 mean: "<<curvature<<endl;
	cout<<"1 gauss: "<<gaussian<<endl;	
	
	set_cells_around_singular_vertex(dt,v,1,cellSets);
	vertices = ring_vertices(dt,cellSets[0],v,1);
	curvature = mean_curvature2(vertices,v);
	//curvature = mean_curvature(vertices,v);
	gaussian = gaussian_curvature(vertices,v);
	//principal_curvatures(vertices,v,k1,k2);
	Vector n2 = normal_curvature(vertices,v);


	cout<<"2 mean: "<<curvature<<endl;
	cout<<"3 gauss: "<<gaussian<<endl;

	ring_points.push_back(v->point() );
	ring_points.push_back( CGAL::ORIGIN + (200*n1 + ( v->point() - CGAL::ORIGIN ) ) ); // for grey
	ring_points.push_back( CGAL::ORIGIN + (200*n2 + ( v->point() - CGAL::ORIGIN ) ) ); // for red
	drawMesh(dt,ring_points);


	getchar();
}

void test_for_delaunay(Delaunay &dt){
	vector< pair<Point,int> > centroids;
	vector<pair<Point,Point>> edges;
	vector<Point> points;
	vector<Facet> facets;	
	Vertex v;
	singularities(dt,points,centroids,v);
}