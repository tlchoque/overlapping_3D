
Vector average_without_point( vector<Vertex> &star, Vector x ){
	Vector m = Vector(0,0,0);
	for(unsigned int i = 0; i < star.size(); ++i)
		m = m + (star[i]->point() - CGAL::ORIGIN );
	return m - x;
}

Vector star_vertex_average(Delaunay &dt, Vertex v){
	list<Vertex> incident;
	dt.finite_incident_vertices(v,back_inserter(incident));	
	Vector avg = Vector(0,0,0);
	int n = 0;
	for(std::list<Vertex>::iterator it = incident.begin(); it != incident.end(); it++){
		avg = avg + ( (*it)->point() - CGAL::ORIGIN );
		++n;
	}
	return avg/n;
}

Vector average_around_vertex(Delaunay &dt, Vertex v, Vector x , int &n){
	list<Vertex> incident;
	dt.finite_incident_vertices(v,back_inserter(incident));	
	Vector avg = Vector(0,0,0);
	n = 0;
	for(std::list<Vertex>::iterator it = incident.begin(); it != incident.end(); it++){
		if( !(*it)->info().is_boundary() ) continue;
		avg = avg + ( (*it)->point() - CGAL::ORIGIN );
		++n;
	}
	return avg - x;
}

Vector v_value(vector<Vertex> &star, Vector x2, Vector x1){
	for(unsigned int i = 0; i < star.size();++i){
		Vector t = star[i]->point() - CGAL::ORIGIN;
		if( t == x1 ){
			int next_i = ( i + 1 )%star.size();
			size_t prev_i = ( i + star.size() - 1 )%star.size();
			//cout<< star.size()<<" :size "<<prev_i<<" * "<<next_i<<endl;
			Vector e_prev = ( star[prev_i]->point() - CGAL::ORIGIN ) - x2;
			Vector e_next = ( star[next_i]->point() - CGAL::ORIGIN ) - x2;
			//cout<<e_prev<<" ** "<<e_next<<endl;
			return e_prev - e_next;
		}
	}
	return Vector(0,0,0);
}

int unrepeated_label(Vertex v1 , Vertex v2 ){
	for(unsigned int i = 0 ; v1->info().m_regions_around.size(); ++i){
		int r = v1->info().m_regions_around[i];
		if( r == v2->info().m_regions_around[0] || r == v2->info().m_regions_around[1] ) continue;
		return r;
	}
	return -1;
}

bool compare_labels_triple_vertices(Vertex v1 , Vertex v2 ){
	std::sort(v1->info().m_regions_around.begin(), v1->info().m_regions_around.end());
    std::sort(v2->info().m_regions_around.begin(), v2->info().m_regions_around.end());
    return v1->info().m_regions_around == v2->info().m_regions_around;
}

void double_node_relaxation(vector<Vertex> &star, Vertex v){
	Vector x = v->point() - CGAL::ORIGIN;
	Vector dxs = middle_point(star) - x;
	Vector n = average_normal(star,x);
	Vector newV = x + dxs - (dxs*n)*n;
	v->set_point(CGAL::ORIGIN + newV);
}

void triple_node_relaxation(Delaunay &dt, Vertex v){
	Vector V,n1,n2,t,dx,newV;
	double normal;
	V = v->point() - CGAL::ORIGIN;
	vector<Vertex> star1 = v->info().m_s_stars[0];
	vector<Vertex> star2 = v->info().m_s_stars[1];
	n1 = normal_edges(star1,V);
	n2 = normal_edges(star2,V);
	normal = norm( CGAL::cross_product( n1 , n2 ) );
	if( normal > 0 ){
		t = CGAL::cross_product( n1 , n2 )/normal;
		dx = star_vertex_average(dt,v) - V;
		//cout<<"dx: "<<dx<<endl;
		newV = V + dx;
		v->set_point(CGAL::ORIGIN + newV);
		//cout<<"newV: "<<newV<<endl;
	}
}

void double_edge_relaxation(Delaunay &dt, Vertex v1 , Vertex v2){
	Vector A,A1,A2,V1,V2,v,avg1,avg2,xs1,xs2,dx1,dx2,n,newV1,newV2;
	size_t n1,n2;
	double normA,h;
	int label = common_label2(v1,v2);
		
	bool draw = 0;
	vector<Point> test;
	test.push_back(v1->point());
	test.push_back(v2->point());

	if(draw){
		drawMesh(dt,test);
		getchar();
	}

	vector<Vertex> star1 = v1->info().m_s_stars[index_label_in_star(v1,label)];
	vector<Vertex> star2 = v2->info().m_s_stars[index_label_in_star(v2,label)];

	V1 = v1->point() - CGAL::ORIGIN;
	V2 = v2->point() - CGAL::ORIGIN;		

	A1 = normal_edges(star1,V1);
	A2 = normal_edges(star2,V2);

	v = v_value( star2 , V2 , V1 );
	avg1 = average_without_point(star1,V2);
	avg2 = average_without_point(star2,V1);

	n1 = star1.size();
	n2 = star2.size();

	xs1 = ( avg2 + n2*avg1 )/( n1*n2 - 1 );
	xs2 = ( xs1 + avg2 )/n2;

	dx1 = xs1 - V1;
	dx2 = xs2 - V2;

	A = A1 + A2 + CGAL::cross_product( v , dx1 - dx2 );
		
	normA = norm(A);
	if( normA > 0){
		n = A/normA;
		h = -( dx1*A1 + dx2*A2 + dx2*CGAL::cross_product( v , dx1 ) )/normA;
		newV1 = V1 + dx1 + h*n;
		newV2 = V2 + dx2 + h*n;
		v1->set_point(CGAL::ORIGIN + newV1);
		v2->set_point(CGAL::ORIGIN + newV2);

		if(draw){
			cout<<"V1: "<<V1<<" V2: "<<V2<<endl;
			cout<<"newV1: "<<newV1<<" newV2: "<<newV2<<endl;
			test.clear();
			test.push_back(CGAL::ORIGIN + newV1);
			test.push_back(CGAL::ORIGIN + newV2);
			drawMesh(dt,test);
			getchar();			
		}
	}
}

void triple_edge_relaxation(Delaunay &dt, Vertex v1 , Vertex v2){
	Vector A_r1,A_r2,A1_r1,A2_r1,A1_r2,A2_r2,V1,V2,v_r1,v_r2,dx1,dx2,newV1,newV2,avg1,avg2,xs1,xs2;
	double	det,norm_A_r1,norm_A_r2,g_r1,g_r2,h_r1,h_r2;
	int n1,n2;

	int r1 = v1->info().m_regions_around[0];
	int r2 = v1->info().m_regions_around[1];

	vector<Vertex> star1_r1 = v1->info().m_s_stars[ 0 ];
	vector<Vertex> star2_r1 = v2->info().m_s_stars[index_label_in_star(v2,r1)];

	vector<Vertex> star1_r2 = v1->info().m_s_stars[ 1 ];
	vector<Vertex> star2_r2 = v2->info().m_s_stars[index_label_in_star(v2,r2)];

	V1 = v1->point() - CGAL::ORIGIN;
	V2 = v2->point() - CGAL::ORIGIN;		

	bool draw = 1;
	vector<Point> test;
	test.push_back(v1->point());
	test.push_back(v2->point());

	if(draw){
		drawMesh(dt,test);
		getchar();
	}

	A1_r1 = normal_edges(star1_r1,V1);
	A2_r1 = normal_edges(star2_r1,V2);
	v_r1 = v_value(star2_r1, V2 , V1 );

	A1_r2 = normal_edges(star1_r2,V1);
	A2_r2 = normal_edges(star2_r2,V2);
	v_r2 = v_value(star2_r2, V2 , V1 );

	avg1 = average_around_vertex(dt,v1,V2,n1);
	avg2 = average_around_vertex(dt,v2,V1,n2);

	xs1 = ( avg2 + n2*avg1 )/( n1*n2 - 1 );
	xs2 = ( xs1 + avg2 )/n2;

	dx1 = xs1 - V1;
	dx2 = xs2 - V2;

	A_r1 = A1_r1 + A2_r1 + CGAL::cross_product( v_r1 , dx1 - dx2 );
	A_r2 = A1_r2 + A2_r2 + CGAL::cross_product( v_r2 , dx1 - dx2 );

	norm_A_r1 = norm(A_r1);
	norm_A_r2 = norm(A_r2);

	det = norm_A_r1*norm_A_r1*norm_A_r2*norm_A_r2 - pow( A_r1*A_r2 , 2 );
	if(det > 0 && CGAL::cross_product(A_r1,A_r2) != Vector(0,0,0) ){

		g_r1 = - dx1*A1_r1 - dx2*A2_r1 - dx2*CGAL::cross_product( v_r1 , dx1 );
		g_r2 = - dx1*A1_r2 - dx2*A2_r2 - dx2*CGAL::cross_product( v_r2 , dx1 );
		h_r1 = ( norm_A_r2*norm_A_r2 *g_r1 - A_r1*A_r2*g_r2 )/det;
		h_r2 = ( - A_r1*A_r2*g_r1 + norm_A_r1*norm_A_r1*g_r2 )/det;
		newV1 = V1 + dx1 + h_r1*A_r1 + h_r2*A_r2;
		newV2 = V2 + dx2 + h_r1*A_r1 + h_r2*A_r2;
		v1->set_point(CGAL::ORIGIN + newV1);
		v2->set_point(CGAL::ORIGIN + newV2);

		Vector Dx1,Dx2;
		Dx1 = newV1 - V1;
		Dx2 = newV2 - V2;

		double vol1 = Dx1*A1_r1 + Dx2*A2_r1 + Dx2*CGAL::cross_product(v_r1,Dx1);
		double vol2 = Dx1*A1_r2 + Dx2*A2_r2 + Dx2*CGAL::cross_product(v_r2,Dx1);

		if(draw){
			cout<<"total volume : "<<vol1 + vol2<<endl;
			cout<<"A_r1: "<<A_r1<<" A_r2: "<<A_r2<<endl;
			cout<<"h_r1: "<<h_r1<<" h_r2: "<<h_r2<<endl;
			cout<<"dx1: "<<dx1<<" dx2: "<<dx2<<endl;
			cout<<"V1: "<<V1<<" V2: "<<V2<<endl;
			cout<<"newV1: "<<newV1<<" newV2: "<<newV2<<endl;
			test.clear();
			test.push_back(CGAL::ORIGIN + newV1);
			test.push_back(CGAL::ORIGIN + newV2);
			drawMesh(dt,test);
			getchar();			
		}
	}
}

void mixed_edge_relaxation(Delaunay &dt, Vertex v1 , Vertex v2,int &count){
	if( !( v1->info().m_regions_around.size() == 3 && v2->info().m_regions_around.size() == 2 ) ) return;

	Vector A_r1,A_r2,A1_r1,A2_r1,A1_r2,A2_r2,V1,V2,v_r1,v_r2,dx1,dx2,newV1,newV2,avg1,avg2,xs1,xs2;
	double	det,norm_A_r1,norm_A_r2,g_r1,g_r2,h_r1,h_r2;
	int n1,n2;

	//v1 is the triple node and v2 is doble node
	int r1 = common_label2(v1,v2);
	int r2 = unrepeated_label(v1,v2);

	vector<Vertex> star1_r1 = v1->info().m_s_stars[index_label_in_star(v1,r1)];
	vector<Vertex> star2_r1 = v2->info().m_s_stars[index_label_in_star(v2,r1)];

	vector<Vertex> star1_r2 = v1->info().m_s_stars[index_label_in_star(v1,r2)];
	//vector<Vertex> star2_r2 = v2->info().m_s_stars[index_label_in_star(v2,r2)];

	V1 = v1->point() - CGAL::ORIGIN;
	V2 = v2->point() - CGAL::ORIGIN;

	bool draw = 0;
	
	vector<Point> test;
	test.push_back(v1->point());
	test.push_back(v2->point());

	/*if(draw && count > 33 ){
		drawMesh(dt,test);
		getchar();
	}*/

	A1_r1 = normal_edges(star1_r1,V1);
	A2_r1 = normal_edges(star2_r1,V2);
	v_r1 = v_value(star2_r1, V2, V1 );

	A1_r2 = normal_edges(star1_r2,V1);
	A2_r2 = Vector(0,0,0);
	v_r2 = Vector(0,0,0);

	avg1 = average_around_vertex(dt,v1,V2,n1);
	//avg2 = average_around_vertex(dt,v2,V1,n2);

	n2 = (int)star2_r1.size();
	avg2 = average_without_point(star2_r1,V1);

	xs1 = ( avg2 + n2*avg1 )/( n1*n2 - 1 );
	xs2 = ( xs1 + avg2 )/n2;

	dx1 = xs1 - V1;
	dx2 = xs2 - V2;

	//check this part

	A_r1 = A1_r1 + A2_r1 + CGAL::cross_product( v_r1 , dx1 - dx2 );
	A_r2 = A1_r2 + A2_r2 + CGAL::cross_product( v_r2 , dx1 - dx2 );

	norm_A_r1 = norm(A_r1);
	norm_A_r2 = norm(A_r2);

	double normCross = norm( CGAL::cross_product(A_r1,A_r2) );

	det = norm_A_r1*norm_A_r1*norm_A_r2*norm_A_r2 - pow( A_r1*A_r2 , 2 );
	if( det > 0 ){
		++count;
		g_r1 = - dx1*A1_r1 - dx2*A2_r1 - dx2*CGAL::cross_product( v_r1 , dx1);
		g_r2 = - dx1*A1_r2 - dx2*A2_r2 - dx2*CGAL::cross_product( v_r2 , dx1);
		h_r1 = ( norm_A_r2*norm_A_r2 *g_r1 - A_r1*A_r2*g_r2 )/det;
		h_r2 = ( - A_r1*A_r2*g_r1 + norm_A_r1*norm_A_r1*g_r2 )/det;
		newV1 = V1 + dx1 + h_r1*A_r1 + h_r2*A_r2;
		newV2 = V2 + dx2 + h_r1*A_r1 + h_r2*A_r2;

		/*if( norm(V1 - newV1) > 5) return;
		if( norm(V2 - newV2) > 5) return;*/

		//if(draw && count > 29){
		if(draw){
			/*for(unsigned int i = 0; i < v1->info().m_regions_around.size(); ++i)
				cout<<v1->info().m_regions_around[i]<<" ";
			cout<<endl;
			for(unsigned int i = 0; i < v2->info().m_regions_around.size(); ++i)
				cout<<v2->info().m_regions_around[i]<<" ";
			cout<<endl;*/
			/*cout<<"cross produt: "<<CGAL::cross_product(A_r1,A_r2)<<endl;
			cout<<"norm : "<<normCross<<endl;*/
			cout<<"label r1: "<<r1<<"  label r2: "<<r2<<endl;	

			Vector u1 = unit_vector(A1_r1);
			Vector u2 = unit_vector(A2_r1);
			Vector u3 = unit_vector(A1_r2);

			/*test.push_back(CGAL::ORIGIN + V1+ 20*u1);
			test.push_back(CGAL::ORIGIN + V2+ 20*u2);
			test.push_back(CGAL::ORIGIN + V1+ 20*u3);*/

			test.push_back(CGAL::ORIGIN + xs1);
			test.push_back(CGAL::ORIGIN + xs2);

			drawMesh(dt,test);
			getchar();

			/*vertexRegion(dt,v1,test,1);
			getchar();*/
		}

		v1->set_point(CGAL::ORIGIN + newV1);
		v2->set_point(CGAL::ORIGIN + newV2);

		Vector Dx1,Dx2;
		Dx1 = newV1 - V1;
		Dx2 = newV2 - V2;

		double vol1 = Dx1*A1_r1 + Dx2*A2_r1 + Dx2*CGAL::cross_product(v_r1,Dx1);
		double vol2 = Dx1*A1_r2 + Dx2*A2_r2 + Dx2*CGAL::cross_product(v_r2,Dx1);

		//if(draw && count > 29 ){
		if(draw){
			cout<<" count "<<count<<endl;

			/*for(unsigned int i = 0; i < v1->info().m_regions_around.size(); ++i)
				cout<<v1->info().m_regions_around[i]<<" ";
			cout<<endl;
			for(unsigned int i = 0; i < v2->info().m_regions_around.size(); ++i)
				cout<<v2->info().m_regions_around[i]<<" ";
			cout<<endl;*/

			cout<<"vol1 : "<<vol1 <<endl;
			cout<<"vol2 : "<<vol2 <<endl;
			cout<<"total volume : "<<vol1 + vol2<<endl;
			//cout<<"label r1: "<<r1<<"  label r2: "<<r2<<endl;	

			cout<<"avg1: "<<avg1<<" avg2: "<<avg2<<endl;
			cout<<"A1_r1: "<<A1_r1<<" A2_r1: "<<A2_r1<<endl;
			cout<<"A1_r2: "<<A1_r2<<" A2_r2: "<<A2_r2<<endl;
			cout<<"v_r1: "<<v_r1<<" v_r2: "<<v_r2<<endl;
			cout<<"A_r1: "<<A_r1<<" A_r2: "<<A_r2<<endl;
			cout<<"g_r1: "<<g_r1<<" g_r2: "<<g_r2<<endl;			
			cout<<"h_r1: "<<h_r1<<" h_r2: "<<h_r2<<endl;
			cout<<"dx1: "<<dx1<<" dx2: "<<dx2<<endl;
			cout<<"V1: "<<V1<<" V2: "<<V2<<endl;
			cout<<"xs1: "<<xs1<<" xs2: "<<xs2<<endl;
			cout<<"newV1: "<<newV1<<" newV2: "<<newV2<<endl;
			//test.clear();

			

			test.push_back(CGAL::ORIGIN + newV1);
			test.push_back(CGAL::ORIGIN + newV2);
			drawMesh(dt,test);
			getchar();			
		}
	}
}

void vertex_relaxation( Delaunay &dt, vector<Vertex> &boundaryVertices ){
	for ( unsigned int i = 0 ; i < boundaryVertices.size() ; ++i){
		Vertex v = boundaryVertices[i] ;
		if( v->info().is_corner() || v->info().m_regions_around.size() > 3  || v->info().is_singular_2() )		continue;
		if( v->info().m_regions_around.size() == 2 ){
			vector<Vertex> star = v->info().m_s_stars[0];
			double_node_relaxation( star, v );
		}
		else{
			//cout<<"size: "<<v->info().m_regions_around.size()<<endl;
			triple_node_relaxation(dt,v);
		}
	}
}

void edge_relaxation( Delaunay &dt, vector<Edge> &edges){
	Vertex v1,v2;
	int count = 0;
	for ( unsigned int i = 0; i < edges.size(); ++i){
		v1 = edges[i].first->vertex(edges[i].second);
		v2 = edges[i].first->vertex(edges[i].third);
		if( ( v1->info().is_corner() || v1->info().m_regions_around.size() > 3 || v1->info().is_singular_2() ) || 
			( v2->info().is_corner() || v2->info().m_regions_around.size() > 3 || v1->info().is_singular_2() ) )		continue;

		if( v1->info().m_regions_around.size() == 2 && v2->info().m_regions_around.size() == 2){
			double_edge_relaxation(dt,v1,v2);
			continue;
		}
		if( v1->info().m_regions_around.size() == 3 && v2->info().m_regions_around.size() == 3){
			if( compare_labels_triple_vertices(v1,v2) )
				triple_edge_relaxation(dt,v1,v2);
			continue;
		}
		/*if( v1->info().m_regions_around.size() == 3 )
			mixed_edge_relaxation(dt,v1,v2,count);
		else
			mixed_edge_relaxation(dt,v2,v1,count);*/
	}
}

void kuprat(Delaunay &dt, int iteration){
	mark_boundary(dt);
	vector<Vertex> boundaryVertices;
	vector<Edge> edges;
	boundary_vertices(dt,boundaryVertices);

	edges_of_boundary(dt,edges);
	initiate_process(dt,boundaryVertices);

	for( int i = 0; i < iteration; ++i){
		vertex_relaxation(dt,boundaryVertices);
		edge_relaxation(dt,edges);
	}
	restart_vertex_state(dt);
	cout<<"end kuprat"<<endl;


	/*drawMesh(dt);
	getchar();*/
}