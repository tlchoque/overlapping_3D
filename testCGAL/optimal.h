//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues> 
//
//using Eigen::MatrixXd;
//using Eigen::Matrix3d;
//using Eigen::EigenSolver;

void setOrientationOptimal(int &idx1, int &idx2,int & idx3, int flag ){
	vector<int> list(3);
	list[0]=idx1;
	list[1]=idx2;
	list[2]=idx3;
	sort(list.begin(),list.end());
	idx1 = list[0];
	if(flag == 0 || flag==2){
		idx2 = list[2];
		idx3 = list[1];
	}
	else{
		idx2 = list[1];
		idx3 = list[2];
	}
}

double normSquare(Vector v){//norm-1
	return v*v;
	//return  abs(v.x()) + abs(v.y()) + abs(v.z());
}

Vector unit_normal_facet(Delaunay &dt,Facet &f){
	int idx1 = (f.second+1)%4;
	int idx2 = (f.second+2)%4;
	int idx3 = (f.second+3)%4;
	setOrientationOptimal(idx1,idx2,idx3,f.second);
	return unit_normal_triangle( f.first->vertex(idx1)->point(), f.first->vertex(idx2)->point(), f.first->vertex(idx3)->point());
}

double areaFacet(Delaunay &dt,Facet &f){
	int idx1 = (f.second+1)%4;
	int idx2 = (f.second+2)%4;
	int idx3 = (f.second+3)%4;
	setOrientationOptimal(idx1,idx2,idx3,f.second);
	return triangle_area(f.first->vertex(idx1)->point(), f.first->vertex(idx2)->point(), f.first->vertex(idx3)->point());
}

double volumeAroundVertex(Delaunay &dt,Vertex &v){
	list<Cell> incidentCells;
	dt.finite_incident_cells(v,back_inserter(incidentCells));
	double sum=0;
	for(std::list<Cell>::iterator it = incidentCells.begin(); it != incidentCells.end(); it++)
		sum+=abs( dt.tetrahedron(*it).volume() ) ;
	return sum;
}

double sumNormSquare(Delaunay &dt,Facet &f){
	Vertex v0 = f.first->vertex( (f.second) );
	double sum=0;
	for(int i = 1; i < 4 ; ++i)
		//sum+=normSquare(rest( f.first->vertex( (f.second+i)%4 )->point() , v0->point() ));
		sum+=normSquare( f.first->vertex( (f.second+i)%4 )->point() - v0->point() );
	return sum;
}

Point smooth(Delaunay &dt, Vertex &v){
	list<Cell> incidentCells;
	dt.finite_incident_cells(v,back_inserter(incidentCells));
	Vector total=Vector(0,0,0);
	for(std::list<Cell>::iterator it = incidentCells.begin(); it != incidentCells.end(); it++){ 
		Cell c = *it;
		Facet f = Facet(c,c->index(v));
		double sumNS= sumNormSquare(dt,f);
		double S = areaFacet(dt,f);
		Vector n = unit_normal_facet(dt,f);
		/*Point cenc = CGAL::centroid(dt.triangle(f));
		cenc = sum(cenc,scalar(n,50));
		verticesToDraw.push_back(cenc);
		test.push_back(f);

		drawMesh(dt,test,verticesToDraw);
		cout<<"are: "<<S<<" sum: "<<sumNS<<endl;
		getchar();*/
		n = Vector( sumNS * S *n.x(), sumNS * S *n.y() , sumNS * S *n.z() );
		total = total + n;
	}
	total = Vector(total.x()/3,total.y()/3,total.z()/3);
	double ohmio = volumeAroundVertex(dt,v);
	Point x = Point( v->point().x() - total.x() / (2*ohmio) , v->point().y() - total.y() / (2*ohmio), v->point().z() - total.z() / (2*ohmio) );
	return x;
}

// algorithm for boundary smoothing 

void getOrientation(int &idx1, int &idx2,int & idx3, int flag , int vindex){// 
	vector<int> list(3);
	list[0]=idx1;
	list[1]=idx2;
	list[2]=idx3;
	sort(list.begin(),list.end());
	idx1 = list[0];
	if(flag == 0 || flag==2){
		int cen = list[2];
		list[2] = list[1];
		list[1] = cen;
	}
	int index=0;
	for(int i = 0; i < 3; ++i){
		if(list[i] == vindex)
			index = i;
	}
	// order the orientation
	if(index == 0){// the vindex is located at 1st
		idx1 = list[1];
		idx2 = list[2];
	}
	else if(index == 1){// the vindex is located at 2nd
		idx1 = list[2];
		idx2 = list[0];
	}
	else{// the vindex is located at 3er
		idx1 = list[0];
		idx2 = list[1];
	}
	idx3=vindex;
}

vector<Vertex> getOrderedVertices(Delaunay &dt, Vertex &v,int label){
	map<Vertex,Vertex> vv;
	deque<void*> cpyStar = v->info().m_star[label];
	for(unsigned int i=0;i < cpyStar.size();++i){
		Facet f =*(Facet *)cpyStar[i];
		Cell c = f.first;
		if(c->info().m_label==label){
			f = dt.mirror_facet(f);
			c = f.first;
		}
		int idx1 = (f.second+1)%4;
		int idx2 = (f.second+2)%4;
		int idx3 = (f.second+3)%4;

		getOrientation(idx1,idx2,idx3,f.second,c->index(v));// f.second != v
		vv[c->vertex(idx1)] = c->vertex(idx2);
		//cout<<c->vertex(idx1)->point()<<" // "<<c->vertex(idx2)->point()<<endl;
		if(i == 0)
			vv[v] = c->vertex(idx1);
	}
	vector<Vertex> orderedVertices;// the first is repeated at last
	Vertex it = v;

	for(unsigned int i=0;i < cpyStar.size();++i){
		it = vv[it];
		orderedVertices.push_back(it);
	}
	return orderedVertices;
}

//Vector norm_tangent_plane(Vertex &v,vector<Vertex> &ov){
//	Vector sum = Vector(0,0,0);
//	for(unsigned int i = 0; i < ov.size() - 1; ++i){
//		double S = triangle_area( v->point(),ov[i]->point(),ov[i+1]->point() );
//		Vector n = unit_normal_triangle( v->point(),ov[i]->point(),ov[i+1]->point() );
//		sum = Vector(sum.x() + S*n.x(), sum.y() + S*n.y(), sum.z() + S*n.z());
//	}
//	sum = unit_vector(sum);
//	return sum;
//}

void getOrthogonals(Delaunay &dt,Vector &n, Vector &a,Vector &b){//n->norm a, b orthogonal vectors in plane
	//first orhogonal vector
	a = Vector(0,1,-n.y()/n.z());
	a = unit_vector(a);
	Vector aux = Vector(1,0, -n.x()/n.z());
	//1st projection
	double p1 = aux*n;
	Vector proj1 = Vector( p1*n.x(), p1*n.y(),p1*n.z() );
	//2nd projection
	double p2 = aux*a;
	Vector proj2 = Vector( p2*a.x(), p2*a.y(),p2*a.z() );
	//last orhogonal vector
	b = Vector( aux.x() - proj1.x() - proj2.x(), aux.y() - proj1.y() - proj2.y() , aux.z() - proj1.z() - proj2.z());
	b = unit_vector(b);
}

double volumeAroundVertexBoundary(Delaunay &dt,Vertex &v,int label){
	list<Cell> incidentCells;
	dt.finite_incident_cells(v,back_inserter(incidentCells));
	double sum=0;
	for(std::list<Cell>::iterator it = incidentCells.begin(); it != incidentCells.end(); it++){
		if( (*it)->info().m_label == label )
			sum+=abs( dt.tetrahedron(*it).volume() ) ;
	}
	return sum;
}

Point smoothB1(Delaunay &dt, Vertex &v, vector<Vertex> &vertices,int label){
	Vector n = norm_tangent_plane(vertices,v);	
	Vector s,t;
	getOrthogonals(dt,n,s,t);
	double E=0,F=0,G=0,H=0,I=0;
	cout<<"1 ort vector: "<<s<<endl;
	cout<<"2 ort vector: "<<t<<endl;


	/*Point Y1 = rest( vertices[0]->point() ,v->point() );
	Point Y2 = rest( vertices[1]->point() ,v->point() );
	Point crossY = cross(Y1,Y2);*/

	for(unsigned int i = 0; i < vertices.size(); ++i){
		int next_i = (i+1)%vertices.size();

		Vector Yi = vertices[i]->point() - v->point() ;
		Vector Yi1 = vertices[next_i]->point() - v->point() ;
		Vector sumYi = Yi + Yi1;
		Vector crossYi = CGAL::cross_product(Yi,Yi1);
		double ssumYi = s * sumYi ;
		double tsumYi = t * sumYi ;
		double scrossYi = s * crossYi ;
		double tcrossYi = t * crossYi ;
		double Yisquare = Yi * Yi;
		double Yi1square =  Yi1 * Yi1;
		double YiYi1 = Yi * Yi1;
		double term = Yisquare + Yi1square + YiYi1;	

		E += ssumYi * scrossYi;
		F += tsumYi * tcrossYi; 
		G += ssumYi * tcrossYi + tsumYi * scrossYi;
		H += term * scrossYi;
		I += term * tcrossYi;	

		/*double tcrossY = inner( t , crossY );
		F += tsumYi * tcrossY; */
	}
	
	vector<Facet> test;// for drawing
	vector<Point> verticesToDraw;//for drawing

	list<Cell> incidentCells;
	dt.finite_incident_cells(v,back_inserter(incidentCells));
	Vector total = Vector(0,0,0);
	for(std::list<Cell>::iterator it = incidentCells.begin(); it != incidentCells.end(); it++){
		Cell c = *it;
		if( c->info().m_label == label){
			Facet f = Facet(c,c->index(v));
			double sumNS = sumNormSquare(dt,f);
			double S = areaFacet(dt,f);
			Vector n = unit_normal_facet(dt,f);

			/*Point cenc = CGAL::centroid(dt.triangle(f));
			cenc = sum(cenc,scalar(n,50));
			verticesToDraw.push_back(cenc);*/

			//cout<<"sum: "<<sumNS<<" area: "<<S<<" unit: "<<n<<endl;
			n = Vector( sumNS * S *n.x() , sumNS * S *n.y() , sumNS * S *n.z() );
			total = total + n;
		}
	}
	double ohmio = volumeAroundVertexBoundary(dt,v,label);

	E =	ohmio/4 - E / 60;
	F =	ohmio/4 - F / 60;
	G = - G / 60;
	H = s*total/ 12 - H / 60;
	I = t*total / 12 - I / 60;

	double det = 4 * E  * F - G * G ;
	double udet = - 2 * H * F + G * I ;
	double vdet = - 2 * E * I +  H * G ;
	double U = udet/det;
	double V = vdet/det;
	Vector x = (v->point() - CGAL::ORIGIN) +  U*s + V*t ;

	Point X = Point(x.x(), x.y(),x.z() );

	cout<<"volume: "<<ohmio<<endl;
	cout<<E<<" "<<F<<" "<<G<<" "<<H<<" "<<I<<endl;
	cout<<"U: "<<U<<" V: "<<V<<endl;
	//---------debug by drawing

	Point d = Point( 50*n.x() + v->point().x() , 50*n.y() + v->point().y(), 50*n.z() + v->point().z() );
	Point e = Point( 50*s.x() + v->point().x() , 50*s.y() + v->point().y(), 50*s.z() + v->point().z() );
	Point f = Point( 50*t.x() + v->point().x() , 50*t.y() + v->point().y(), 50*t.z() + v->point().z() );

	verticesToDraw.push_back(v->point());
	verticesToDraw.push_back(d);	
	verticesToDraw.push_back(e);	
	verticesToDraw.push_back(f);
	verticesToDraw.push_back(X);

	vertexRegion(dt,v,verticesToDraw,label);	
	//drawMesh(dt,test,verticesToDraw);
	getchar();

	//---------	

	return X;
}

Point smoothB2(Delaunay &dt, Vertex &v, Vector &d,vector<Vertex> &vertices,int label){
	double A=0,B=0;
	for(unsigned int i = 0; i < vertices.size() ; ++i){
		int next_i = (i+1)%vertices.size();
		Vector Yi =  vertices[i]->point() - v->point() ;
		Vector Yi1 = vertices[next_i]->point() - v->point() ;
		Vector sumYi = Yi + Yi1;
		Vector crossYi = CGAL::cross_product(Yi,Yi1);
		double dsumYi = d * sumYi ;
		double dcrossYi = d * crossYi;

		double Yisquare = Yi * Yi;
		double Yi1square = Yi1 * Yi1;
		double YiYi1 = Yi * Yi1 ;

		A += dsumYi * dcrossYi;
		B += (Yisquare + Yi1square + YiYi1) * dcrossYi;
	}

	list<Cell> incidentCells;
	dt.finite_incident_cells(v,back_inserter(incidentCells));
	Vector total = Vector(0,0,0);
	for(std::list<Cell>::iterator it = incidentCells.begin(); it != incidentCells.end(); it++){
		Cell c = *it;
		if( c->info().m_label == label){
			Facet f = Facet(c,c->index(v));
			double sumNS= sumNormSquare(dt,f);
			double S = areaFacet(dt,f);
			Vector n = unit_normal_facet(dt,f);
			//cout<<"sum: "<<sumNS<<" area: "<<S<<" unit: "<<n<<endl;
			n = Vector( sumNS * S *n.x() / 3, sumNS * S *n.y() / 3 , sumNS * S *n.z() / 3 );
			total = total + n;
		}
	}
	double ohmio = volumeAroundVertexBoundary(dt,v,label);
	A = ohmio/4 - A/60;
	B = d*total/12 - B/60;
	double f = - B /(2 * A);
	//Point x = sum(v->point(), scalar( d, f ) );
	Vector x = (v->point() - CGAL::ORIGIN) + f*d ;
	return Point(x.x(),x.y(),x.z());
}

double fewDecimals(double &x){
	double dec= 10;
	double y = floor(x * dec)/dec;
	return y;
}

bool pairCompare(const std::pair<double, Vector>& firstElem, const std::pair<double, Vector>& secondElem) {
	return firstElem.first > secondElem.first;
}

//Point smoothBoundary(Delaunay &dt, Vertex &v){
//	int label=0;
//	for(unsigned int i=1;i < v->info().m_star.size();++i){
//		//cout<<"size star i: "<<v->info().m_star[i].size()<<endl;
//		if( v->info().m_star[i].size() > v->info().m_star[label].size() )
//			label = i;
//	}
//	//cout<<"label: "<<label<<" size star: "<<v->info().m_star[label].size()<<endl;
//	vector<Vertex> vertices = getOrderedVertices(dt,v,label);
//	cout<<"vertices size: "<<vertices.size()<<endl;
//	Matrix3d M,aux;
//	M << 0,0,0,
//		0,0,0,
//		0,0,0;
//
//	for(unsigned int i = 0; i < vertices.size() ;++i){
//		int next_i = (i+1)%vertices.size();
//		double S = triangle_area( v->point(),vertices[i]->point(),vertices[next_i]->point() );
//		Vector n = unit_normal_triangle( v->point(),vertices[i]->point(),vertices[next_i]->point() );
//		double xy = S * n.x() * n.y();
//		double xz = S * n.x() * n.z();
//		double yz = S * n.y() * n.z();
//		//cout<<"unit normal: "<<n<<endl;
//		aux<< S*n.x()*n.x(), xy				, xz,
//			  xy		   , S*n.y()*n.y()  , yz,
//			  xz		   , yz				, S*n.z()*n.z();
//		M = M + aux;
//	}
//	//cout<<M<<endl;
//	EigenSolver<Matrix3d> es(M);
//	vector<pair<double,Vector>> eigens;
//	for(unsigned int i = 0; i < 3 ; ++i){
//		eigens.push_back( make_pair( es.eigenvalues()[i].real(), Vector(es.eigenvectors().col(i).x().real() , es.eigenvectors().col(i).y().real() , es.eigenvectors().col(i).z().real() ) ) );
//	}
//	std::sort(eigens.begin(), eigens.end(), pairCompare);
//
//	//----------depurate by  drawing
//
//	//vector<Point> verticesToDraw;//for drawing
//	//vector<Facet> test;// for drawing
//
//	//verticesToDraw.push_back(v->point());
//
//	//----------
//
//	for(unsigned int i = 0; i < 3 ; ++i){
//		cout<<eigens[i].first<<" , ";
//		eigens[i] = make_pair( fewDecimals(eigens[i].first) , eigens[i].second );
//		//verticesToDraw.push_back(sum(v->point(),scalar( Point(eigens[i].second.x(), eigens[i].second.y(), eigens[i].second.z()) , eigens[i].first/1000 ) ) );
//	}
//	cout<<endl;
//
//	/*vertexRegion(dt,v,verticesToDraw);
//	getchar();*/
//
//
//	Point x = v->point();
//	if( eigens[0].first > eigens[1].first ){
//		x = smoothB1(dt,v,vertices,label);
//	}
//	else if( eigens[1].first > eigens[2].first ){
//		x = smoothB2(dt,v,eigens[2].second,vertices,label);
//	}
//
//	return x;
//}

void preMesh(Delaunay &dt){
	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++)
		it->info().m_label = 1;
}

//void optimalDelaunay(Delaunay &dt){
//	preMesh(dt);
//	vector<Facet> test;// for drawing
//	vector<Point> verticesToDraw;//for drawing
//
//	mark_boundary(dt);
//	vector<Vertex> boundaryVertices;
//	boundary_vertices(dt,boundaryVertices);
//	for(int i = 0; i < 1; ++i){
//		for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) {
//			Vertex v = vi;
//			if(v->info().m_boundary){
//				Point vs = smoothBoundary(dt,v);
//				cout<<"yes boundary: "<<v->point()<<" ** "<<vs<<endl;
//				if(vs != v->point() ){
//					/*map<Point,int> mpLabel;
//					beforeMove(dt,v,mpLabel,vs);*/
//					dt.move(v,vs);	
//					//afterMove(dt,v,mpLabel);
//					getchar();
//				}
//			}
//			else{
//				Point vs = smooth(dt,v);
//				cout<<"no  boundary: "<<v->point()<<" -- "<<vs<<endl;
//				map<Point,int> mpLabel;
//				dt.move(v,vs);	
//				getchar();
//			}
//		} 
//	}
//}