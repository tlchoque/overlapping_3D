
void generation(Delaunay & dt, int  numPoints, int windowWidth, int windowHeight,int windowDepth ){
	vector<Point> rPoints;
	int max = 10;
	Vector *p = new Vector[max]; 
	Vector *p_Temp = NULL;   
	int nv = 0;
	double x, y, z;
	bool ok = false;
	randomize();
	nv = 0;
	p = new Vector[max];
	while (nv != numPoints){
		do{
			ok = true;
			x = (double)(rand()%windowWidth);
			y = (double)(rand()%windowHeight);
			z = (double)(rand()%windowDepth);
			for(int n_Cpt = 0; n_Cpt <= nv; ++n_Cpt){
				if( (x == p[n_Cpt].x() ) && (y == p[n_Cpt].y() ) && (z == p[n_Cpt].z() ) ) {ok = false;break;}
			}
		}while(!ok);
		if (nv >= max){
			max = max * 2;            
			p_Temp = new Vector[max]; 
			for (int i = 0; i < nv; ++i) {
				p_Temp[i] = p[i];  
			}
			delete []p;  
			p = p_Temp; 
		}   
		p[nv] = Vector(x,y,z); 
		nv++;
	}
	for (int i = 0; i < nv; ++i) 
		rPoints.push_back(Point(p[i].x() ,p[i].y() ,p[i].z() ));   
	delete []p_Temp; 

	dt.insert(rPoints.begin(),rPoints.end());
}

void load_points2(std::vector<Point>& rPoints){
	rPoints.push_back(Point(25,21,30)); 
	//rPoints.push_back(Point(30,30,30)); 
	rPoints.push_back(Point(15,55.98,30));  
	rPoints.push_back(Point(45,55.98,30)); 
	rPoints.push_back(Point(0,30,30));   
	rPoints.push_back(Point(60,30,30)); 
	rPoints.push_back(Point(15,4.019,30));   
	rPoints.push_back(Point(45,4.019,30));   

	rPoints.push_back(Point(30,4.019,30)); 

	rPoints.push_back(Point(30,30,0));   
	rPoints.push_back(Point(30,30,60)); 


	/*Point a,b,c;
	a= Point(30,30,60);
	b= Point(0,30,30);
	c = Point(15,4.019,30);
	Point cen = CGAL::centroid(a,b,c);
	rPoints.push_back(cen); 

	Point q = Point(30,30,30);*/

	//computed recently
	//rPoints.push_back(Point(26.25,23.504,52.5)); 
	//rPoints.push_back(Point(22.5,17.0095,45)); 
 }

void load_points(std::vector<Point>& rPoints){
	rPoints.push_back(Point(24,23,30)); 
	rPoints.push_back(Point(15,55.98,30));  
	rPoints.push_back(Point(45,55.98,30)); 
	//rPoints.push_back(Point(30,55.98,30));
	/*rPoints.push_back(Point(30,55.98,20));  
	rPoints.push_back(Point(30,55.98,40));*/
	rPoints.push_back(Point(0,30,30));   
	rPoints.push_back(Point(60,26,30)); 
	//rPoints.push_back(Point(7,15,30));   

	rPoints.push_back(Point(15.4,4.019,30));   

	rPoints.push_back(Point(45,4.019,30));   
	rPoints.push_back(Point(30,4.019,30)); //this is not included

	/*rPoints.push_back(Point(30,30,0));   
	rPoints.push_back(Point(30,30,60)); */

	rPoints.push_back(Point(30,15.3,56)); 
	rPoints.push_back(Point(30,25.3,56));

	rPoints.push_back(Point(30,21.3,5)); 
	rPoints.push_back(Point(37, 30 ,7));
	//rPoints.push_back(Point(25,25,5)); 
 }

void forManifoldVertex(Delaunay &dt,vector<pair<Point,int>> &centroids){
	int k=0;
	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){		
		//red --> 0,2,3,4,8,9,11
		//blue --> 1,6,10,19
		//gray--> 5,7,12,13,14,15,16,17,18,19

		/*if( k == 2 || k == 3 || k == 4 || k == 8  || k == 11 || k == 9)	it->info().m_label=1;
		else if( k == 1 || k == 6 || k == 10 || k == 0 )	it->info().m_label=2;
		else it->info().m_label=0;*/

		if( k == 2 || k == 5 || k == 7 || k == 12  || k == 15 || k == 17  || k == 1 || k == 8 || k == 18 || k == 14 || k == 11   )	it->info().m_label=0;
		else if(  k == 10 ||  k == 9  || k == 19  || k == 4  )	it->info().m_label = 1;
		else {
			it->info().m_label=2;				
		}		
		++k;	

		/*it->info().m_label=k%9;
			++k;*/

		//it->info().m_label = k%10;
		//cout<<"inof: "<<it->info().m_label<<endl;
		//k++;
		//cout<<"k: "<<k<<" "<<dt.tetrahedron(it)<<endl;
		Point c = CGAL::centroid(it->vertex(0)->point(), it->vertex(1)->point(),it->vertex(2)->point(),it->vertex(3)->point() );
		centroids.push_back(make_pair(c,it->info().m_label ));
	}
	cout<<"k: "<<k<<endl;
}
void create_star(Delaunay &dt){
	vector<Point> points;
	vector< pair<Point,int> > centroids;
	load_points(points);
	dt.insert(points.begin(),points.end());
	forManifoldVertex(dt,centroids);
}