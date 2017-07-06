
#include <fstream>
#include <sstream>
#include <string>

void saveMesh(Delaunay &dt, const char* filename,string file){
	vector<vector<Cell>> cellByLabel(nroLabels);
	vector<vector<Vertex>> vertexByLabel(nroLabels);

	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++)
		vi->info().m_index=-1;

	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		if(it->info().m_label  < 0) cout<<"label "<<it->info().m_label<<endl;
		cellByLabel[it->info().m_label].push_back(it);
	}
	for(unsigned int i = 0; i < cellByLabel.size(); ++i){
		for(unsigned int j = 0; j < cellByLabel[i].size(); ++j){
			Vertex v0 = cellByLabel[i][j]->vertex(0);
			Vertex v1 = cellByLabel[i][j]->vertex(1);
			Vertex v2 = cellByLabel[i][j]->vertex(2);
			Vertex v3 = cellByLabel[i][j]->vertex(3);
			if(v0->info().m_index != i){
				vertexByLabel[i].push_back(v0);
				v0->info().m_index=i;
			}
			if(v1->info().m_index != i){
				vertexByLabel[i].push_back(v1);
				v1->info().m_index=i;
			}
			if(v2->info().m_index != i){
				vertexByLabel[i].push_back(v2);
				v2->info().m_index=i;
			}
			if(v3->info().m_index != i){
				vertexByLabel[i].push_back(v3);
				v3->info().m_index=i;
			}
		}
	}
	
	for(unsigned int i = 0; i < nroLabels; ++i){
		std::ofstream os;
		std::stringstream ss;
		ss<<filename<<"/"<<file<<i<<".vtk";
		//cout<<ss.str()<<endl;
		os.open( ss.str() );
		os << "# vtk DataFile Version 1.0\nUnstructured Grid Example\nASCII\n\nDATASET UNSTRUCTURED_GRID" << std::endl;
		size_t n,m;
		n = vertexByLabel[i].size();
		m = cellByLabel[i].size();
		os << "POINTS" << ' ' << n << ' ' << "float\n";
		int k = 0;	
		for(unsigned int j = 0; j < vertexByLabel[i].size(); ++j){
			vertexByLabel[i][j]->info().m_index=k;
			os << vertexByLabel[i][j]->point()<< "\n";
			++k;
		}		

		os << "CELLS" << ' ' << m << ' ' << m*5<< '\n';
		for(unsigned int j = 0; j < cellByLabel[i].size(); ++j){
			os <<"4";
            for(unsigned int n=0; n < 4; ++n)
                os << ' ' << cellByLabel[i][j]->vertex(n)->info().m_index;
            os << '\n';
		}

		os << "CELL_TYPES" << ' ' << m << '\n';
		for(unsigned int j = 0; j < cellByLabel[i].size(); ++j)
			os << "10 ";
		os << "\n\n";

		os.close();
	}	

	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++)
		vi->info().m_index=-1;
	cout<<"save finished"<<endl;
}

void savePoints(vector<Point> &points, const char* filename){
	ofstream os(filename);
	os << "# vtk DataFile Version 1.0\nUnstructured Grid Example\nASCII\n\nDATASET UNSTRUCTURED_GRID" << std::endl;
	os << "POINTS" << ' ' << points.size() << ' ' << "float\n";

	for(unsigned int i = 0; i < points.size(); ++i){
		os << points[i]<< "\n";
	}	

	os << "CELLS" << ' ' << points.size() << ' ' << points.size()*2<< '\n';
	for(unsigned int i = 0; i < points.size(); ++i){
		os <<"1"<< ' ' << i<< '\n';
	}
	os << "CELL_TYPES" << ' ' << points.size() << '\n';
	for(unsigned int i = 0; i < points.size(); ++i){
		os << "1 ";
	}
	os << "\n\n";
	os.close();
}

void savePoints2(vector<Point> &circles, const char* filename,double r){
	actor1->GetProperty()->SetOpacity(1);
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkFloatArray> scales = vtkSmartPointer<vtkFloatArray>::New();
    scales->SetName("scales");
    vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");
	for(unsigned int i = 0; i < circles.size(); ++i){	
		points->InsertNextPoint( circles[i].x() , circles[i].y() , circles[i].z() ); // sphere in circle
        scales->InsertNextValue( r );
		colorsp->InsertNextTypedTuple( *(colorsArray + 4) );
	}
	vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->SetPoints(points);
    grid->GetPointData()->AddArray(scales);    
    grid->GetPointData()->SetActiveScalars("scales"); // !!!to set radius first
    grid->GetPointData()->AddArray(colorsp);	
	
	vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
	//polygonSource->GeneratePolygonOff();
	//sphereSource->SetNumberOfSides(100);	

	vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    glyph3D->SetInputData(grid);
    glyph3D->SetSourceConnection(sphereSource->GetOutputPort());

	vtkSmartPointer<vtkPolyDataWriter> writer =    vtkSmartPointer<vtkPolyDataWriter>::New();
	//writer->SetInput(glyph2D->GetOutputPort() );
	writer->SetInputConnection(glyph3D->GetOutputPort());
	writer->SetFileName(filename);
	writer->Write();
}

int map_edges(map<pair<int,int>,int> &m, pair<int,int> & p){
    std::map<pair<int,int>,int>::const_iterator it = m.find( p );
    if (it == m.end())
        return -1;
    return it->second;
}

void save_edges_of_vertices_star(Delaunay &dt,vector<Vertex> vertices, const char* filename){
	vtkSmartPointer<vtkFloatArray> linescales = vtkSmartPointer<vtkFloatArray>::New();
    linescales->SetName("linescales");

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();

	int k=0;		
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) { 
		points->InsertNextPoint(vi->point().x(), vi->point().y(), vi->point().z());
		vi->info().m_index=k;
		++k;
	} 
	vector<Cell> involved_cells;
	for( unsigned int i = 0; i < vertices.size(); ++i){
		Vertex v = vertices[i];
		list<Cell> cells;
		dt.finite_incident_cells(v,back_inserter(cells));
		for(std::list<Cell>::iterator  ite = cells.begin(); ite != cells.end(); ite++){ 
			Cell c = *ite;
			involved_cells.push_back(c);
		}
	}

	map<pair<int,int>,int> map_e;
	for( unsigned int i = 0; i < involved_cells.size(); ++i){
		Cell c = involved_cells[i];
		int source,target,map_p,map_q;
		pair<int,int> p,q;

		source = c->vertex(0)->info().m_index;
		target = c->vertex(1)->info().m_index;
		p = make_pair(source,target);
		q = make_pair(target,source);
		map_p = map_edges(map_e,p);
		map_q = map_edges(map_e,q);
		if( map_p == -1 && map_q == -1){
			line->GetPointIds()->SetId(0,source); 
			line->GetPointIds()->SetId(1,target);
			lines->InsertNextCell(line);
			linescales->InsertNextValue(1);
			map_e[p] = 1;
			map_e[q] = 1;
		}	

		source = c->vertex(0)->info().m_index;
		target = c->vertex(2)->info().m_index;
		p = make_pair(source,target);
		q = make_pair(target,source);
		map_p = map_edges(map_e,p);
		map_q = map_edges(map_e,q);
		if( map_p == -1 && map_q == -1){
			line->GetPointIds()->SetId(0,source); 
			line->GetPointIds()->SetId(1,target);
			lines->InsertNextCell(line);
			linescales->InsertNextValue(1);
			map_e[p] = 1;
			map_e[q] = 1;
		}	

		source = c->vertex(0)->info().m_index;
		target = c->vertex(3)->info().m_index;
		p = make_pair(source,target);
		q = make_pair(target,source);
		map_p = map_edges(map_e,p);
		map_q = map_edges(map_e,q);
		if( map_p == -1 && map_q == -1){
			line->GetPointIds()->SetId(0,source); 
			line->GetPointIds()->SetId(1,target);
			lines->InsertNextCell(line);
			linescales->InsertNextValue(1);
			map_e[p] = 1;
			map_e[q] = 1;
		}	

		source = c->vertex(1)->info().m_index;
		target = c->vertex(2)->info().m_index;
		p = make_pair(source,target);
		q = make_pair(target,source);
		map_p = map_edges(map_e,p);
		map_q = map_edges(map_e,q);
		if( map_p == -1 && map_q == -1){
			line->GetPointIds()->SetId(0,source); 
			line->GetPointIds()->SetId(1,target);
			lines->InsertNextCell(line);
			linescales->InsertNextValue(1);
			map_e[p] = 1;
			map_e[q] = 1;
		}	

		source = c->vertex(1)->info().m_index;
		target = c->vertex(3)->info().m_index;
		p = make_pair(source,target);
		q = make_pair(target,source);
		map_p = map_edges(map_e,p);
		map_q = map_edges(map_e,q);
		if( map_p == -1 && map_q == -1){
			line->GetPointIds()->SetId(0,source); 
			line->GetPointIds()->SetId(1,target);
			lines->InsertNextCell(line);
			linescales->InsertNextValue(1);
			map_e[p] = 1;
			map_e[q] = 1;
		}	

		source = c->vertex(2)->info().m_index;
		target = c->vertex(3)->info().m_index;
		p = make_pair(source,target);
		q = make_pair(target,source);
		map_p = map_edges(map_e,p);
		map_q = map_edges(map_e,q);
		if( map_p == -1 && map_q == -1){
			line->GetPointIds()->SetId(0,source); 
			line->GetPointIds()->SetId(1,target);
			lines->InsertNextCell(line);
			linescales->InsertNextValue(1);
			map_e[p] = 1;
			map_e[q] = 1;
		}	
	}

	vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
	linesPolyData->SetPoints(points);
	linesPolyData->SetLines(lines);	
	linesPolyData->GetCellData()->SetScalars(linescales);	

	vtkSmartPointer<vtkPolyDataWriter> writer2 =    vtkSmartPointer<vtkPolyDataWriter>::New();
	writer2->SetInputData(linesPolyData);	
	//writer2->SetFileName("multimaterial/edges_1.vtk");
	writer2->SetFileName(filename);
	writer2->Write();


	//end code
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++)
		vi->info().m_index=-1;
}

void visite_cells2(Delaunay &dt,Cell &c, Vertex &v,int &label,vector<pair<Point,Point>> &edges){
	vector<Cell> visitedCells;
	visitedCells.push_back(c);
	c->info().visited = true;
	Cell vc,c0,c1,c2;
	while(visitedCells.size() != 0){
		vc = visitedCells.back();
		visitedCells.pop_back();
		int index = vc->index(v);
		for(unsigned int i = 1; i < 4 ;++i){
			c0 = vc->neighbor( (index + i)%4 );
			Point a = CGAL::centroid(vc->vertex(0)->point(), vc->vertex(1)->point(),vc->vertex(2)->point(), vc->vertex(3)->point());
			Point b = CGAL::centroid(c0->vertex(0)->point(), c0->vertex(1)->point(),c0->vertex(2)->point(), c0->vertex(3)->point());
			if(!c0->info().visited && c0->info().m_label == label ){
				visitedCells.push_back(c0);
				c0->info().visited = true;				
				edges.push_back(make_pair(a,b));
			}
			else if(c0->info().m_label == label){
				if( ! ( std::find(edges.begin(), edges.end(), make_pair(a,b)  )!= edges.end() ) )
					edges.push_back(make_pair(a,b));
				if( ! ( std::find(edges.begin(), edges.end(), make_pair(b,a)  )!= edges.end() ) )
					edges.push_back(make_pair(b,a));
			}
		}
	}
}

void get_edges_graph(Delaunay &dt,Vertex v, vector<pair<Point,Point>> &edges, int label){
	list<Cell> incident;
	dt.finite_incident_cells(v,back_inserter(incident));
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		if(! (*it)->info().visited && (*it)->info().m_label == label ) {
			visite_cells2(dt, (*it) ,v,(*it)->info().m_label,edges);
			(*it)->info().visited=1;
		}
	}	
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++)
		(*it)->info().visited=0;
}

void save_edges_graph(Delaunay &dt, vector<pair<Point,Point>> &edges, const char* filename){
	vtkSmartPointer<vtkFloatArray> linescales = vtkSmartPointer<vtkFloatArray>::New();
    linescales->SetName("linescales");

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();

	int k=0;		
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) { 
		points->InsertNextPoint(vi->point().x(), vi->point().y(), vi->point().z());
		vi->info().m_index=k;
		++k;
	} 
	vector<pair<int,int>> indexs;
	for(unsigned int i = 0; i < edges.size(); ++i){
		points->InsertNextPoint(edges[i].first.x(), edges[i].first.y(), edges[i].first.z());
		int f=k;
		++k;
		points->InsertNextPoint(edges[i].second.x(), edges[i].second.y(), edges[i].second.z());
		int s=k;
		++k;
		indexs.push_back(make_pair(f,s));
	}

	for(unsigned int i = 0; i < indexs.size(); ++i){
		line->GetPointIds()->SetId(0,indexs[i].first); 
		line->GetPointIds()->SetId(1,indexs[i].second);
		lines->InsertNextCell(line);
		linescales->InsertNextValue(10);
	}

	vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
	linesPolyData->SetPoints(points);
	linesPolyData->SetLines(lines);	
	linesPolyData->GetCellData()->SetScalars(linescales);	

	vtkSmartPointer<vtkPolyDataWriter> writer2 =    vtkSmartPointer<vtkPolyDataWriter>::New();
	writer2->SetInputData(linesPolyData);	
	//writer2->SetFileName("multimaterial/edges_1.vtk");
	writer2->SetFileName(filename);
	writer2->Write();
	//end code
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++)
		vi->info().m_index=-1;
}

void save_spheres_graph(Delaunay & dt, vector<Point> &centroids, const char* filename){
	vtkSmartPointer<vtkPoints> pointsVertex =   vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");

	vtkSmartPointer<vtkFloatArray> scales = vtkSmartPointer<vtkFloatArray>::New();
    scales->SetName("scales");

	for(unsigned int i=0;i<centroids.size();++i){
		pointsVertex->InsertNextPoint (centroids[i].x(), centroids[i].y(), centroids[i].z());
		//colorsp->InsertNextTypedTuple( *(colorsArray + centroids[i].second) );
		colorsp->InsertNextTypedTuple( grey );
		scales->InsertNextValue(0.001); 
	}

	vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(pointsVertex);
	polydata->GetPointData()->AddArray(scales); 
	polydata->GetPointData()->SetActiveScalars("scales");
	polydata->GetPointData()->SetScalars(colorsp);

	vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
	glyph3D->SetColorModeToColorByScalar();
	glyph3D->SetScaleModeToScaleByScalar();
	glyph3D->SetScaleFactor(0.02);
	glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
	glyph3D->SetInputData(polydata);		
	//glyph3D->ScalingOn();
	glyph3D->Update();

	vtkSmartPointer<vtkPolyDataWriter> writer3 =    vtkSmartPointer<vtkPolyDataWriter>::New();
	writer3->SetInputConnection(glyph3D->GetOutputPort());
	writer3->SetFileName(filename);
	writer3->Write();
}

void save_partial_mesh(Delaunay &dt,vector<Cell> &cells, const char* filename,string file){
	cout<<"size cells : "<<cells.size()<<endl;
	vector<vector<Cell>> cellByLabel(nroLabels);
	vector<vector<Vertex>> vertexByLabel(nroLabels);

	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++)
		vi->info().m_index=-1;

	for(unsigned int i = 0; i < cells.size(); ++i){
		Cell it = cells[i];
		if(it->info().m_label  < 0) cout<<"label "<<it->info().m_label<<endl;
		cellByLabel[it->info().m_label].push_back(it);
	}
	for(unsigned int i = 0; i < cellByLabel.size(); ++i){
		for(unsigned int j = 0; j < cellByLabel[i].size(); ++j){
			Vertex v0 = cellByLabel[i][j]->vertex(0);
			Vertex v1 = cellByLabel[i][j]->vertex(1);
			Vertex v2 = cellByLabel[i][j]->vertex(2);
			Vertex v3 = cellByLabel[i][j]->vertex(3);
			if(v0->info().m_index != i){
				vertexByLabel[i].push_back(v0);
				v0->info().m_index=i;
			}
			if(v1->info().m_index != i){
				vertexByLabel[i].push_back(v1);
				v1->info().m_index=i;
			}
			if(v2->info().m_index != i){
				vertexByLabel[i].push_back(v2);
				v2->info().m_index=i;
			}
			if(v3->info().m_index != i){
				vertexByLabel[i].push_back(v3);
				v3->info().m_index=i;
			}
		}
	}
	
	for(unsigned int i = 0; i < nroLabels; ++i){
		std::ofstream os;
		std::stringstream ss;
		ss<<filename<<"/"<<file<<i<<".vtk";
		//ss<<"identification/tony"<<file<<i<<".vtk";
		
		//cout<<ss.str()<<endl;
		os.open( ss.str() );
		os << "# vtk DataFile Version 1.0\nUnstructured Grid Example\nASCII\n\nDATASET UNSTRUCTURED_GRID" << std::endl;
		size_t n,m;
		n = vertexByLabel[i].size();
		m = cellByLabel[i].size();
		os << "POINTS" << ' ' << n << ' ' << "float\n";
		int k = 0;	
		for(unsigned int j = 0; j < vertexByLabel[i].size(); ++j){
			vertexByLabel[i][j]->info().m_index=k;
			os << vertexByLabel[i][j]->point()<< "\n";
			++k;
		}		

		os << "CELLS" << ' ' << m << ' ' << m*5<< '\n';
		for(unsigned int j = 0; j < cellByLabel[i].size(); ++j){
			os <<"4";
            for(unsigned int n=0; n < 4; ++n)
                os << ' ' << cellByLabel[i][j]->vertex(n)->info().m_index;
            os << '\n';
		}

		os << "CELL_TYPES" << ' ' << m << '\n';
		for(unsigned int j = 0; j < cellByLabel[i].size(); ++j)
			os << "10 ";
		os << "\n\n";

		os.close();
	}	

	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++)
		vi->info().m_index=-1;
	cout<<"save finished"<<endl;
}

vector<Cell> cells_has_vertex(Delaunay &dt, vector<Vertex> &vertices){
	vector<Cell> cells_v;
	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		for( unsigned int i = 0; i < vertices.size(); ++i){
			Vertex v = vertices[i];
			if( (*it).has_vertex(v) ){
				cells_v.push_back(it);
				break;
			}
		}
	}
	return cells_v;
}

vector<Point> get_centroids_by_label(Delaunay &dt,Vertex v , int label){
	list<Cell> incident;
	vector<Point> centroids;
	dt.finite_incident_cells(v,back_inserter(incident));
	for(std::list<Cell>::iterator it = incident.begin(); it != incident.end(); it++){
		Cell c= *it;
		if( c->info().m_label == label){
			Point p = CGAL::centroid(c->vertex(0)->point(), c->vertex(1)->point(),c->vertex(2)->point(),c->vertex(3)->point() );
			centroids.push_back( p);
		}
	}
	return centroids;
}

void save_vertex_star(Delaunay &dt, vector<Vertex> &vertices, const char* filename,string file){
	std::stringstream c_dir,e_dir;
	c_dir<<filename<<"/"<<file<<".vtk";
	e_dir<<filename<<"/"<<file<<"_edges.vtk";
	Point p = Point(24,23,30);
	Vertex v;
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) { 
		if( vi->point() == p){
			v=vi;break;
		}
	} 
	/*vector<Vertex> vertices;
	vertices.push_back(v);*/

	vector<Cell> cells_v = cells_has_vertex(dt,vertices);	
	save_partial_mesh(dt,cells_v,filename,file);	// original
	save_edges_of_vertices_star(dt,vertices,( e_dir.str() ).c_str());

	/*for(unsigned int i = 0; i < 3; ++i){
		std::stringstream ss;
		ss<<"identification/original/edges"<<i<<".vtk";
		std::stringstream ee;
		ee<<"identification/original/spheres"<<i<<".vtk";
		vector<pair<Point,Point>> edges;
		get_edges_graph(dt,v,edges,i);
		vector<Point> centroids = get_centroids_by_label(dt,v,i);
		save_edges_graph(dt,edges,( ss.str() ).c_str() ) ;
		save_spheres_graph(dt,centroids,( ee.str() ).c_str()  );
	}*/
	cout<<"finishing save star ..."<<endl;
	//getchar();
}