
#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

#include <iostream>
#include <vector>

using namespace std; 

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkDataSetMapper.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLine.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSphereSource.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkProperty.h>

void drawMesh(Delaunay & dt, vector<Facet> &boundary, vector<Point> &v){
	unsigned char grey[3] = {192,192,192};
	unsigned char blue[3] = {0, 0, 255};
	unsigned char red[3] = {255, 0, 0};	
	unsigned char green[3] = {0, 255, 0};
	unsigned char purple[3] = {153,0,153};
	unsigned char bluewhite[3] = {0, 255, 255};
	unsigned char brown[3] = {153, 76, 0};	
	unsigned char black[3] = {0, 0, 0};
	unsigned char pink[3] = {255, 102, 255};
	unsigned char *colorsArray[9] = {grey,red,blue,green,purple,bluewhite,brown,black,pink};

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int k=0;		
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) { 
		points->InsertNextPoint(vi->point().x(), vi->point().y(), vi->point().z());
		vi->info().m_index=k;
		++k;
	} 

	vtkSmartPointer<vtkUnsignedCharArray> colors =  vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid2 =    vtkSmartPointer<vtkUnstructuredGrid>::New();
	unstructuredGrid2->SetPoints(points);
	vtkSmartPointer<vtkCellArray> tetrahedrons = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();


	//just lines of an label
	/*vtkSmartPointer<vtkUnsignedCharArray> colorLines = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorLines->SetNumberOfComponents(3);
	colorLines->SetName("ColorsLines");

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();*/

	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		if(it->info().m_label > 0 ){
			//cout<<"cell to draw"<<endl;
			colors->InsertNextTypedTuple( *(colorsArray + it->info().m_label) );
			tetra->GetPointIds()->SetId(0, it->vertex(0)->info().m_index);
			tetra->GetPointIds()->SetId(1, it->vertex(1)->info().m_index);
			tetra->GetPointIds()->SetId(2, it->vertex(2)->info().m_index);
			tetra->GetPointIds()->SetId(3, it->vertex(3)->info().m_index);
			tetrahedrons->InsertNextCell(tetra);
		}
	}

	unstructuredGrid2->SetCells(VTK_TETRA, tetrahedrons);
	unstructuredGrid2->GetCellData()->SetScalars(colors);

	//draw all lines
	vtkSmartPointer<vtkUnsignedCharArray> colorLines = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorLines->SetNumberOfComponents(3);
	colorLines->SetName("ColorsLines");

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();

	for(Finite_edges_iterator ei=dt.finite_edges_begin();ei!=dt.finite_edges_end(); ei++){
		Vertex vs = ei->first->vertex(ei->second);
		Vertex vt = ei->first->vertex(ei->third);
		line->GetPointIds()->SetId(0,vs->info().m_index); 
		line->GetPointIds()->SetId(1,vt->info().m_index);
		lines->InsertNextCell(line);
		colorLines->InsertNextTypedTuple(*(colorsArray + 7));
	}

	vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
	linesPolyData->SetPoints(points);
	linesPolyData->SetLines(lines);
	linesPolyData->GetCellData()->SetScalars(colorLines);
 
	//draw surface
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
	for(unsigned int i=0;i<boundary.size();++i){
		triangle->GetPointIds()->SetId(0, boundary[i].first->vertex( (boundary[i].second+1)%4 )->info().m_index);
		triangle->GetPointIds()->SetId(1, boundary[i].first->vertex( (boundary[i].second+2)%4 )->info().m_index);
		triangle->GetPointIds()->SetId(2, boundary[i].first->vertex( (boundary[i].second+3)%4 )->info().m_index);
		triangles->InsertNextCell(triangle);
	}
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	polydata->SetPolys(triangles);
	vtkSmartPointer<vtkPolyDataMapper> mapper3 =  vtkSmartPointer<vtkPolyDataMapper>::New();
	
	//draw points 
	vtkSmartPointer<vtkPoints> pointsVertex =   vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");
	for(unsigned int i=0;i<v.size();++i){
		pointsVertex->InsertNextPoint (v[i].x(), v[i].y(), v[i].z());
		colorsp->InsertNextTypedTuple( *(colorsArray + i%9 )  );
	}

	vtkSmartPointer<vtkPolyData> pointsPolydata =   vtkSmartPointer<vtkPolyData>::New(); 
	pointsPolydata->SetPoints(pointsVertex);

	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =    vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexFilter->SetInputData(pointsPolydata);
	vertexFilter->Update();

	vtkSmartPointer<vtkPolyData> polydatapc =    vtkSmartPointer<vtkPolyData>::New();
	polydatapc->ShallowCopy(vertexFilter->GetOutput());
	polydatapc->GetPointData()->SetScalars(colorsp);
	vtkSmartPointer<vtkPolyDataMapper> mapper4 =    vtkSmartPointer<vtkPolyDataMapper>::New();
	// Create a mapper and actor
	vtkSmartPointer<vtkDataSetMapper> mapper =   vtkSmartPointer<vtkDataSetMapper>::New();
	vtkSmartPointer<vtkPolyDataMapper> mapper2 =  vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(unstructuredGrid2->GetProducerPort());
	mapper2->SetInputConnection(linesPolyData->GetProducerPort());
	mapper3->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputData(unstructuredGrid2);
	mapper2->SetInputData(linesPolyData);
	mapper3->SetInputData(polydata);
	mapper4->SetInputData(polydatapc);
#endif
 
	vtkSmartPointer<vtkActor> actor1 =  vtkSmartPointer<vtkActor>::New();
	actor1->SetMapper(mapper);
 
  	vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
	actor2->SetMapper(mapper2);

	vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
	actor3->SetMapper(mapper3);

	vtkSmartPointer<vtkActor> actor4 = vtkSmartPointer<vtkActor>::New();
	actor4->SetMapper(mapper4);
	actor4->GetProperty()->SetPointSize(10);

	vtkSmartPointer<vtkRenderer> renderer =  vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =  vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
 
	renderer->AddActor(actor1);//tetrahedrons
	//renderer->AddActor(actor2);//lines
	//renderer->AddActor(actor3);//surface
	renderer->AddActor(actor4);//points

	renderer->SetBackground(1, 1, 1);
	renderWindow->Render();
	//renderWindowInteractor->Start();
}

void drawSurface(Polyhedron &P, vector<Vector> &vertices){
	unsigned char grey[3] = {192,192,192};
	unsigned char red[3] = {255, 0, 0};
	unsigned char blue[3] = {0, 0, 255};
	unsigned char green[3] = {0, 255, 0};
	unsigned char purple[3] = {153,0,153};
	unsigned char bluewhite[3] = {0, 255, 255};
	unsigned char brown[3] = {153, 76, 0};	
	unsigned char black[3] = {0, 0, 0};
	unsigned char pink[3] = {255, 102, 255};
	unsigned char *colorsArray[9] = {grey,red,blue,green,purple,bluewhite,brown,black,pink};

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int k=0;
	for ( Vertex_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v){
		points->InsertNextPoint(v->point().x(), v->point().y(), v->point().z() );
		v->m_index=k;
		++k;
	}
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
	for ( Facet_iterator f = P.facets_begin(); f != P.facets_end(); ++f){
		triangle->GetPointIds()->SetId(0, f->halfedge()->vertex()->m_index);
		triangle->GetPointIds()->SetId(1, f->halfedge()->next()->vertex()->m_index);
		triangle->GetPointIds()->SetId(2, f->halfedge()->next()->next()->vertex()->m_index );
		triangles->InsertNextCell(triangle);
	}
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	polydata->SetPolys(triangles);
	vtkSmartPointer<vtkPolyDataMapper> mapper1 =  vtkSmartPointer<vtkPolyDataMapper>::New();

	// draw points
	vtkSmartPointer<vtkPoints> pointsVertex =   vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");
	/*for(int i=0;i<v.size();++i){
		pointsVertex->InsertNextPoint (v[i]->point().x(), v[i]->point().y(), v[i]->point().z());
		colorsp->InsertNextTypedTuple( *(colorsArray + i%9) );
	}*/
	for(unsigned int i=0;i<vertices.size();++i){
		pointsVertex->InsertNextPoint (vertices[i].x(), vertices[i].y(), vertices[i].z());
		colorsp->InsertNextTypedTuple( *(colorsArray + i%9) );
	}

	vtkSmartPointer<vtkPolyData> pointsPolydata =   vtkSmartPointer<vtkPolyData>::New(); 
	pointsPolydata->SetPoints(pointsVertex);

	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =    vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexFilter->SetInputData(pointsPolydata);
	vertexFilter->Update();

	vtkSmartPointer<vtkPolyData> polydatapc =    vtkSmartPointer<vtkPolyData>::New();
	polydatapc->ShallowCopy(vertexFilter->GetOutput());
	polydatapc->GetPointData()->SetScalars(colorsp);
	vtkSmartPointer<vtkPolyDataMapper> mapper2 =    vtkSmartPointer<vtkPolyDataMapper>::New();

#if VTK_MAJOR_VERSION <= 5
	mapper1->SetInputConnection(polydata->GetProducerPort());
#else
	mapper1->SetInputData(polydata);
	mapper2->SetInputData(polydatapc);
#endif
	vtkSmartPointer<vtkActor> actor1 =  vtkSmartPointer<vtkActor>::New();
	actor1->SetMapper(mapper1);

	vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
	actor2->SetMapper(mapper2);
	actor2->GetProperty()->SetPointSize(10);

	vtkSmartPointer<vtkRenderer> renderer =  vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =  vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
 
	renderer->AddActor(actor1);
	renderer->AddActor(actor2);
	renderer->SetBackground(.3, .6, .3);
	renderWindow->Render();
	renderWindowInteractor->Start();
}

void drawMesh2(Delaunay & dt, vector<Facet> &boundary, vector<Point> &v){
	unsigned char grey[3] = {192,192,192};
	unsigned char red[3] = {255, 0, 0};
	unsigned char blue[3] = {0, 0, 255};
	unsigned char green[3] = {0, 255, 0};
	unsigned char purple[3] = {153,0,153};
	unsigned char bluewhite[3] = {0, 255, 255};
	unsigned char brown[3] = {153, 76, 0};	
	unsigned char black[3] = {0, 0, 0};
	unsigned char pink[3] = {255, 102, 255};
	unsigned char *colorsArray[9] = {grey,red,blue,green,purple,bluewhite,brown,black,pink};

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int k=0;		
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) { 
		points->InsertNextPoint(vi->point().x(), vi->point().y(), vi->point().z());
		vi->info().m_index=k;
		++k;
	} 

	vtkSmartPointer<vtkUnsignedCharArray> colors =  vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid2 =    vtkSmartPointer<vtkUnstructuredGrid>::New();
	unstructuredGrid2->SetPoints(points);
	vtkSmartPointer<vtkCellArray> tetrahedrons = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();


	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		//it->info().m_label=0;
		if(it->info().m_label > 0){
			colors->InsertNextTypedTuple( *(colorsArray + it->info().m_label) );
			tetra->GetPointIds()->SetId(0, it->vertex(0)->info().m_index);
			tetra->GetPointIds()->SetId(1, it->vertex(1)->info().m_index);
			tetra->GetPointIds()->SetId(2, it->vertex(2)->info().m_index);
			tetra->GetPointIds()->SetId(3, it->vertex(3)->info().m_index);
			tetrahedrons->InsertNextCell(tetra);
		}
	}

	unstructuredGrid2->SetCells(VTK_TETRA, tetrahedrons);
	unstructuredGrid2->GetCellData()->SetScalars(colors);

	 
	//draw surface
	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
	for(unsigned int i=0;i<boundary.size();++i){
		triangle->GetPointIds()->SetId(0, boundary[i].first->vertex( (boundary[i].second+1)%4 )->info().m_index);
		triangle->GetPointIds()->SetId(1, boundary[i].first->vertex( (boundary[i].second+2)%4 )->info().m_index);
		triangle->GetPointIds()->SetId(2, boundary[i].first->vertex( (boundary[i].second+3)%4 )->info().m_index);
		triangles->InsertNextCell(triangle);
	}
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	polydata->SetPolys(triangles);
	vtkSmartPointer<vtkPolyDataMapper> mapper3 =  vtkSmartPointer<vtkPolyDataMapper>::New();
	
	//draw points 
	vtkSmartPointer<vtkPoints> pointsVertex =   vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");
	for(unsigned int i=0;i<v.size();++i){
		pointsVertex->InsertNextPoint (v[i].x(), v[i].y(), v[i].z());
		colorsp->InsertNextTypedTuple( *(colorsArray + i%9) );
	}

	vtkSmartPointer<vtkPolyData> pointsPolydata =   vtkSmartPointer<vtkPolyData>::New(); 
	pointsPolydata->SetPoints(pointsVertex);

	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =    vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexFilter->SetInputData(pointsPolydata);
	vertexFilter->Update();

	vtkSmartPointer<vtkPolyData> polydatapc =    vtkSmartPointer<vtkPolyData>::New();
	polydatapc->ShallowCopy(vertexFilter->GetOutput());
	polydatapc->GetPointData()->SetScalars(colorsp);
	vtkSmartPointer<vtkPolyDataMapper> mapper4 =    vtkSmartPointer<vtkPolyDataMapper>::New();
	// Create a mapper and actor
	vtkSmartPointer<vtkDataSetMapper> mapper =   vtkSmartPointer<vtkDataSetMapper>::New();


	//draw sphere 

	vtkSmartPointer<vtkSphereSource> sphereSource =  vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->SetCenter(30, 30, 30);
	sphereSource->SetRadius(2.5);

	vtkSmartPointer<vtkPolyDataMapper> mapper5 =  vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper5->SetInputConnection(sphereSource->GetOutputPort());
 
	vtkSmartPointer<vtkActor> actorSphere = 	vtkSmartPointer<vtkActor>::New();
	actorSphere->SetMapper(mapper5);

#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(unstructuredGrid2->GetProducerPort());
	mapper2->SetInputConnection(linesPolyData->GetProducerPort());
	mapper3->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputData(unstructuredGrid2);
	mapper3->SetInputData(polydata);
	mapper4->SetInputData(polydatapc);
#endif
 
	vtkSmartPointer<vtkActor> actor1 =  vtkSmartPointer<vtkActor>::New();
	actor1->SetMapper(mapper);
 
	vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
	actor3->SetMapper(mapper3);

	vtkSmartPointer<vtkActor> actor4 = vtkSmartPointer<vtkActor>::New();
	actor4->SetMapper(mapper4);
	actor4->GetProperty()->SetPointSize(10);

	vtkSmartPointer<vtkRenderer> renderer =  vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =  vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
 
	renderer->AddActor(actor1);//tetrahedrons
	//renderer->AddActor(actor2);//lines
	//renderer->AddActor(actorSphere);
	//renderer->AddActor(actor3);//surface
	renderer->AddActor(actor4);//points

	//renderer->SetBackground(.3, .6, .3);
	renderer->SetBackground(1, 1, 1);
	renderWindow->Render();
	renderWindowInteractor->Start();
}