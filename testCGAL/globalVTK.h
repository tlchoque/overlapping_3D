#define vtkRenderingCore_AUTOINIT 4(vtkRenderingOpenGL2, vtkInteractionStyle,vtkRenderingVolumeOpenGL2,vtkRenderingFreeType)

#include <Windows.h>// for threads

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
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkElevationFilter.h>
#include <vtkVectorText.h>
#include <vtkCommand.h>
#include <vtkAutoInit.h>

#include "vtkShrinkFilter.h"
#include "vtkStructuredGridReader.h"
#include "vtkExtractEdges.h"
#include "vtkGenericCell.h"
#include "vtkPolyhedron.h"
#include "vtkCubeSource.h"
#include "vtkIdList.h"
#include "vtkDataArray.h"
#include "vtkPointLocator.h"

#include <vtkRegularPolygonSource.h>
#include <vtkFloatArray.h>
#include <vtkGlyph3D.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPolyDataWriter.h>


void* handler = 0;
void* mtx;

#define VTK_CRT(type, name) \
    vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

vtkSmartPointer<vtkRenderer> renderer =  vtkSmartPointer<vtkRenderer>::New();
vtkSmartPointer<vtkRenderWindow> renderWindow =  vtkSmartPointer<vtkRenderWindow>::New();	
vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
vtkInteractorStyleTrackballCamera *styleDisplay =         vtkInteractorStyleTrackballCamera::New();
vtkSmartPointer<vtkRenderWindowInteractor> irenDisplay = vtkSmartPointer<vtkRenderWindowInteractor>::New();

vtkSmartPointer<vtkDataSetMapper> mapper1 =  vtkSmartPointer<vtkDataSetMapper>::New();//for tetrahedrons
vtkSmartPointer<vtkPolyDataMapper> mapper2 =  vtkSmartPointer<vtkPolyDataMapper>::New();//for triangles
vtkSmartPointer<vtkPolyDataMapper> mapper3 =  vtkSmartPointer<vtkPolyDataMapper>::New();//for points
vtkSmartPointer<vtkDataSetMapper> mapper4 =  vtkSmartPointer<vtkDataSetMapper>::New();//for set cells around a point
vtkSmartPointer<vtkDataSetMapper> mapper5 =  vtkSmartPointer<vtkDataSetMapper>::New();//for set cells or edges


vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
vtkSmartPointer<vtkActor> actor4 = vtkSmartPointer<vtkActor>::New();
vtkSmartPointer<vtkActor> actor5 = vtkSmartPointer<vtkActor>::New();

unsigned char grey[3] = {192,192,192};//0
unsigned char red[3] = {255, 0, 0};	//1
unsigned char blue[3] = {0, 0, 255};//2
unsigned char green[3] = {0, 255, 0};//3
unsigned char purple[3] = {153,0,153};//4
unsigned char bluewhite[3] = {0, 255, 255};//5
unsigned char brown[3] = {153, 76, 0};	//6
unsigned char black[3] = {0, 0, 0}; //7
unsigned char pink[3] = {255, 102, 255};//8
unsigned char orange[3] = {249, 186, 59};//9
unsigned char *colorsArray[10] = {grey,red,blue,green,purple,bluewhite,brown,black,pink,orange};

class CommandSubclass2 : public vtkCommand{
public:
    vtkTypeMacro(CommandSubclass2, vtkCommand);
    static CommandSubclass2 *New()    {
        return new CommandSubclass2;
    }
    void Execute(vtkObject *caller, unsigned long vtkNotUsed(eventId), 
        void *vtkNotUsed(callData))    {
        vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);
        iren->Render();
    }
};

unsigned long __stdcall displayVTK(void* param){
    renderWindow->SetSize(800, 600);
    renderWindow->AddRenderer(renderer);
    renderWindow->Render();
    irenDisplay->SetRenderWindow(renderWindow);
    irenDisplay->SetInteractorStyle(styleDisplay);// move slowly
    irenDisplay->Initialize();
    vtkSmartPointer<CommandSubclass2> timerCallback = vtkSmartPointer<CommandSubclass2>::New();
    irenDisplay->AddObserver ( vtkCommand::TimerEvent, timerCallback );
    irenDisplay->CreateRepeatingTimer(100);
    irenDisplay->Start();
    return 0;
}

void drawMesh(Delaunay & dt, vector<Point> &vertices,bool sa);

void drawMesh(Delaunay & dt, vector<Facet> &boundary, vector<Point> &v,bool sa);

void drawMesh(Delaunay & dt, bool sa);

void begin_VTK(Delaunay & dt){
	unsigned long id_thread;
    VTK_CRT(vtkVectorText, text);
    text->SetText("...");
    VTK_CRT(vtkElevationFilter, elevation);
    elevation->SetInputConnection(text->GetOutputPort());
	mapper1->SetInputConnection(elevation->GetOutputPort());
    mapper1->Update();	
	
    actor1->SetMapper(mapper1);
	//actor1->GetProperty()->SetOpacity(0.7);
	actor2->SetMapper(mapper2);
	actor3->SetMapper(mapper3);
	actor3->GetProperty()->SetPointSize(20);
	actor4->SetMapper(mapper4);

	actor5->SetMapper(mapper5);
	actor5->GetProperty()->SetLineWidth(4);

    renderer->AddActor(actor1);
	renderer->AddActor(actor2);
	renderer->AddActor(actor3);
	renderer->AddActor(actor4);
	renderer->AddActor(actor5);
    renderer->SetBackground(1, 1, 1);	
    handler = CreateThread(0, 0, displayVTK, 0, 0, &id_thread);
    if(!handler){
        printf("Cannot create thread. Error code = %d\n", GetLastError());
        getchar();
        return;
    }	
}

void drawMesh(Delaunay & dt){
	vector<Facet> boundary;
	vector<Point> v;
	drawMesh(dt,boundary,v,0);
}

void drawMesh(Delaunay & dt, vector<Facet> &boundary) {
	vector<Point> v;
	drawMesh(dt, boundary, v, 0);
}

void drawMesh(Delaunay & dt, vector<Point> &v ){
	vector<Facet> boundary;
	drawMesh(dt,boundary,v,0);
}

void drawMesh(Delaunay & dt, vector<Facet> &boundary, vector<Point> &v) {
	drawMesh(dt, boundary, v, 0);
}

void drawMesh(Delaunay & dt, vector<Point> &vertices,bool sa){
	vector<Facet> boundary;
	drawMesh(dt,boundary,vertices,sa);
}

void drawMesh(Delaunay & dt, bool sa){
	vector<Point> vertices;
	vector<Facet> boundary;
	drawMesh(dt,boundary,vertices,sa);
}

void drawMesh(Delaunay & dt, vector<Facet> &boundary, vector<Point> &v,bool sa){
	actor1->GetProperty()->SetOpacity(0.3);
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
		if( sa ){
			if(it->info().m_solution > 0 ){
				colors->InsertNextTypedTuple( *(colorsArray + it->info().m_solution) );//changed
				tetra->GetPointIds()->SetId(0, it->vertex(0)->info().m_index);
				tetra->GetPointIds()->SetId(1, it->vertex(1)->info().m_index);
				tetra->GetPointIds()->SetId(2, it->vertex(2)->info().m_index);
				tetra->GetPointIds()->SetId(3, it->vertex(3)->info().m_index);
				tetrahedrons->InsertNextCell(tetra);
			}
		}
		else{
			if(it->info().m_label > 0){
				colors->InsertNextTypedTuple( *(colorsArray + it->info().m_label) );//changed
				tetra->GetPointIds()->SetId(0, it->vertex(0)->info().m_index);
				tetra->GetPointIds()->SetId(1, it->vertex(1)->info().m_index);
				tetra->GetPointIds()->SetId(2, it->vertex(2)->info().m_index);
				tetra->GetPointIds()->SetId(3, it->vertex(3)->info().m_index);
				tetrahedrons->InsertNextCell(tetra);
			}
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
	
	//draw points 
	vtkSmartPointer<vtkPoints> pointsVertex =   vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");
	for(unsigned int i=0;i<v.size();++i){
		pointsVertex->InsertNextPoint (v[i].x(), v[i].y(), v[i].z());
		colorsp->InsertNextTypedTuple( *(colorsArray + i%10 )  );
	}
	vtkSmartPointer<vtkPolyData> pointsPolydata =   vtkSmartPointer<vtkPolyData>::New(); 
	pointsPolydata->SetPoints(pointsVertex);

	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =    vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexFilter->SetInputData(pointsPolydata);
	vertexFilter->Update();

	vtkSmartPointer<vtkPolyData> polydatapc =    vtkSmartPointer<vtkPolyData>::New();
	polydatapc->ShallowCopy(vertexFilter->GetOutput());
	polydatapc->GetPointData()->SetScalars(colorsp);
	
	//restart vertices
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) 
		vi->info().m_index=-1;
	
	mapper1->SetInputData(unstructuredGrid2);
	mapper2->SetInputData(polydata);
	mapper3->SetInputData(polydatapc);
	
	// for generate vtk file
	/*vtkSmartPointer<vtkUnstructuredGridWriter> writer =  vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writer->SetInputData(unstructuredGrid2);
	writer->SetFileName("carp/original/tauro.vtk");
	writer->Write();*/

}

void aroundVertex(Delaunay & dt,Vertex & v,Vertex & w,int label){
	actor4->GetProperty()->SetOpacity(1);

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int k = 0;
	points->InsertNextPoint(v->point().x(), v->point().y(), v->point().z());
	v->info().m_index=k;++k;

	list<Vertex> vertices;
	dt.finite_incident_vertices(v,back_inserter(vertices));
	for(std::list<Vertex>::iterator iv = vertices.begin(); iv != vertices.end(); iv++){
		Vertex vi = *iv;
		points->InsertNextPoint(vi->point().x(), vi->point().y(), vi->point().z());
		vi->info().m_index=k;
		++k;
	}

	if(w->info().m_index == -1){
		points->InsertNextPoint(w->point().x(), w->point().y(), w->point().z());
		w->info().m_index=k;++k;
	}

	list<Vertex> verticesW;
	dt.finite_incident_vertices(w,back_inserter(verticesW));
	for(std::list<Vertex>::iterator iv = verticesW.begin(); iv != verticesW.end(); iv++){
		Vertex vi = *iv;
		if(vi->info().m_index == -1){
			points->InsertNextPoint(vi->point().x(), vi->point().y(), vi->point().z());
			vi->info().m_index=k;
			++k;
		}
	}
	
	//cells
	vtkSmartPointer<vtkUnsignedCharArray> colors =  vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid2 =  vtkSmartPointer<vtkUnstructuredGrid>::New();
	unstructuredGrid2->SetPoints(points);
	vtkSmartPointer<vtkCellArray> tetrahedrons = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();

	list<Cell> cells;
	dt.finite_incident_cells(v,back_inserter(cells));
	for(std::list<Cell>::iterator  ite = cells.begin(); ite != cells.end(); ite++){ 
		Cell it = *ite;
		if(it->info().m_label == label ){
			colors->InsertNextTypedTuple( *(colorsArray + it->info().m_label) );
			tetra->GetPointIds()->SetId(0, it->vertex(0)->info().m_index);
			tetra->GetPointIds()->SetId(1, it->vertex(1)->info().m_index);
			tetra->GetPointIds()->SetId(2, it->vertex(2)->info().m_index);
			tetra->GetPointIds()->SetId(3, it->vertex(3)->info().m_index);
			tetrahedrons->InsertNextCell(tetra);
		}
	}

	list<Cell> cellsW;
	dt.finite_incident_cells(w,back_inserter(cellsW));
	for(std::list<Cell>::iterator  ite = cellsW.begin(); ite != cellsW.end(); ite++){ 
		Cell it = *ite;
		if(it->info().m_label == label ){
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

	vtkSmartPointer<vtkShrinkFilter> shrink =vtkSmartPointer<vtkShrinkFilter>::New();
	shrink->SetInputData( unstructuredGrid2 );
	shrink->SetShrinkFactor( 1 );
	
	//restart index
	for(std::list<Vertex>::iterator iv = vertices.begin(); iv != vertices.end(); iv++)
		(*iv)->info().m_index=-1;
	for(std::list<Vertex>::iterator iv = verticesW.begin(); iv != verticesW.end(); iv++)
		(*iv)->info().m_index=-1;

	//points
	vtkSmartPointer<vtkPoints> pointsVertex =   vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");
	
	vector<Point> vpoints;
	vpoints.push_back(v->point());
	vpoints.push_back(w->point());
	for(unsigned int i=0;i<vpoints.size();++i){
		pointsVertex->InsertNextPoint (vpoints[i].x(), vpoints[i].y(), vpoints[i].z());
		colorsp->InsertNextTypedTuple( *(colorsArray + 8 )  );
	}
	vtkSmartPointer<vtkPolyData> pointsPolydata =   vtkSmartPointer<vtkPolyData>::New(); 
	pointsPolydata->SetPoints(pointsVertex);

	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =    vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexFilter->SetInputData(pointsPolydata);
	vertexFilter->Update();

	vtkSmartPointer<vtkPolyData> polydatapc =    vtkSmartPointer<vtkPolyData>::New();
	polydatapc->ShallowCopy(vertexFilter->GetOutput());
	polydatapc->GetPointData()->SetScalars(colorsp);
	//end points

	//mapper->SetInputData(unstructuredGrid2);	
	mapper3->SetInputData(polydatapc);
	mapper4->SetInputConnection(shrink->GetOutputPort());
}

//for detect bad oriented tetras
void aroundVertex2(Delaunay & dt,Vertex & v,Vertex & w,int label){
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	list<Vertex> vertices;int k = 0;
	points->InsertNextPoint(v->point().x(), v->point().y(), v->point().z());
	v->info().m_index=k;++k;

	dt.finite_incident_vertices(v,back_inserter(vertices));
	for(std::list<Vertex>::iterator iv = vertices.begin(); iv != vertices.end(); iv++){
		Vertex vi = *iv;
		points->InsertNextPoint(vi->point().x(), vi->point().y(), vi->point().z());
		vi->info().m_index=k;
		++k;
	}

	list<Vertex> verticesW;
	if(w->info().m_index == -1){
		points->InsertNextPoint(w->point().x(), w->point().y(), w->point().z());
		w->info().m_index=k;++k;
	}

	dt.finite_incident_vertices(w,back_inserter(verticesW));
	for(std::list<Vertex>::iterator iv = verticesW.begin(); iv != verticesW.end(); iv++){
		Vertex vi = *iv;
		if(vi->info().m_index == -1){
			points->InsertNextPoint(vi->point().x(), vi->point().y(), vi->point().z());
			vi->info().m_index=k;
			++k;
		}
	}
	
	//cells
	vtkSmartPointer<vtkUnsignedCharArray> colors =  vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid2 =    vtkSmartPointer<vtkUnstructuredGrid>::New();
	unstructuredGrid2->SetPoints(points);
	vtkSmartPointer<vtkCellArray> tetrahedrons = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();

	list<Cell> cells;
	dt.finite_incident_cells(v,back_inserter(cells));
	int count=0;
	for(std::list<Cell>::iterator  ite = cells.begin(); ite != cells.end(); ite++){ 
		Cell it = *ite;
		if(it->info().m_label > -1 ){
			if(it->info().m_testCell == 1)
				colors->InsertNextTypedTuple( *(colorsArray + 8) );
			else
				colors->InsertNextTypedTuple( *(colorsArray + it->info().m_label) );

			tetra->GetPointIds()->SetId(0, it->vertex(0)->info().m_index);
			tetra->GetPointIds()->SetId(1, it->vertex(1)->info().m_index);
			tetra->GetPointIds()->SetId(2, it->vertex(2)->info().m_index);
			tetra->GetPointIds()->SetId(3, it->vertex(3)->info().m_index);
			tetrahedrons->InsertNextCell(tetra);
		}
	}

	list<Cell> cellsW;
	dt.finite_incident_cells(w,back_inserter(cellsW));
	for(std::list<Cell>::iterator  ite = cellsW.begin(); ite != cellsW.end(); ite++){ 
		Cell it = *ite;
		if(it->info().m_label > -1 ){
			if(it->info().m_testCell == 1)
				colors->InsertNextTypedTuple( *(colorsArray + 8) );
			else
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

	vtkSmartPointer<vtkShrinkFilter> shrink =vtkSmartPointer<vtkShrinkFilter>::New();
	shrink->SetInputData( unstructuredGrid2 );
	shrink->SetShrinkFactor( 1 );
	
	//restart index
	for(std::list<Vertex>::iterator iv = vertices.begin(); iv != vertices.end(); iv++){
		Vertex vi = *iv;
		vi->info().m_index=-1;
	}
	for(std::list<Vertex>::iterator iv = verticesW.begin(); iv != verticesW.end(); iv++){
		Vertex vi = *iv;
		vi->info().m_index=-1;
	}

	//points
	vtkSmartPointer<vtkPoints> pointsVertex =   vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");
	
	vector<Point> vpoints;
	vpoints.push_back(v->point());
	vpoints.push_back(w->point());

	/*for(std::list<Vertex>::iterator iv = vertices.begin(); iv != vertices.end(); iv++){
		Vertex vi = *iv;
		vpoints.push_back(Point(vi->point().x(), vi->point().y(), vi->point().z()));
	}*/

	for(unsigned int i=0;i<vpoints.size();++i){
		pointsVertex->InsertNextPoint (vpoints[i].x(), vpoints[i].y(), vpoints[i].z());
		colorsp->InsertNextTypedTuple( *(colorsArray + 6 )  );
	}


	vtkSmartPointer<vtkPolyData> pointsPolydata =   vtkSmartPointer<vtkPolyData>::New(); 
	pointsPolydata->SetPoints(pointsVertex);

	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =    vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexFilter->SetInputData(pointsPolydata);
	vertexFilter->Update();

	vtkSmartPointer<vtkPolyData> polydatapc =    vtkSmartPointer<vtkPolyData>::New();
	polydatapc->ShallowCopy(vertexFilter->GetOutput());
	polydatapc->GetPointData()->SetScalars(colorsp);
	
	//mapper->SetInputData(unstructuredGrid2);
	//mapper1->SetInputConnection(shrink->GetOutputPort());
	mapper3->SetInputData(polydatapc);
	mapper4->SetInputData(unstructuredGrid2);
}

//draw 1-neighbord of vertex and points
void vertexRegion(Delaunay & dt,Vertex & v,vector<Point> &vertexPoints,int label){
	actor4->GetProperty()->SetOpacity(0.4);

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int k = 0;
	points->InsertNextPoint(v->point().x(), v->point().y(), v->point().z());
	v->info().m_index=k;++k;

	list<Vertex> vertices;
	dt.finite_incident_vertices(v,back_inserter(vertices));
	for(std::list<Vertex>::iterator iv = vertices.begin(); iv != vertices.end(); iv++){
		Vertex vi = *iv;
		points->InsertNextPoint(vi->point().x(), vi->point().y(), vi->point().z());
		vi->info().m_index=k;
		++k;
	}

	vtkSmartPointer<vtkUnsignedCharArray> colors =  vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid2 =  vtkSmartPointer<vtkUnstructuredGrid>::New();
	unstructuredGrid2->SetPoints(points);
	vtkSmartPointer<vtkCellArray> tetrahedrons = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();

	list<Cell> cells;
	dt.finite_incident_cells(v,back_inserter(cells));
	for(std::list<Cell>::iterator  ite = cells.begin(); ite != cells.end(); ite++){ 
		Cell it = *ite;
		/*if(it->info().m_label == label ){
			colors->InsertNextTypedTuple( *(colorsArray + it->info().m_label) );
			tetra->GetPointIds()->SetId(0, it->vertex(0)->info().m_index);
			tetra->GetPointIds()->SetId(1, it->vertex(1)->info().m_index);
			tetra->GetPointIds()->SetId(2, it->vertex(2)->info().m_index);
			tetra->GetPointIds()->SetId(3, it->vertex(3)->info().m_index);
			tetrahedrons->InsertNextCell(tetra);
		}*/
		colors->InsertNextTypedTuple( *(colorsArray + it->info().m_label) );
		tetra->GetPointIds()->SetId(0, it->vertex(0)->info().m_index);
		tetra->GetPointIds()->SetId(1, it->vertex(1)->info().m_index);
		tetra->GetPointIds()->SetId(2, it->vertex(2)->info().m_index);
		tetra->GetPointIds()->SetId(3, it->vertex(3)->info().m_index);
		tetrahedrons->InsertNextCell(tetra);
	}
	
	unstructuredGrid2->SetCells(VTK_TETRA, tetrahedrons);
	unstructuredGrid2->GetCellData()->SetScalars(colors);

	vtkSmartPointer<vtkShrinkFilter> shrink =vtkSmartPointer<vtkShrinkFilter>::New();
	shrink->SetInputData( unstructuredGrid2 );
	shrink->SetShrinkFactor( 1 );
	
	//restart index
	for(std::list<Vertex>::iterator iv = vertices.begin(); iv != vertices.end(); iv++){
		Vertex vi = *iv;
		vi->info().m_index=-1;
	}

	vtkSmartPointer<vtkPoints> pointsVertex =   vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");
		
	/*pointsVertex->InsertNextPoint (v->point().x(),v->point().y(), v->point().z());
	colorsp->InsertNextTypedTuple( *(colorsArray)  );*/

	for(unsigned int i=0;i<vertexPoints.size();++i){
		pointsVertex->InsertNextPoint (vertexPoints[i].x(), vertexPoints[i].y(), vertexPoints[i].z());
		colorsp->InsertNextTypedTuple( *(colorsArray + i%10 )  );
	}

	vtkSmartPointer<vtkPolyData> pointsPolydata =   vtkSmartPointer<vtkPolyData>::New(); 
	pointsPolydata->SetPoints(pointsVertex);

	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =    vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexFilter->SetInputData(pointsPolydata);
	vertexFilter->Update();

	vtkSmartPointer<vtkPolyData> polydatapc =    vtkSmartPointer<vtkPolyData>::New();
	polydatapc->ShallowCopy(vertexFilter->GetOutput());
	polydatapc->GetPointData()->SetScalars(colorsp);	
	
	mapper1->SetInputConnection(shrink->GetOutputPort());
	mapper3->SetInputData(polydatapc);
	//mapper4->SetInputData(unstructuredGrid2);

	//mapper4->SetInputConnection(shrink->GetOutputPort());
}

void two_ring(Delaunay & dt,Vertex & v,vector<Point> &vertexPoints,int label){
	//points
	vector<Vertex> restartIndex;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	list<Vertex> vertices;int k = 0;
	points->InsertNextPoint(v->point().x(), v->point().y(), v->point().z());
	v->info().m_index=k;++k;
	restartIndex.push_back(v);

	vector<Vertex> vertices2ring;
	dt.finite_incident_vertices(v,back_inserter(vertices));
	for(std::list<Vertex>::iterator iv = vertices.begin(); iv != vertices.end(); iv++){
		Vertex vi = *iv;
		points->InsertNextPoint(vi->point().x(), vi->point().y(), vi->point().z());
		vi->info().m_index = k;
		++k;
		vertices2ring.push_back(vi);
		restartIndex.push_back(vi);//restart
	}

	for(unsigned int i=0; i<vertices2ring.size();++i){
		list<Vertex> verticesRing;
		dt.finite_incident_vertices(vertices2ring[i],back_inserter(verticesRing));
		for(std::list<Vertex>::iterator iv = verticesRing.begin(); iv != verticesRing.end(); iv++){
			Vertex vi = *iv;
			if(vi->info().m_index == -1){
				points->InsertNextPoint(vi->point().x(), vi->point().y(), vi->point().z());
				vi->info().m_index=k;
				++k;
				restartIndex.push_back(vi);//restart
				//vertexPoints.push_back(vi->point());
			}
		}
	}

	//cells
	vtkSmartPointer<vtkUnsignedCharArray> colors =  vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid2 =    vtkSmartPointer<vtkUnstructuredGrid>::New();
	unstructuredGrid2->SetPoints(points);
	vtkSmartPointer<vtkCellArray> tetrahedrons = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
	
	list<Cell> cells;
	dt.finite_incident_cells(v,back_inserter(cells));
	for(std::list<Cell>::iterator ite = cells.begin(); ite != cells.end(); ite++){ 		
		Cell it = *ite;
		if(it->info().m_label == label ){
			colors->InsertNextTypedTuple( *(colorsArray + it->info().m_label) );
			tetra->GetPointIds()->SetId(0, it->vertex(0)->info().m_index);
			tetra->GetPointIds()->SetId(1, it->vertex(1)->info().m_index);
			tetra->GetPointIds()->SetId(2, it->vertex(2)->info().m_index);
			tetra->GetPointIds()->SetId(3, it->vertex(3)->info().m_index);
			tetrahedrons->InsertNextCell(tetra);			
		}
	}

	for(unsigned int i=0; i<vertices2ring.size();++i){
		list<Cell> cellsRing;
		dt.finite_incident_cells(vertices2ring[i],back_inserter(cellsRing));
		for(std::list<Cell>::iterator  ite = cellsRing.begin(); ite != cellsRing.end(); ite++){ 	
			Cell it = *ite;
			if(it->info().m_label == label ){
				colors->InsertNextTypedTuple( *(colorsArray + it->info().m_label) );
				tetra->GetPointIds()->SetId(0, it->vertex(0)->info().m_index);
				tetra->GetPointIds()->SetId(1, it->vertex(1)->info().m_index);
				tetra->GetPointIds()->SetId(2, it->vertex(2)->info().m_index);
				tetra->GetPointIds()->SetId(3, it->vertex(3)->info().m_index);
				tetrahedrons->InsertNextCell(tetra);
			}
		}
	}

	unstructuredGrid2->SetCells(VTK_TETRA, tetrahedrons);
	unstructuredGrid2->GetCellData()->SetScalars(colors);

	vtkSmartPointer<vtkShrinkFilter> shrink =vtkSmartPointer<vtkShrinkFilter>::New();
	shrink->SetInputData( unstructuredGrid2 );
	shrink->SetShrinkFactor( 1 );	
	
	//restart index	
	for(unsigned int i = 0 ; i < restartIndex.size(); ++i)
		restartIndex[i]->info().m_index=-1;

	//draw points
	vtkSmartPointer<vtkPoints> pointsVertex =   vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");
		
	vertexPoints.push_back(v->point());
	for(unsigned int i=0;i<vertexPoints.size();++i){
		pointsVertex->InsertNextPoint (vertexPoints[i].x(), vertexPoints[i].y(), vertexPoints[i].z());
		colorsp->InsertNextTypedTuple( *(colorsArray + i%10 )  );
	}

	vtkSmartPointer<vtkPolyData> pointsPolydata =   vtkSmartPointer<vtkPolyData>::New(); 
	pointsPolydata->SetPoints(pointsVertex);

	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =    vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexFilter->SetInputData(pointsPolydata);
	vertexFilter->Update();

	vtkSmartPointer<vtkPolyData> polydatapc =    vtkSmartPointer<vtkPolyData>::New();
	polydatapc->ShallowCopy(vertexFilter->GetOutput());
	polydatapc->GetPointData()->SetScalars(colorsp);
	//end draw points

	mapper3->SetInputData(polydatapc);
	mapper4->SetInputConnection(shrink->GetOutputPort());
}

void drawCells(Delaunay & dt, vector<Cell> & cells){	
	vector<Vertex> restartIndex;
	int k = 0;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	for(unsigned int i = 0; i < cells.size(); ++i){
		Cell c = cells[i];
		for( unsigned int j = 0; j < 4; ++j){
			Vertex v = c->vertex(j);
			if(v->info().m_index == -1){
				points->InsertNextPoint(v->point().x(), v->point().y(), v->point().z());
				v->info().m_index=k;
				++k;
				restartIndex.push_back(v);
			}
		}
	}
	vtkSmartPointer<vtkUnsignedCharArray> colors =  vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid2 =    vtkSmartPointer<vtkUnstructuredGrid>::New();
	unstructuredGrid2->SetPoints(points);
	vtkSmartPointer<vtkCellArray> tetrahedrons = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();

	for(unsigned int i = 0; i < cells.size(); ++i){
		Cell c = cells[i];
		if(c->info().m_label > -1 ){
			//colors->InsertNextTypedTuple( *(colorsArray + c->info().m_label + 5) );
			colors->InsertNextTypedTuple( *(colorsArray + 8) );
			tetra->GetPointIds()->SetId(0, c->vertex(0)->info().m_index);
			tetra->GetPointIds()->SetId(1, c->vertex(1)->info().m_index);
			tetra->GetPointIds()->SetId(2, c->vertex(2)->info().m_index);
			tetra->GetPointIds()->SetId(3, c->vertex(3)->info().m_index);
			tetrahedrons->InsertNextCell(tetra);
		}
	}
	unstructuredGrid2->SetCells(VTK_TETRA, tetrahedrons);
	unstructuredGrid2->GetCellData()->SetScalars(colors);

	vtkSmartPointer<vtkShrinkFilter> shrink =vtkSmartPointer<vtkShrinkFilter>::New();
	shrink->SetInputData( unstructuredGrid2 );
	shrink->SetShrinkFactor( 1 );	

	for(unsigned int i = 0 ; i < restartIndex.size(); ++i)
		restartIndex[i]->info().m_index=-1;

	mapper5->SetInputConnection(shrink->GetOutputPort());
}

//for drawing picture for thesis
void drawMesh_star(Delaunay & dt){
	vector<pair<Point,int>> centroids;
	vector<Facet> boundary;
	vector<pair<Point,Point>> edges;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int k=0;		
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
		/*cout<<f<<" "<<s<<endl;
		cout<<edges[i].first<<" "<<edges[i].second<<endl;*/
	}

	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =    vtkSmartPointer<vtkUnstructuredGrid>::New();
	unstructuredGrid->SetPoints(points);

	vtkSmartPointer<vtkUnsignedCharArray> colors =  vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	vtkSmartPointer<vtkCellArray> tetrahedrons = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();

	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++){
		if(it->info().m_label > -1 ){ 
			colors->InsertNextTypedTuple( *(colorsArray + it->info().m_label) );
			tetra->GetPointIds()->SetId(0, it->vertex(0)->info().m_index);
			tetra->GetPointIds()->SetId(1, it->vertex(1)->info().m_index);
			tetra->GetPointIds()->SetId(2, it->vertex(2)->info().m_index);
			tetra->GetPointIds()->SetId(3, it->vertex(3)->info().m_index);
			tetrahedrons->InsertNextCell(tetra);
		}
	}

	unstructuredGrid->SetCells(VTK_TETRA, tetrahedrons);
	unstructuredGrid->GetCellData()->SetScalars(colors);

	vtkSmartPointer<vtkShrinkFilter> shrink =vtkSmartPointer<vtkShrinkFilter>::New();
	shrink->SetInputData( unstructuredGrid );
	shrink->SetShrinkFactor( 0.7 );

	//draw all lines

	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid2 =    vtkSmartPointer<vtkUnstructuredGrid>::New();
	unstructuredGrid->SetPoints(points);

	vtkSmartPointer<vtkUnsignedCharArray> colorLines = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorLines->SetNumberOfComponents(3);
	colorLines->SetName("ColorsLines");

	vtkSmartPointer<vtkFloatArray> linescales = vtkSmartPointer<vtkFloatArray>::New();
    linescales->SetName("linescales");

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();

	for(Finite_edges_iterator ei=dt.finite_edges_begin();ei!=dt.finite_edges_end(); ei++){
		Vertex vs = ei->first->vertex(ei->second);
		Vertex vt = ei->first->vertex(ei->third);
		line->GetPointIds()->SetId(0,vs->info().m_index); 
		line->GetPointIds()->SetId(1,vt->info().m_index);
		lines->InsertNextCell(line);
		linescales->InsertNextValue(1);
		colorLines->InsertNextTypedTuple(*(colorsArray + 7));
	}

	//for(unsigned int i = 0; i < indexs.size(); ++i){
	//	line->GetPointIds()->SetId(0,indexs[i].first); 
	//	line->GetPointIds()->SetId(1,indexs[i].second);
	//	lines->InsertNextCell(line);
	//	linescales->InsertNextValue(10);
	//	colorLines->InsertNextTypedTuple(purple);
	//	//cout<<indexs[i].first<<" "<<indexs[i].second<<endl;
	//}

	unstructuredGrid2->SetCells(VTK_LINE, lines);
	unstructuredGrid2->GetCellData()->SetScalars(colorLines);

	vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
	linesPolyData->SetPoints(points);
	linesPolyData->SetLines(lines);	
	linesPolyData->GetCellData()->SetScalars(linescales);	
	
	//draw surface

	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid3 =    vtkSmartPointer<vtkUnstructuredGrid>::New();
	unstructuredGrid->SetPoints(points);

	vtkSmartPointer<vtkUnsignedCharArray> triangleColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	triangleColors->SetNumberOfComponents(3);
	triangleColors->SetName("Colors");	

	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
	for(unsigned int i=0;i<boundary.size();++i){
		triangle->GetPointIds()->SetId(0, boundary[i].first->vertex( (boundary[i].second+1)%4 )->info().m_index);
		triangle->GetPointIds()->SetId(1, boundary[i].first->vertex( (boundary[i].second+2)%4 )->info().m_index);
		triangle->GetPointIds()->SetId(2, boundary[i].first->vertex( (boundary[i].second+3)%4 )->info().m_index);
		triangles->InsertNextCell(triangle);
		triangleColors->InsertNextTypedTuple(red);
	}

	unstructuredGrid3->SetPoints(points);
	unstructuredGrid3->SetCells(VTK_TRIANGLE, triangles);
	unstructuredGrid3->GetCellData()->SetScalars(triangleColors);	
	//draw points 



	vtkSmartPointer<vtkPoints> pointsVertex =   vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");

	vtkSmartPointer<vtkFloatArray> scales = vtkSmartPointer<vtkFloatArray>::New();
    scales->SetName("scales");

	for(unsigned int i=0;i<centroids.size();++i){
		if( centroids[i].second > -1 ){
			pointsVertex->InsertNextPoint (centroids[i].first.x(), centroids[i].first.y(), centroids[i].first.z());
			//colorsp->InsertNextTypedTuple( *(colorsArray + centroids[i].second) );
			colorsp->InsertNextTypedTuple( grey );
			scales->InsertNextValue(0.001); 
		}
	}

	//pointsVertex->InsertNextPoint (25, 21, 30);
	////colorsp->InsertNextTypedTuple( *(colorsArray + centroids[i].second) );
	//colorsp->InsertNextTypedTuple( grey );
	//scales->InsertNextValue(0.001); 
	
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
	
	//-----------for writing

	/*vtkSmartPointer<vtkUnstructuredGridWriter> writer =  vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writer->SetInputData(unstructuredGrid);
	writer->SetFileName("multimaterial/region_2.vtk");
	writer->Write();	*/

	/*vtkSmartPointer<vtkPolyDataWriter> writer2 =    vtkSmartPointer<vtkPolyDataWriter>::New();
	writer2->SetInputData(linesPolyData);	
	writer2->SetFileName("multimaterial/edges_1.vtk");
	writer2->Write();*/
	
	/*vtkSmartPointer<vtkUnstructuredGridWriter> writer =  vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writer->SetInputData(unstructuredGrid3);
	writer->SetFileName("triangle_4.vtk");
	writer->Write();*/

	/*vtkSmartPointer<vtkPolyDataWriter> writer3 =    vtkSmartPointer<vtkPolyDataWriter>::New();
	writer3->SetInputConnection(glyph3D->GetOutputPort());
	writer3->SetFileName("multimaterial/s_region_1.vtk");
	writer3->Write();*/

	//// Create a mapper and actor
	//vtkSmartPointer<vtkDataSetMapper> mapper =   vtkSmartPointer<vtkDataSetMapper>::New();
	//vtkSmartPointer<vtkDataSetMapper> mapper2 =   vtkSmartPointer<vtkDataSetMapper>::New();
	//vtkSmartPointer<vtkDataSetMapper> mapper3 =   vtkSmartPointer<vtkDataSetMapper>::New();
	//vtkSmartPointer<vtkPolyDataMapper> mapper4 =    vtkSmartPointer<vtkPolyDataMapper>::New();
	//vtkSmartPointer<vtkDataSetMapper> mapper5 =   vtkSmartPointer<vtkDataSetMapper>::New();//triangles

	//mapper1->SetInputData(unstructuredGrid);

	mapper1->SetInputConnection(shrink->GetOutputPort());

	mapper5->SetInputData(linesPolyData);

	actor5->SetMapper(mapper5);
	actor5->GetProperty()->SetLineWidth(4);
 
}

void draw_circles(vector<pair<Point,double>> & circles){	
	//actor1->GetProperty()->SetOpacity(1);
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkFloatArray> scales = vtkSmartPointer<vtkFloatArray>::New();
    scales->SetName("scales");
    vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");
	//cout<<"drawing circles size: "<<circles.size()<<endl;
	for(unsigned int i = 0; i < circles.size(); ++i){	
		points->InsertNextPoint( circles[i].first.x() , circles[i].first.y() , circles[i].first.z()  ); // sphere in circle
        scales->InsertNextValue( 2*circles[i].second );
		//scales->InsertNextValue( circles[i].second );
		colorsp->InsertNextTypedTuple( *(colorsArray + 2) );
	}

	vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->SetPoints(points);
    grid->GetPointData()->AddArray(scales);    
    grid->GetPointData()->SetActiveScalars("scales"); // !!!to set radius first
    grid->GetPointData()->AddArray(colorsp);	

	vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
	polygonSource->GeneratePolygonOff();
	polygonSource->SetNumberOfSides(100);	

	vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    glyph3D->SetInputData(grid);
    glyph3D->SetSourceConnection(polygonSource->GetOutputPort());
	mapper5->SetInputConnection(glyph3D->GetOutputPort());
}

void draw_spheres(vector<pair<Point,double>> & circles){	
	actor5->GetProperty()->SetOpacity(0.3);
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkFloatArray> scales = vtkSmartPointer<vtkFloatArray>::New();
    scales->SetName("scales");
    vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	//cout<<"drawing circles size: "<<circles.size()<<endl;
	for(unsigned int i = 0; i < circles.size(); ++i){	
		points->InsertNextPoint( circles[i].first.x() , circles[i].first.y() , circles[i].first.z()  ); // sphere in circle
        scales->InsertNextValue( 2*circles[i].second );
		//scales->InsertNextValue( circles[i].second );
	}

	vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->SetPoints(points);
    grid->GetPointData()->AddArray(scales);    
    grid->GetPointData()->SetActiveScalars("scales"); // !!!to set radius first
    grid->GetPointData()->AddArray(colorsp);	

	vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
	polygonSource->GeneratePolygonOff();
	polygonSource->SetNumberOfSides(100);	


	vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	polydata->GetPointData()->AddArray(scales); 
	polydata->GetPointData()->SetActiveScalars("scales");

	vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
   // glyph3D->SetInputData(grid);
   // glyph3D->SetSourceConnection(polygonSource->GetOutputPort());
	
	glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
	glyph3D->SetInputData(polydata);		

	mapper5->SetInputConnection(glyph3D->GetOutputPort());
}