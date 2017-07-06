#define PI 3.14159265
#define HOLE 8

const int nroLabels = 9;

#include <iostream>
using namespace std;

#include "Imesh3D.h"// structure
#include "IOimesh3D.h"//  read and write files

void print(vector<Vertex> &vertices){
	for(unsigned int i = 0; i < vertices.size(); ++i){
		cout<<vertices[i]->point() <<"  ";
	}
	cout<<endl;
}
vector<Point> vtop(vector<Vertex> &vertices){
	vector<Point> points;
	for(unsigned int i = 0; i < vertices.size(); ++i){
		points.push_back(vertices[i]->point() );
	}
	return points;
}


#include "globalVTK.h"// draw mesh
#include "outVTK.h"//  read and write files
#include "additional.h"// operations in mesh

//for example mesh
#include "generation.h"
#include "segmentation.h"

//algorithms
#include "verify.h"
#include "operators.h"
#include "movement.h"
#include "singular_vertices.h"
#include "criteria.h"
#include "join.h" // join two components 

#include "Hausdorff.h"


#include "winding.h"

#include "point_insertion.h"
#include "repairing.h"
#include "repairing_2.h"
#include "repairing_mesh.h"

#include "simulated.h"
#include "simulated_2.h"

#include "test_curvature.h"
#include "histogram.h"

#include "pick_removal.h"
#include "optimal.h"

#include "relabeling.h"
#include "relabeling_2.h"

#include "tsur.h"
#include "smoothing.h"
#include "ODT.h"



