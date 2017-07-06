#include <r3d.h>
#include <r2d.h>
#include <v3d.h>
#include <v2d.h>
#include <iostream>
#include <vector>

#include <tests/test_helpers.h>
#include <tests/utest.h>
#include <math.h>

// numerical tolerances for pass/warn/fail tests
#define TOL_WARN 1.0e-8
#define TOL_FAIL 1.0e-4

// minimum volume allowed for test polyhedra 
#define MIN_VOL 1.0e-8

// number of trials per test
#define NUM_TRIALS 10

// for recursive tests, the stack capacity and 
// maximum recursion depth allowed
#define STACK_SIZE 64
#define MAX_DEPTH 16

// order of polynomial integration for all tests 
#define POLY_ORDER 2


// --                 unit tests                    -- //

void test_overlap() {

	// Intersects pairs of tetrahedra, for timing purposes only.

#undef NUM_TRIALS
#define NUM_TRIALS 1

	// variables: the polyhedra and their moments
	r3d_poly poly; 
	r3d_rvec3 verts[4];
	r3d_plane faces[4];
	//r3d_real om[R3D_NUM_MOMENTS(POLY_ORDER)];
	vector<r3d_real> om(R3D_NUM_MOMENTS(POLY_ORDER));

	// for timing only, no assertions
	r3d_int trial;
	printf("Intersecting %d pairs of tetrahedra.\n", NUM_TRIALS);
	for(trial = 0; trial < NUM_TRIALS; ++trial) {

		// generate the first random tet
		rand_tet_3d(verts, MIN_VOL);
		r3d_init_tet(&poly, verts);

		for(unsigned int i = 0 ; i < 4; ++i){
			cout<<verts[i].x<<" "<<verts[i].y<<" "<<verts[i].z<<endl;
		}
		//r3d_print(&poly);

		// generate the second random tet
		rand_tet_3d(verts, MIN_VOL);
		r3d_tet_faces_from_verts(faces, verts);

		for(unsigned int i = 0 ; i < 4; ++i){
			cout<<verts[i].x<<" "<<verts[i].y<<" "<<verts[i].z<<endl;
		}

		// clip the first tet against the faces of the second
		r3d_clip(&poly, faces, 4);

		// find the moments (up to quadratic order) of the clipped poly
		r3d_reduce(&poly, om, POLY_ORDER);

		for(unsigned int i = 0 ; i < om.size(); ++i){
			cout<<om[i]<<endl;;
		}
	}
	getchar();
}

r3d_rvec3 point2vec(Point &p){
	r3d_rvec3 tmp;
	tmp.x = p.x();
	tmp.y = p.y();
	tmp.z = p.z();
	return tmp;
}

double overlap_cells( Cell a, Cell b){
	r3d_poly poly; 
	r3d_rvec3 verts[4];
	r3d_plane faces[4];
	vector<r3d_real> om(R3D_NUM_MOMENTS(POLY_ORDER));
	for(unsigned int i = 0 ; i < 4; ++i)
		verts[i] = point2vec( a->vertex(i)->point() );
	//rand_tet_3d(verts, MIN_VOL);
	r3d_init_tet(&poly, verts);

	for(unsigned int i = 0 ; i < 4; ++i)
		verts[i] = point2vec( b->vertex(i)->point() );

	r3d_tet_faces_from_verts(faces, verts);
	r3d_clip(&poly, faces, 4);
	r3d_reduce(&poly, om, POLY_ORDER);
	
	/*for(unsigned int i = 0 ; i < om.size(); ++i){
			cout<<om[i]<<endl;;
	}*/
	return om[0];
}


void begin_overlapping(Delaunay &dt){
	for(Finite_cells_iterator it = dt.finite_cells_begin(); it != dt.finite_cells_end(); it++)
		it->info().m_previous_label = it->info().m_label;
}


vector<Cell> involved_cells(Delaunay &dt, Cell c, vector<Vertex> &vertices){
	vector<Vertex> involved_vertices;
	for(unsigned int i = 0 ; i < vertices.size(); ++i){
		Vertex v = vertices[i];
		if( dt.side_of_sphere(c,v->point(),false) == CGAL::ON_BOUNDED_SIDE ){
			involved_vertices.push_back(v);
		}
	}

	vector<Cell> involved_cells;
	for(unsigned int i = 0 ; i < involved_vertices.size(); ++i){
		Vertex w = involved_vertices[i];
		list<Cell> cells_w;
		dt.finite_incident_cells(w,back_inserter(cells_w));
		for(std::list<Cell>::iterator it = cells_w.begin(); it != cells_w.end(); it++){ 
			Cell cw = *it;
			if( cw->info().m_label == c->info().m_label && cw->info().m_overlap_state == false ){
				cw->info().m_overlap_state = true;
				involved_cells.push_back(cw);
			}
		}
	}

	for(unsigned int i = 0 ; i < involved_cells.size(); ++i)
		involved_cells[i]->info().m_overlap_state = false;

	return involved_cells;
}

void end_overlapping(Delaunay &dt_old, Delaunay &dt_new, double num_vertices){
	double overlapping=0,total=0;
	/*for(Finite_cells_iterator it = dt_new.finite_cells_begin(); it != dt_new.finite_cells_end(); it++){
		total += CGAL::volume(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point(), it->vertex(3)->point() );
		if( it->info().m_previous_label == it->info().m_label )
			overlapping += CGAL::volume(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point(), it->vertex(3)->point() );
	}*/	

	//get the new vertices
	vector<Vertex> new_vertices;
	int count = 0;
	for (Finite_vertices_iterator vi = dt_new.finite_vertices_begin(); vi != dt_new.finite_vertices_end(); vi++) { 
		++count;
		if( count > num_vertices){
			new_vertices.push_back(vi);
		}
	}

	//relate the old and new vertices
	for (Finite_vertices_iterator v_old = dt_old.finite_vertices_begin(),  v_new = dt_new.finite_vertices_begin(); 
		v_old != dt_old.finite_vertices_end(); v_old++, v_new ++) { 
		v_old->info().m_vertex_in_other_mesh = v_new;
	}

	//get the deleted cells
	vector<Cell> deleted_cells;
	for(Finite_cells_iterator it = dt_old.finite_cells_begin(); it != dt_old.finite_cells_end(); it++){
		total += CGAL::volume(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point(), it->vertex(3)->point() );
		Cell c_old = it;
		Cell c_new;
		Vertex v0 = c_old->vertex(0)->info().m_vertex_in_other_mesh;
		Vertex v1 = c_old->vertex(1)->info().m_vertex_in_other_mesh;
		Vertex v2 = c_old->vertex(2)->info().m_vertex_in_other_mesh;
		Vertex v3 = c_old->vertex(3)->info().m_vertex_in_other_mesh;
		if(!dt_new.is_cell(v0,v1,v2,v3,c_new) ){
			deleted_cells.push_back(c_old);
		}
		else{
			if( c_new->info().m_label == c_old->info().m_label ){
				overlapping += CGAL::volume(c_old->vertex(0)->point(), c_old->vertex(1)->point(), c_old->vertex(2)->point(), c_old->vertex(3)->point() );
			}
		}
	}

	//get the overlapping
	for(unsigned int i = 0 ; i < deleted_cells.size(); ++i){
		Cell del = deleted_cells[i];
		vector<Cell> involved = involved_cells(dt_new,del,new_vertices);
		for(unsigned int j = 0 ; j < involved.size();++ j ){
			Cell inv = involved[j];
			overlapping +=overlap_cells(del,inv);
		}
	}

	double difference = total - overlapping;
	cout<<"total area: "<<total<<endl;
	cout<<"partial area: "<<difference<<endl;
	cout<<"ratio: "<<difference/total<<endl;
}
