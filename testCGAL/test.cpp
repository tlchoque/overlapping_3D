
#include "headers.h"
//#include "overlap.h"
const int windowWidth = 1000;//1000
const int windowHeight = 900;//900
const int windowDepth = 500;//900
const int numPoints = 2000;//1400-> remains 1 vertex or 2000 last 1200

int main() {
	//test_overlap();
	Delaunay dt;
	begin_VTK(dt);	

	generation(dt,numPoints,windowWidth,windowHeight,windowDepth);
	segmentation(dt,windowWidth,windowHeight,windowDepth);	

	std::stringstream ss,original,manifold,relabeling,smoothing,manifold_histo,relabeling_histo,smoothing_histo;
	string file = "hyena";

	ss<<"D:/project/samples/meshes3D/"<<file<<".txt";
	//ss<<"E:/thesis project/meshes/meshes3D/"<<"hyena_rep3"<<".txt";//important

	original<<file<<"/original";
	manifold<<file<<"/manifold";
	relabeling<<file<<"/relabeling";
	smoothing<<file<<"/smoothing";

	manifold_histo<<file<<"/manifold/histo.txt";
	relabeling_histo<<file<<"/relabeling/histo.txt";
	smoothing_histo<<file<<"/smoothing/histo.txt";

	/*ifstream input(ss.str(), std::ifstream::in);
	input>>(input,dt);*/

	
	/*obj_file(dt,( manifold.str() ).c_str(),"original" );
	getchar();*/
		
	double old = dt.number_of_vertices();

	//Delaunay old_dt = dt;
	//begin_overlapping(dt);	
	//create_star(dt);
	//save_vertex_star(dt);
	
	repairing(dt);	
	simulated_annealing_whithout_points(dt);

	drawMesh(dt);
	getchar();

	/*ofstream output("carp/manifold/carp_rep1.txt");
	output<<(output,dt); 
	output.close();*/

	off_file(dt,( manifold.str() ).c_str(),"post_test" );

	drawMesh(dt);
	getchar();


	//check_lalbel(dt);
	//saveMesh(dt,( manifold.str() ).c_str(),"rep2");
	
	/*get_inserted_points(dt,old,"hyena/manifold/rep2_new_points15.vtk",1.5);
	get_inserted_points(dt,old,"hyena/manifold/rep2_new_points2.vtk",2);
	get_inserted_points(dt,old,"hyena/manifold/rep2_new_points3.vtk",3);
	get_inserted_points(dt,old,"hyena/manifold/rep2_new_points4.vtk",4);
	get_inserted_points(dt,old,"hyena/manifold/rep2_new_points5.vtk",5);
	get_inserted_points(dt,old,"hyena/manifold/rep2_new_points7.vtk",7);*/

	//end_overlapping(old_dt,dt,old_num_vertices);


	/*mesh_curvature_histogram(dt,( manifold_histo.str() ).c_str() );
	saveMesh(dt,( manifold.str() ).c_str(),"manifold");*/

	/*ofstream output("hyena/manifold/hyena.txt");
	output<<(output,dt);
	output.close();
	drawMesh(dt);
	getchar();*/
				
	/*saveMesh(dt,( manifold.str() ).c_str(),"manifold");
	end_verification(dt,( manifold.str() ).c_str() );
	mesh_curvature_histogram(dt,( manifold_histo.str() ).c_str() );
	get_inserted_points(dt,old,"carp/manifold/new_points1.vtk",1);
	get_inserted_points(dt,old,"carp/manifold/new_points2.vtk",2);
	get_inserted_points(dt,old,"carp/manifold/new_points4.vtk",4);
	get_inserted_points(dt,old,"carp/manifold/new_points5.vtk",5);
	get_inserted_points(dt,old,"carp/manifold/new_points7.vtk",7);

	cout<<"finished repairing save"<<endl;*/

	/*begin_verification(dt);
	relabeling_to_smooth_2(dt);*/

	/*saveMesh(dt,( relabeling.str() ).c_str(),"relabeling");	
	end_verification(dt,( relabeling.str() ).c_str() );*/

	/*mesh_curvature_histogram(dt,( relabeling_histo.str() ).c_str() );
	cout<<"relabeling finished"<<endl;*/
	
	/*ofstream output("chest/manifold3/manifold.txt");
	output<<(output,dt);
	output.close();*/

	//relabeling_to_smooth(dt);
	//remove_boundary_picks(dt,120);

	/*saveMesh(dt,( relabeling.str() ).c_str(),"relabeling");
	mesh_curvature_histogram(dt,( relabeling_histo.str() ).c_str() );*/

	//saveMesh(dt,"hyena/manifold2");

	/*ofstream output("hyena/manifold/manifold.txt");
	output<<(output,dt);
	output.close();*/
	
	/*tsur(dt,5);
	saveMesh(dt,( smoothing.str() ).c_str(),"smoothing");
	mesh_curvature_histogram(dt,( smoothing_histo.str() ).c_str() );*/

	/*drawMesh(dt);
	getchar();
	CloseHandle(handler);
    Sleep(10);*/
	return 0;
}