#include <iostream>
#include <vector>
#include <deque>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>// not
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

struct Ccell3D{
	int m_label;
	double m_minAngle;
	unsigned char m_color[3]; 
	bool m_auxiliar;// compute volume from facets and function getCells
	int m_state;//the cell has been relabeled = 1
	bool m_testCell;
	std::vector<double> m_dangles;
	bool m_visited;

	int m_original_label;

	int m_old_solution;
	int m_new_solution;

	int m_solution;
	
	bool m_hole_state;
	bool m_hole_visited;
	bool m_path_state;

	bool m_relabeled;

	//for overlapping
	int m_previous_label;
	bool m_overlap_state;
	//bool m_overlap_aux;


	bool visited;

	bool m_boundary;
	Ccell3D(){
		m_label= -1;// change to  -1
		m_original_label = -1;
		 
		m_color[0]=255;
		m_color[1]=255;
		m_color[2]=255;
		m_auxiliar = false;
		m_state=0;
		m_testCell=0;
		vector<double> d(6);
		m_dangles = d;

		m_visited = false;// this cell has not visited for singular test

		m_old_solution = -1;
		m_new_solution = -1;

		m_hole_state = false; //is not a hole
		m_hole_visited = false;
		m_path_state = false;

		m_relabeled = false;
		m_original_label = -1;


		m_previous_label = -1;
		m_overlap_state = false;
		//m_overlap_aux = false;

		visited = 0;

		m_boundary = 0;
	}
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Triangulation_cell_base_with_info_3<Ccell3D,K> Fb;

template < class GT, class Vb = CGAL::Triangulation_vertex_base_3<GT> >
class My_vertex_base
  : public Vb
{

public:
  typedef typename Vb::Vertex_handle  Vertex_handle;
  typedef typename Vb::Cell_handle    Cell_handle;
  typedef typename Vb::Point          Point;

  struct Info{
		int m_index;
		int m_state;//0:delaunay 1:boundary 2:singular
		bool m_auxiliar;
		bool m_boundary;// boundary vertex
		double m_tolerance;
		//for smoothing
		std::vector<vector<Vertex_handle>> m_s_stars;
		double m_s_greater_edge;
		double m_s_lower_angle;

		std::deque<std::deque<void*>> m_star;

		// updating
		int m_type;
		bool is_corner()		const { return m_type == 0; }
		bool is_internal()	const { return m_type == 1; }
		bool is_boundary()	const { return m_type == 2 || m_regions_around.size() > 1; }

		vector<int> m_regions_around;
		bool has_two_regions()		const { return m_regions_around.size() == 2; }
		bool has_more_than_two_regions()		const { return m_regions_around.size() > 2; }

		vector<int> m_singular_regions_around;
		bool is_singular()		const { return m_singular_regions_around.size() > 0; }

		Point m_smoothed_point;
		
		//update 2
		vector<pair<int,int>> m_singular_regions_and_subgraphs;
		bool is_singular_2()		const { return m_singular_regions_and_subgraphs.size() > 0; }		

		double m_curvature_sa;
		std::vector<vector<Vertex_handle>> m_sa_stars;// for simulated annealing
		std::vector<double> m_curvatures;
		double m_feature;

		Vertex_handle m_vertex_in_other_mesh;
		Info(){
			m_tolerance = INT_MAX;//initializa the vector star		
			m_index=-1;
			m_state=0;
			m_auxiliar = false;
			m_boundary= false;
			vector<vector<Vertex_handle>> stars(2);
			m_s_stars = stars;

			std::deque<std::deque<void*>> m_star2(nroLabels);
			m_star = m_star2;


			//updating
			m_type = -1;

			m_feature = false;
		}
	};

   Info _info;

  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef My_vertex_base<GT, Vb2>                        Other;
  };
  My_vertex_base() {}
  My_vertex_base(const Point& p)
    : Vb(p) {}
  My_vertex_base(const Point& p, Cell_handle c)
    : Vb(p, c) {}

  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
};

//typedef CGAL::Triangulation_vertex_base_with_info_3<Cvertex3D<K>, K>    Vb;
typedef CGAL::Triangulation_data_structure_3<My_vertex_base<K>,Fb>    Tds;
//typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds> Delaunay;

typedef K::Point_3								Point;
typedef K::Vector_3								Vector;
typedef Delaunay::Finite_cells_iterator			Finite_cells_iterator;
typedef Delaunay::Finite_edges_iterator			Finite_edges_iterator;
typedef Delaunay::Finite_facets_iterator		Finite_facets_iterator;
typedef Delaunay::Finite_vertices_iterator		Finite_vertices_iterator ;
typedef Delaunay::Vertex_handle					Vertex;
typedef Delaunay::Cell_handle					Cell;
typedef Delaunay::Face_circulator				Face_circulator;
typedef Delaunay::Facet_circulator				Facet_circulator;
typedef Delaunay::Cell_circulator				Cell_circulator;
typedef Delaunay::Facet							Facet;
typedef Delaunay::Edge							Edge;



