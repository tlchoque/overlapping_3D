
std::ostream & operator<< (std::ostream& os, Delaunay &dt)
{
  /*typedef Triangulation_3<GT, Tds>                 Triangulation;
  typedef typename Triangulation::size_type        size_type;
  typedef typename Triangulation::Vertex_handle    Vertex_handle;
  typedef typename Triangulation::Vertex_iterator  Vertex_iterator;
  typedef typename Triangulation::Cell_iterator    Cell_iterator;
  typedef typename Triangulation::Edge_iterator    Edge_iterator;
  typedef typename Triangulation::Facet_iterator   Facet_iterator;*/

  // outputs dimension and number of vertices
  Delaunay::size_type n = dt.number_of_vertices();
  if (CGAL::is_ascii(os))
    os << dt.dimension() << std::endl << n << std::endl;
  else{
	CGAL::write(os, dt.dimension());
    CGAL::write(os, n);
  }

  if (n == 0)
    return os;

  std::vector<Vertex> TV(n+1);
  Delaunay::size_type i = 0;

  // write the vertices

  for (Delaunay::Vertex_iterator it = dt.vertices_begin(), end = dt.vertices_end();
         it != end; ++it)
    TV[i++] = it;

  CGAL_triangulation_assertion( i == n+1 );
  CGAL_triangulation_assertion( dt.is_infinite(TV[0]) );

  CGAL::Unique_hash_map<Vertex, std::size_t > V;

  V[dt.infinite_vertex()] = 0;
  for (i=1; i <= n; i++) {
    os << *TV[i];
    V[TV[i]] = i;
    if (CGAL::is_ascii(os))
        os << std::endl;
  }

    // asks the tds for the combinatorial information
  dt.tds().print_cells(os, V);


  // write the non combinatorial information on the cells
  // using the << operator of Cell
  // works because the iterator of the tds traverses the cells in the
  // same order as the iterator of the triangulation
  switch ( dt.dimension() ) {
  case 3:
    {
      for(Delaunay::Cell_iterator it = dt.cells_begin(), end = dt.cells_end(); it != end; ++it) {
        os << *it; // other information
		os << (double)(it->info().m_color[0])/255<<" "<<(double)it->info().m_color[1]/255<<" "<<(double)it->info().m_color[2]/255<<" "<<it->info().m_label<<" ";
        if(CGAL::is_ascii(os))
          os << std::endl;
      }
      break;
    }
  case 2:
    {
      for(Delaunay::Facet_iterator it = dt.facets_begin(), end = dt.facets_end(); it != end; ++it) {
        os << *((*it).first); // other information
        if(CGAL::is_ascii(os))
          os << std::endl;
      }
      break;
    }
  case 1:
    {
      for(Delaunay::Edge_iterator it = dt.edges_begin(), end = dt.edges_end(); it != end; ++it) {
        os << *((*it).first); // other information
        if(CGAL::is_ascii(os))
          os << std::endl;
      }
      break;
    }
  }
  return os ;
}

std::istream &operator>> (std::istream& is, Delaunay &dt)
  // reads
  // the dimension
  // the number of finite vertices
  // the non combinatorial information on vertices (point, etc)
  // the number of cells
  // the cells by the indices of their vertices in the preceding list
  // of vertices, plus the non combinatorial information on each cell
  // the neighbors of each cell by their index in the preceding list of cells
  // when dimension < 3 : the same with faces of maximal dimension
{
  dt.tds().clear(); // infinite vertex deleted
  dt.set_infinite_vertex(dt.tds().create_vertex());

  std::size_t n;
  int d;
  if(CGAL::is_ascii(is))
     is >> d >> n;
  else {
    CGAL::read(is, d);
    CGAL::read(is, n);
  }

  if(!is) return is;
  dt.tds().set_dimension(d);

  //std::map< std::size_t, Vertex > V;
  std::vector< Vertex > V(n+1);
  V[0] = dt.infinite_vertex();

  // the infinite vertex is numbered 0

  for (std::size_t i = 1; i <= n; i++) {
    V[i] = dt.tds().create_vertex();
    if(!(is >> *V[i])) return is;
  }

  //std::map< std::size_t, Cell > C;
  std::vector< Cell > C;

  std::size_t m;
  dt.tds().read_cells(is, V, m, C);

  double r,g,b;int label;

  for (std::size_t j=0 ; j < m; j++){
	  is>>r>>g>>b>>label;
	C[j]->info().m_color[0] = (int)(255*r);
	C[j]->info().m_color[1] = (int)(255*g);
	C[j]->info().m_color[2] = (int)(255*b);
	//C[j]->info().m_label = label;
	
	//Point a,b,c,d,p;
	//a=C[j]->vertex(0)->point();
	//b=C[j]->vertex(1)->point();
	//c=C[j]->vertex(2)->point();
	//d=C[j]->vertex(3)->point();
	//p=Point(254,338,2);
	//if( !dt.is_infinite(C[j])   && (a== p || b== p || c==p|| d==p) ){
	//	C[j]->info().m_label = 0;
	//	C[j]->info().m_original_label = 0;
	//}
	//else{
	//C[j]->info().m_label = label;
	//C[j]->info().m_original_label = label;
	//}

	C[j]->info().m_label = label;
	C[j]->info().m_original_label = label;
    if(!(is >> *(C[j]))) return is;
  }

  std::cout<<"read finished"<<endl;

  //CGAL_triangulation_assertion( dt.is_valid(false) );
  return is;
}
