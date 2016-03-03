#ifndef CGAL_TDS_INCREMENTAL_BUILDER_2_H
#define CGAL_TDS_INCREMENTAL_BUILDER_2_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/array.h>
#include <set>
#include <list>
#include <vector>

namespace CGAL {

template <class Triangulation_2> 
class Triangulation_incremental_builder_2
{
  typedef Triangulation_2 T;
  typedef typename T::Point Point;
  typedef typename T::Triangulation_data_structure Tds;
  typedef typename Tds::Vertex_handle Vertex_handle;
  typedef typename Tds::Face_handle Face_handle;
  typedef typename Tds::Edge Edge;
  typedef typename std::pair<Vertex_handle,Vertex_handle> Vh_pair;

public:
  Triangulation_incremental_builder_2(T &t) : _t(t) {}

  void begin_triangulation(int dim  = 2)
  {
    _t.tds().clear();
    _t.tds().set_dimension(dim); 
  }

  void end_triangulation()
  {
    Tds &tds = _t.tds();
    Vertex_handle vinf = 0;
    // deal with  boundaries
    if ( !edge_map.empty())
      {
	vinf = tds.create_vertex();
	std::map<Vh_pair, Edge> inf_edge_map;
	while (!edge_map.empty())
	  {
	    Face_handle fh = edge_map.begin()->second.first;
	    int ih = edge_map.begin()->second.second;
	    Face_handle fn = tds.create_face(vinf, 
					     fh->vertex(tds.cw(ih)), 
					     fh->vertex(tds.ccw(ih)));
	    vinf->set_face(fn);
	    tds.set_adjacency(fn, 0, fh, ih);
	    set_adjacency(fn, 1, inf_edge_map);
	    set_adjacency(fn, 2, inf_edge_map);
	    edge_map.erase(edge_map.begin());
	}
	CGAL_triangulation_assertion(inf_edge_map.empty());
      }
    reorient_faces();
    _t.set_infinite_vertex(vinf);
  }

  Vertex_handle add_vertex(const Point &p)
  {
    Vertex_handle v = _t.tds().create_vertex();
    v->set_point(p);
    return v;
  }

  Face_handle add_face(Vertex_handle vh0, Vertex_handle vh1,
		       Vertex_handle vh2)
  {
    Face_handle fh = _t.tds().create_face();
    fh->set_vertex(0, vh0);
    fh->set_vertex(1, vh1);
    fh->set_vertex(2, vh2);
    vh0->set_face(fh);
    vh1->set_face(fh);
    vh2->set_face(fh);
    for (std::size_t ih  = 0; ih < 3; ++ih)
      set_adjacency(fh, ih, edge_map);
    return fh;
  }

private:
  std::map<Vh_pair, Edge> edge_map;
  T &_t;

  void
  reorient_faces()
  {
    Tds &tds = _t.tds();
    typedef typename Tds::Face_iterator Face_iterator;
    // reorient the faces of a triangulation 
    // needed for example in off_file_input
    // because the genus is not known, the number of faces 
    std::set<Face_handle> oriented_set;
    std::stack<Face_handle>  st;
    Face_iterator fit = tds.faces_begin();
    int nf  = std::distance(tds.faces_begin(),tds.faces_end());
    
    while (static_cast<int>(oriented_set.size()) != nf)
      {
	while ( oriented_set.find(fit) != oriented_set.end())
	  {
	    ++fit; // find a germ for  non oriented components 
	  }
	// orient component
	oriented_set.insert(fit);
	st.push(fit);
	while ( ! st.empty()) {
	  Face_handle fh = st.top();
	  st.pop();
	  for(int ih = 0 ; ih < 3 ; ++ih){
	    Face_handle fn = fh->neighbor(ih);
	    if (oriented_set.find(fn) == oriented_set.end()){
	  int in = fn->index(fh);
	  if (fn->vertex(tds.cw(in)) != fh->vertex(tds.ccw(ih)))
	    fn->reorient();
	  oriented_set.insert(fn);
	  st.push(fn);
	    }
	  }
	}
	
      }
    return;
  }
  
  void
  set_adjacency(Face_handle fh, 
		int ih, 
		typename std::map<Vh_pair, Edge>& edge_map)
  {
    Tds &tds = _t.tds();
    // set adjacency to (fh,ih) using the the map edge_map
    // or insert (fh,ih) in edge map
    auto vhcw  =  fh->vertex(tds.cw(ih));
    auto vhccw =  fh->vertex(tds.ccw(ih)); 
    auto vhp =  vhcw < vhccw ?  
		       std::make_pair(vhcw, vhccw) 
		       : std::make_pair(vhccw, vhcw) ;
    auto emapit = edge_map.find(vhp);
    if (emapit == edge_map.end()) {// not found, insert edge
      edge_map.insert(std::make_pair(vhp, typename Tds::Edge(fh,ih)));
    }
    else
      { //found set adjacency and erase
	auto e = emapit->second;
	tds.set_adjacency( fh,ih, e.first, e.second);
	edge_map.erase(emapit);
      } 
  }
  
};


} //namespace CGAL

#endif // CGAL_TDS_INCREMENTAL_BUILDER_2_H //
