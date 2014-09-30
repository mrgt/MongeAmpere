#ifndef MA_VORONOI_INTERSECTION_HPP
#define MA_VORONOI_INTERSECTION_HPP

#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <MA/predicates.hpp>
#include <boost/variant.hpp>

namespace MA
{

  enum EdgeType { DELAUNAY, POLYGON };

  template <class Polygon, class DT, class Traits>
  struct Pgon_intersector
  {
    typedef typename CGAL::Kernel_traits<typename DT::Point>::Kernel K;
    typedef typename CGAL::Line_2<K> Line;
    typedef typename CGAL::Point_2<K> Point;
    typedef typename DT::Vertex_handle Vertex_handle;
    typedef boost::variant<Vertex_handle,size_t> Pgon_edge;
    typedef std::pair<Pgon_edge,Pgon_edge> Pgon_vertex;
    typedef std::vector<Pgon_vertex> Pgon;

    EdgeType edge_type (const Pgon_edge &e) const
    {
      if (e.type() == typeid(Vertex_handle))
	return DELAUNAY;
      return POLYGON;
    }

    Pgon_edge
    common_edge(const Pgon_vertex &f, 
		const Pgon_vertex &g) const
    {
      if (f.first == g.first)
	return f.first;
      if (f.first == g.second)
	return f.first;
      return f.second;
    }

    Line
    edge_to_line (Vertex_handle v, const Pgon_edge &e) const
    {
      typename Traits::Construct_dual construct_dual;
      if (const size_t *i = boost::get<size_t> (&e))
	return _polygon.edge(*i).supporting_line();
      const Vertex_handle *w = boost::get<Vertex_handle> (&e);
      assert(w);
      return construct_dual(v->point(), (*w)->point());
    }
    
    Point
    vertex_to_point (Vertex_handle v, const Pgon_vertex &E) const
    {
      Line l1 = edge_to_line(v,E.first), l2 = edge_to_line(v,E.second);
      return line_line_intersection(l1,l2);
    }

    bool inside(Vertex_handle v, Vertex_handle w, 
		const Pgon_vertex &E) const
    {
      EdgeType t1 = edge_type(E.first), t2 = edge_type(E.second);
      typename Traits::Side1 side1;
      typename Traits::Side2 side2; 
      typename Traits::Side3 side3;
      if (t1 == POLYGON && t2 == POLYGON)
	{
	  auto i2 = *boost::get<size_t> (&E.second);
	  return side1(v->point(), w->point(), _polygon[i2]);
	}
      if (t1 == DELAUNAY && t2 == DELAUNAY)
	{
	  auto u1 = *boost::get<Vertex_handle> (&E.first);
	  auto u2 = *boost::get<Vertex_handle> (&E.second);
	  return side2(v->point(), u2->point(),
		       u1->point(), w->point());
	}

      Vertex_handle u;
      size_t i;
      if (t1 == POLYGON && t2 == DELAUNAY)
	{
	  i = *boost::get<size_t> (&E.first);
	  u = *boost::get<Vertex_handle> (&E.second);
      	}
      else
      	{
      	  i = *boost::get<size_t> (&E.second);
      	  u = *boost::get<Vertex_handle> (&E.first);
      	}
      return side3(v->point(), u->point(),
		   _polygon.edge(i).source(),
		   _polygon.edge(i).target(),
		   w->point());
    }

    const Polygon &_polygon; 
    const DT &_dt;

    Pgon_intersector(const Polygon &polygon, const DT &dt):
      _polygon(polygon), _dt(dt)
    {}

    void
    operator ()(const Pgon &P,
		Vertex_handle v,
		Vertex_handle w,
		Pgon &R) const
    {
      size_t n = P.size();
      if (n == 0)
	return;
      
      Pgon_edge L (w); // edge corresponding to the bisector of [vw]
      auto S = P[n-1];
      R.clear();
      for (auto E:P)
	{
	  if (inside(v,w,E)) // E inside (negative side of) L
	    {
	      if (!inside(v,w,S)) // S not inside L
		{
		  auto edge = common_edge(S, E);
		  R.push_back(std::make_pair(L, edge));
		}
	      R.push_back(E);
	    }
	  else if (inside(v,w,S)) // S inside L
	    {
	      auto edge = common_edge(S, E);
	      R.push_back(std::make_pair(L, edge));
	    }
	  S = E;
	}
    }
  };

  template <class Polygon, class DT>
  Polygon 
  voronoi_polygon_intersection
                 (const Polygon &P,
		  const DT &dt,
		  const typename DT::Vertex_handle v)
  {
    typedef typename CGAL::Kernel_traits<typename DT::Point>::Kernel K;
    typedef Voronoi_intersection_traits<K> Traits;
    typedef Pgon_intersector<Polygon, DT, Traits> Pgon_isector;
    typedef typename Pgon_isector::Pgon Pgon;

    Pgon R;
    for (size_t i = 0; i < P.size(); ++i)
      R.push_back (std::make_pair(i,(i+1)%P.size()));

    Pgon_isector isector (P,dt);

    auto c = dt.incident_edges (v), done(c);
    do
      {
	Pgon Rl;
	if (dt.is_infinite(c))
	  continue;
	auto w = c->first->vertex(dt.ccw(c->second));
	isector(R, v, w, Rl);
	R = Rl;
      }
    while (++c != done);

    Polygon res;
    std::transform(R.begin(), R.end(), std::back_inserter(res),
		   [&](decltype(R[0]) p)
		   {return isector.vertex_to_point(v, p);});
    return res;
  }
}

#endif
