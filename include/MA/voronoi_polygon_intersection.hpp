#ifndef MA_VORONOI_INTERSECTION_HPP
#define MA_VORONOI_INTERSECTION_HPP

#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <MA/polygon_intersection.hpp>
#include <boost/variant.hpp>

namespace MA
{
  namespace details
  {
    template <class RT, class V>
    bool is_hidden(const RT &rt, const V &v)
    {
      return v->is_hidden();
    }
    
    template <class K, class V>
    bool is_hidden(const CGAL::Delaunay_triangulation_2<K> &dt, 
		   const V &v)
    {
      return false;
    }
    
    template <class K>
    typename CGAL::Line_2<K>
    dual_line(const typename CGAL::Point_2<K> &p,
	      const typename CGAL::Point_2<K> &q)
    {
      return CGAL::bisector(p,q);
    }
    
    template <class BP, class W>
    typename CGAL::Line_2
    <typename CGAL::Kernel_traits<BP>::Kernel>
    dual_line(const typename CGAL::Weighted_point<BP,W> &p,
	      const typename CGAL::Weighted_point<BP,W> &q)
    {
      return CGAL::radical_axis(p,q);
    }

    template <class BP, class W>
    typename CGAL::Point_2
    <typename CGAL::Kernel_traits<BP>::Kernel>
    dual_point(const typename CGAL::Weighted_point<BP,W> &p,
	       const typename CGAL::Weighted_point<BP,W> &q,
	       const typename CGAL::Weighted_point<BP,W> &r)
    {
      return CGAL::weighted_circumcenter(p,q,r);
    }

    template <class K>
    typename CGAL::Point_2<K>
    dual_point(const typename CGAL::Point_2<K> &p,
	       const typename CGAL::Point_2<K> &q,
	       const typename CGAL::Point_2<K> &r)
    {
      return CGAL::circumcenter(p,q,r);
    }
    
    template <class BP, class W>
    bool
    side1(const typename CGAL::Weighted_point<BP,W> &p,
	  const typename CGAL::Weighted_point<BP,W> &q,
	  const BP &E)
    {
      CGAL::Comparison_result c = CGAL::compare_power_distance(p, q, E);
      return (c == CGAL::SMALLER || c == CGAL::EQUAL);
    }

    template <class K>
    bool
    side1(const typename CGAL::Point_2<K> &p,
	  const typename CGAL::Point_2<K> &q,
	  const typename CGAL::Point_2<K> &E)
    {
      CGAL::Comparison_result c = CGAL::compare_distance(E, p, q);
      return (c == CGAL::SMALLER || c == CGAL::EQUAL);
    }

    template <class Point>
    bool
    side2(const Point &v, const Point &u2, 
	  const Point &u3, const Point &w)
    {
      auto E = dual_point(v,u2,u3);
      return side1(v,w,E);
    }

    template <class Point, class Segment>
    bool
    side3(const Point &v, const Point &u,
	  const Segment &S, const Point &w)
    {
      auto L = dual_line(v,u);
      auto E = line_line_intersection(S.supporting_line(),L);
      return side1(v,w,E);
    }

  }

  template <class Polygon, class DT>
  struct Pgon_intersector
  {
    typedef typename CGAL::Kernel_traits<typename DT::Point>::Kernel K;
    typedef typename CGAL::Line_2<K> Line;
    typedef typename CGAL::Point_2<K> Point;
    typedef typename DT::Vertex_handle Vertex_handle;
    typedef boost::variant<Vertex_handle,size_t> Pgon_edge;
    typedef std::pair<Pgon_edge,Pgon_edge> Pgon_vertex;
    typedef std::vector<Pgon_vertex> Pgon;

    enum EdgeType { DELAUNAY, POLYGON };

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
      if (const size_t *i = boost::get<size_t> (&e))
	return _polygon.edge(*i).supporting_line();
      const Vertex_handle *w = boost::get<Vertex_handle> (&e);
      assert(w);
      return details::dual_line(v->point(), (*w)->point());
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
      if (t1 == POLYGON && t2 == POLYGON)
	{
	  auto i2 = *boost::get<size_t> (&E.second);
	  return details::side1(v->point(),
				w->point(),
				_polygon[i2]);
	}
      if (t1 == DELAUNAY && t2 == DELAUNAY)
	{
	  auto u1 = *boost::get<Vertex_handle> (&E.first);
	  auto u2 = *boost::get<Vertex_handle> (&E.second);
	  return details::side2(v->point(),
				u2->point(),
				u1->point(), 
				w->point());
	}
      if (t1 == POLYGON && t2 == DELAUNAY)
	{
	  auto i = *boost::get<size_t> (&E.first);
	  auto u = *boost::get<Vertex_handle> (&E.second);
	  details::side3(v->point(),
			 u->point(),
			 _polygon.edge(i),
			 w->point());
	}
      Point p = vertex_to_point(v, E);
      return details::side1(v->point(), w->point(), p);
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
    if (details::is_hidden(dt, v))
      return Polygon();

    typedef Pgon_intersector<Polygon, DT> Pgon_isector;
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

  
#if 0
  template <class K, class Id>
  struct Pgon_with_id: 
    public std::vector<typename std::pair< typename CGAL::Point_2<K>,
     					   typename std::pair<Id,Id> > >
  {
  };

  template <class Id>
  Id
  find_common_id(const typename std::pair<Id,Id> &f, 
		 const typename std::pair<Id,Id> &g)
  {
    if (f.first == g.first)
      return f.first;
    if (f.first == g.second)
      return f.first;
    return f.second;
  }

  template <class K, class Id>
  void
  polygon_halfplane_intersection2(const Pgon_with_id<K,Id> &P,
				  const typename CGAL::Line_2<K> &L,
				  const Id &Lid,
				  Pgon_with_id<K,Id> &R)
  {
    size_t n = P.size();
    if (n == 0)
      return;

    auto S = P[n-1];
    R.clear();
    for (auto E:P)
      {
	if (inside(L,E.first)) // E inside (negative side of) L
	  {
	    if (!inside(L,S.first)) // S not inside L
	      {
		auto p = segment_line_intersection(S.first, E.first, L);
		auto segid = find_common_id(S.second, E.second);
		R.push_back(std::make_pair(p,
					   std::make_pair(segid, Lid)));
	      }
	    R.push_back(E);
	  }
	else if (inside(L,S.first)) // S inside L
	  {
	    auto p = segment_line_intersection(S.first, E.first, L);
	    auto segid = find_common_id(S.second, E.second);
	    R.push_back(std::make_pair(p, std::make_pair(segid, Lid)));
	  }
	S = E;
      }
  }
  
  template <class K, class DT>
  typename CGAL::Polygon_2<K> 
  voronoi_polygon_intersection
                 (const typename CGAL::Polygon_2<K> &P,
		  const DT &dt,
		  const typename DT::Vertex_handle v)
  {
    typename CGAL::Polygon_2<K> res;
    if (details::is_hidden(dt, v))
      return res;


    typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
    typedef typename CGAL::Cartesian_converter<K,EK> K_to_EK;
    K_to_EK to_EK;

    typedef typename CGAL::Line_2<K> Line;
    typedef Pgon_with_id<EK,Line> Polygon;

    Polygon R;
    size_t N = P.size();
    for (size_t i = 0; i < N; ++i)
      {
	size_t ip = (i+1) % N, im = (i+N-1) % N;
	R.push_back(std::make_pair(to_EK(P[i]),
				   std::make_pair(Line(P[im], P[i]),
						  Line(P[i], P[ip]))));
      }

    auto c = dt.incident_edges (v), done(c);
    do
      {
	Polygon Rl;
	if (dt.is_infinite(c))
	  continue;
	auto w = c->first->vertex(dt.ccw(c->second));
	auto El = details::dual_line(to_EK(w->point()),
				     to_EK(v->point()));
	auto l = details::dual_line(w->point(),
				    v->point());
	polygon_halfplane_intersection2(R, El, l, Rl);
	R = Rl;
      }
    while (++c != done);
	
    std::transform(R.begin(), R.end(), std::back_inserter(res),
		   [](decltype(R[0]) p)
		   {return line_line_intersection(p.second.first,
						  p.second.second);});
    return res;
  }
#endif

}

#endif
