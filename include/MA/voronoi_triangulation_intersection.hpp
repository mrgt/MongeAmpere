#ifndef MA_VORONOI_TRIANGULATION_INTERSECTION
#define MA_VORONOI_TRIANGULATION_INTERSECTION

#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>
#include <MA/predicates.hpp>
#include <queue>

namespace MA
{

  // Given two pairs f and g that contain a common element, the
  // function/common/ returns this element.
  template <class T>
  const T&
  common(const std::pair<T,T> &f, 
	 const std::pair<T,T> &g)
  {
    if (f.first == g.first)
      return f.first;
    if (f.first == g.second)
      return f.first;
    assert((f.second == g.first) || (f.second == g.second));
    return f.second;
  }

  template <class T, class DT, class Traits>
  struct Tri_intersector
  {
    typedef typename CGAL::Kernel_traits<typename DT::Point>::Kernel K;
    typedef typename CGAL::Line_2<K> Line;
    typedef typename CGAL::Point_2<K> Point;

    typedef typename T::Vertex_handle Vertex_handle_T;
    typedef typename DT::Vertex_handle Vertex_handle_DT;
    typedef std::pair<Vertex_handle_T, Vertex_handle_T> Edge_T;
    typedef std::pair<Vertex_handle_DT, Vertex_handle_DT> Edge_DT;

    enum EdgeType { EDGE_T, EDGE_DT };
    struct Pgon_edge
    {
      EdgeType type;
      Edge_T edge_t;
      Edge_DT edge_dt;

      bool
      operator == (const Pgon_edge &other) const
      {
	if (type != other.type)
	  return false;
	if (type == EDGE_T) 
	  return edge_t == other.edge_t;
	else if (type == EDGE_DT) 
	  return edge_dt == other.edge_dt;
      }
    };
    
    static Pgon_edge make_edge_t(Vertex_handle_T a, Vertex_handle_T b)
    {
      Pgon_edge e;
      e.type = EDGE_T;
      e.edge_t.first = a;
      e.edge_t.second = b;
      e.edge_dt.first = 0;
      e.edge_dt.second = 0;
      return e;
    }

    static Pgon_edge make_edge_dt(Vertex_handle_DT a, Vertex_handle_DT b)
    {
      Pgon_edge e;
      e.type = EDGE_DT;
      e.edge_dt.first = a;
      e.edge_dt.second = b;
      e.edge_t.first = 0;
      e.edge_t.second = 0;
      return e;
    }

    typedef std::pair<Pgon_edge,Pgon_edge> Pgon_vertex;

    struct Pgon_vertex_hash
    {
      bool operator() (const Pgon_edge &e)
      {
	std::size_t seed = 0; 
	boost::hash_combine(e.type, seed);
	boost::hash_combine(e.edge_dt.first, seed);
	boost::hash_combine(e.edge_dt.second, seed);
	boost::hash_combine(e.edge_t.first, seed);
	boost::hash_combine(e.edge_t.second, seed);
      }

      bool operator() (const Pgon_vertex &v)
      {
	std::size_t seed = 0;
	boost::hash_combine(seed, v.first);
	boost::hash_combine(seed, v.second);
	return seed;
      }
    };

    typedef std::vector<Pgon_edge> Pgon;

    Line
    edge_to_line (const Pgon_edge &e) const
    {
      if (e.type == EDGE_T)
	{
	  return Line(e.edge_t.first->point(),
		      e.edge_t.second->point());
	}
      // e.type == EDGE_DT
      typename Traits::Construct_dual construct_dual;
      return construct_dual(e.edge_dt.first->point(),
			    e.edge_dt.second->point());
    }
    
    Point
    vertex_to_point (const Pgon_vertex &E) const
    {
      return line_line_intersection(edge_to_line(E.first),
				    edge_to_line(E.second));
    }

    Point
    vertex_to_point (const Pgon_edge &a, const Pgon_edge &b) const
    {
      return vertex_to_point(Pgon_vertex(a,b));
    }

    // The function /inside/ determines whether the point p obtained
    // by intersecting the edges a and b is closer to v than to w, for
    // the weighed square distance i.e. it returns true if
    // ||p-v||^2 - w_v <= ||p-w||^2 - w_w and false if not.
    bool inside(Vertex_handle_DT v, Vertex_handle_DT w, 
		const Pgon_edge &a,
		const Pgon_edge &b) const
    {
      typename Traits::Side3 side3;
      if (a.type == EDGE_T && b.type == EDGE_T)
	{
	  typename Traits::Side1 side1;
	  Vertex_handle_T p = common(a.edge_t, b.edge_t);
	  return side1(v->point(), w->point(), p->point());
	}
      if (a.type == EDGE_DT && b.type == EDGE_DT)
	{
	  typename Traits::Side2 side2; 
	  return side2(v->point(),
		       b.edge_dt.second->point(),
		       a.edge_dt.second->point(),
		       w->point());
	}

      Vertex_handle_DT u;
      Vertex_handle_T p,q;
      if (a.type == EDGE_T && b.type == EDGE_DT)
	{
	  p = a.edge_t.first;
	  q = a.edge_t.second;
	  u = b.edge_dt.second;
      	}
      else
      	{
	  p = b.edge_t.first;
	  q = b.edge_t.second;
	  u = a.edge_dt.second;
      	}
      return side3(v->point(), u->point(),
		   p->point(), q->point(),
		   w->point());
    }

    void
    operator ()(const Pgon &P,
		Vertex_handle_DT v,
		Vertex_handle_DT w,
		Pgon &R) const
    {
      size_t n = P.size();
      if (n == 0)
	return;

      // edge corresponding to the bisector of [vw]
      Pgon_edge L = make_edge_dt(v,w);
      bool prev_inside = inside(v,w,P[n-1], P[0]);
      for (size_t i = 0; i < n; ++i)
	{
	  size_t ii = (i+1)%n;
	  bool cur_inside = inside(v,w,P[i], P[ii]);
	  if (prev_inside)
	    {
	      R.push_back(P[i]);
	      if (cur_inside == false) // do we cross ?
		R.push_back(L);
	    }
	   // if cur_inside == false, we are fully outside and we do
	   // not need to do anything. If cur_inside == true, we
	   // cross.
	  else if (cur_inside == true) 
	    R.push_back(P[i]);
	  prev_inside = cur_inside;
	}
    }
  };


  template <class K, class Tds>
  typename CGAL::Delaunay_triangulation_2<K,Tds>::Vertex_handle
  nearest_vertex(const typename CGAL::Delaunay_triangulation_2<K,Tds> &dt,
		 const typename CGAL::Point_2<K>  &p)
  {
    return dt.nearest_vertex(p);
  }

  template <class K, class Gt, class Tds>
  typename CGAL::Regular_triangulation_2<Gt,Tds>::Vertex_handle
  nearest_vertex(const typename CGAL::Regular_triangulation_2<Gt,Tds> &dt,
		 const typename CGAL::Point_2<K> &p)
  {
    return dt.nearest_power_vertex(p);
  }


  template <class T, class DT, class F>
  void
  voronoi_triangulation_intersection_raw(const T &t,
					 const DT &dt,
					 F out)
  {
    typedef typename CGAL::Kernel_traits<typename DT::Point>::Kernel K;
    typedef Voronoi_intersection_traits<K> Traits;
    typedef typename DT::Vertex_handle Vertex_handle_DT;
    typedef typename T::Face_handle Face_handle_T;
    typedef typename K::Point_2 Point;

    typedef Tri_intersector<T, DT, Traits> Tri_isector;
    typedef typename Tri_isector::Edge_DT Edge_DT;
    typedef typename Tri_isector::Pgon Pgon;
    typedef typename Tri_isector::Pgon_edge Pgon_edge;
    typedef typename Tri_isector::Pgon_vertex Pgon_vertex;
    typedef typename Tri_isector::Pgon_vertex_hash Pgon_vertex_hash;
    
    if (t.number_of_vertices() == 0)
      return;

    // insert seed
    Face_handle_T f = t.finite_faces_begin()++;
    Vertex_handle_DT v = nearest_vertex(dt, f->vertex(0)->point());

    typedef std::pair<Vertex_handle_DT, Face_handle_T> VF_pair;
    std::priority_queue<VF_pair> Q;
    std::set<VF_pair> visited;
    //std::map<Vertex_handle_DT, bool> vertex_visited;
    //for (auto v = dt.finite_vertices_begin();
    //     v != dt.finite_vertices_begin(); ++v)
    // vertex_visited[v] = false;

    Q.push(VF_pair(v,f));
    visited.insert(VF_pair(v,f));
    while (!Q.empty())
      {
	VF_pair vfp = Q.top(); Q.pop();
	Vertex_handle_DT v = vfp.first;
	//vertex_visited[v] = true;
	Face_handle_T f = vfp.second;

	Tri_isector isector;
	
	// convert triangle to Pgon
	Pgon R;
	R.push_back(isector.make_edge_t(f->vertex(0),f->vertex(1)));
	R.push_back(isector.make_edge_t(f->vertex(1),f->vertex(2)));
	R.push_back(isector.make_edge_t(f->vertex(2),f->vertex(0)));

	auto c = dt.incident_vertices (v), done(c);
	do
	  {
	    Pgon Rl;
	    if (dt.is_infinite(c))
	      continue;
	    isector(R, v, c, Rl);
	    R = Rl;
	  }
	while (++c != done);
	
	// propagate to neighbors
	for (auto e:R)
	  {
	    VF_pair p;
	    if (e.type == Tri_isector::EDGE_T)
	      {
		size_t i = f->index(e.edge_t.first);
		size_t j = f->index(e.edge_t.second);

		// with k = third vertex, 
		// i + j + k = 0 + 1 + 2 = 3 => k = 3 - i - j;
		size_t k = 3 - i - j;

		Face_handle_T fn = f->neighbor(k);
		if (t.is_infinite(fn))
		  continue;
		p = VF_pair(v, fn);
	      }
	    else // e.type == EDGE_DT
	      {
		p = VF_pair(e.edge_dt.second, f);
	      }

	    if (visited.find(p) != visited.end())
	      continue;
	    visited.insert(p);
	    Q.push(p);
	  }
	out(R, f, v);
      }

    // for (auto v = dt.finite_vertices_begin();
    // 	 v != dt.finite_vertices_end(); ++v)
    //   {
    // 	if (vertex_visited[v] == false)  
    // 	  std::cerr << "BUG => " << v->info() << "!!\n";
    //   }
  }

  template <class T, class DT, class F>
  void
  voronoi_triangulation_intersection(const T &t,
				     const DT &dt,
				     F out)
  {
    typedef typename CGAL::Kernel_traits<typename DT::Point>::Kernel K;
    typedef Voronoi_intersection_traits<K> Traits;
    typedef CGAL::Polygon_2<K> Polygon;
    typedef Tri_intersector<T, DT, Traits> Tri_isector;
    typedef typename Tri_isector::Pgon Pgon;

    voronoi_triangulation_intersection_raw
      (t, dt,
       [&] (const Pgon &R,
	    typename T::Face_handle f,
	    typename DT::Vertex_handle v)
       {
	 Tri_isector isector;
	 Polygon res;
	 size_t n = R.size();
	 for (size_t i = 0; i < n; ++i)
	   {
	     size_t ii = (i+1)%n;
	     res.push_back(isector.vertex_to_point(R[i], R[ii]));
	   }
	 out(res,f,v);
       });
  }  
}

#endif

