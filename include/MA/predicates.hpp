#ifndef MA_PREDICATES_HPP
#define MA_PREDICATES_HPP

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
    auto E = details::dual_point(v,u2,u3);
    return side1(v,w,E);
  }
  
  template <class Point, class Segment>
  bool
  side3(const Point &v, const Point &u,
	const Segment &S, const Point &w)
  {
    auto L = details::dual_line(v,u);
    auto E = line_line_intersection(S.supporting_line(),L);
    return side1(v,w,E);
  }  
}

#endif

