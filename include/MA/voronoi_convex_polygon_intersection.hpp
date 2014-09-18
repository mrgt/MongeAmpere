#ifndef MA_VORONOI_INTERSECTION_HPP
#define MA_VORONOI_INTERSECTION_HPP

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <MA/polygon_intersection.hpp>
#include <MA/simple_convex_hull_2.hpp>

namespace MA
{
  namespace details
  {
    template <class RT, class V>
    bool is_hidden(const RT &rt, const V &v)
    {
      return v->is_hidden();
    }
    
    template <class V, class K>
    bool is_hidden(const CGAL::Delaunay_triangulation_2<K> &dt,
		   const V &v)
    {
      return false;
    }
    
    template <class K>
    CGAL::Line_2<K>
    construct_dual_line(const CGAL::Point_2<K> &p,
			const CGAL::Point_2<K> &q)
    {
      return CGAL::bisector(p,q);
    }
    
    // template <class RT>
    // typename RT::Geom_traits::Line_2
    // construct_dual_line(const RT &dt, 
    // 			const typename RT::Weighted_point &p,
    // 			const typename RT::Weighted_point &q)
    // {
    //   return dt.geom_traits().construct_radical_axis_2_object()(p,q);
    // }
  }
  
  template <class K, class OutputIterator>
  void
  build_halfplanes_intersection
  (const std::vector<typename CGAL::Line_2<K>> &lines,
   const std::vector<size_t> &ch,
   OutputIterator out)
  {
    for (size_t i = 0; i < ch.size(); ++i)
      {
	size_t ii = (i+1)%ch.size();
	auto l = lines[ch[i]], ll = lines[ch[ii]];
	auto isect = CGAL::intersection(l, ll); 
	assert (isect);
	auto p = boost::get<typename CGAL::Point_2<K>>(&*isect);
	assert (p); // FIXME: could be line
	*out++ = *p;
      }
  }

  template <class K, class DT>
  void
  voronoi_convex_polygon_intersection
        (const std::vector<typename CGAL::Line_2<K>> &Plines,
	 const DT &dt,
	 const typename DT::Vertex_handle v,
	 typename CGAL::Polygon_2<K> &R)
  {
    typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
    typedef typename CGAL::Cartesian_converter<K,EK> K_to_EK;
    K_to_EK to_exact;

    if (details::is_hidden(dt, v))
      {
	R = typename CGAL::Polygon_2<K>();
	return;
      }
    auto c = dt.incident_edges (v), done(c);    
    std::vector<EK::Line_2> Elines;
    auto Alines = Plines;
    for(auto l:Plines)
      Elines.push_back(to_exact(l));

    auto Ap = v->point();
    auto Ep = to_exact(Ap);
    do
      {
	typename CGAL::Polygon_2<K> Rl;
	if (dt.is_infinite(c))
	  continue;
	auto w = c->first->vertex(dt.ccw(c->second));
	auto Aq = w->point();
	auto Eq = to_exact(Aq);
	Elines.push_back(details::construct_dual_line(Ep, Eq));
	Alines.push_back(details::construct_dual_line(Ap, Aq));
      }
    while (++c != done);

    std::vector<size_t> ch;
    halfplanes_intersection(Elines, Ep, ch);
    build_halfplanes_intersection(Alines, ch, std::back_inserter(R));
  }

}

#endif
