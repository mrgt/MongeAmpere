#ifndef MA_POLYGON_INTERSECTION_HPP
#define MA_POLYGON_INTERSECTION_HPP

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

namespace MA {

template <class K>
bool inside(const typename CGAL::Line_2<K> &L,
	    const CGAL::Point_2<K> &E)
{
  CGAL::Oriented_side o = L.oriented_side(E);
  return (o == CGAL::ON_NEGATIVE_SIDE ||
	  o == CGAL::ON_ORIENTED_BOUNDARY);
}

template <class K>
typename CGAL::Point_2<K>
segment_line_intersection (const typename CGAL::Point_2<K> &a, 
			   const typename CGAL::Point_2<K> &b, 
			   const typename CGAL::Line_2<K> &L)
{
  auto isect = CGAL::intersection(CGAL::Segment_2<K>(a,b),L);
  auto p = boost::get<typename CGAL::Point_2<K>>(&*isect);
  assert(p != 0);
  return *p;
}

template <class K>
typename CGAL::Point_2<K>
line_line_intersection (const typename CGAL::Line_2<K> &L,
			const typename CGAL::Line_2<K> &M)
{
  auto isect = CGAL::intersection(L,M);
  auto p = boost::get<typename CGAL::Point_2<K>>(&*isect);
  assert(p != 0);
  return *p;
}


template <class K>
void
polygon_halfplane_intersection(const typename CGAL::Polygon_2<K> &P,
			       const typename CGAL::Line_2<K> &L,
			       typename CGAL::Polygon_2<K> &R)
{
  size_t n = P.size();
  if (n == 0)
    return;
  R = typename CGAL::Polygon_2<K>();
  auto S = P[n-1];
  for (auto E = P.vertices_begin(); E != P.vertices_end(); ++E)
    {
      if (inside(L,*E)) // E inside (negative side of) L
	{
	  if (!inside(L,S)) // S not inside L
	    {
	      R.push_back(segment_line_intersection(S,*E,L));
	    }
	  R.push_back(*E);
	}
      else if (inside(L,S)) // S inside L
	{
	  R.push_back(segment_line_intersection(S,*E,L));
	}
      S = *E;
    }
}

}

#endif
