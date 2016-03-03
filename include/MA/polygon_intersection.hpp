// MongeAmpere++
// Copyright (C) 2014 Quentin Mérigot, CNRS
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#ifndef MA_POLYGON_INTERSECTION_HPP
#define MA_POLYGON_INTERSECTION_HPP

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersection_2.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <MA/predicates.hpp>

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


  // Written by Joseph O'Rourke.
// Last modified: December 1997
// Questions to orourke@cs.smith.edu.
// --------------------------------------------------------------------
// This code is Copyright 1997 by Joseph O'Rourke.  It may be freely
// redistributed in its entirety provided that this copyright notice is
// not removed.
// --------------------------------------------------------------------

namespace detail
{
  enum tInFlag {Unknown, Qin, Pin};
  
  // Advances and prints out an inside vertex if appropriate.
  template<class K>
  size_t advance (size_t a, size_t *aa, size_t n,
		  bool inside,
		  typename CGAL::Polygon_2<K> &poly,
		  const typename K::Point_2 &v)
  {
    if (inside)
      poly.push_back(v);  
    (*aa)++;
    return  (a+1) % n;
  }
}

// P and Q _have_ to be in ccw order.
template<class K>
void polygon_polygon_intersection (const typename CGAL::Polygon_2<K> &P,
				   const typename CGAL::Polygon_2<K> &Q,
				   typename CGAL::Polygon_2<K> &r)
{
  using namespace detail;
  typedef typename K::Segment_2 Segment;
  typedef typename K::Point_2 Point;
  
  size_t n = P.size(), m = Q.size();
  tInFlag inflag = Unknown; // {Pin, Qin, Unknown}: which inside 
  bool FirstPoint = true;   
  size_t a = 0, b = 0;      // indices on P and Q (resp.)
  size_t aa = 0, ba = 0;    // # advances on a & b indices (after 1st inter.) 

  r = CGAL::Polygon_2<K>();
  if (n == 0 || m == 0)
    {
      return;
    }
  
  do
    {
      const size_t a1 = (a + n - 1) % n;
      const size_t b1 = (b + m - 1) % m;
      
      CGAL::Orientation
	cross = CGAL::orientation (P[a] - P[a1], Q[b] - Q[b1]),
	aHB = CGAL::orientation (Q[b1], Q[b], P[a]),
	bHA = CGAL::orientation (P[a1], P[a], Q[b]);
      
      // A = [P[a1],P[a]], B = [Q[b1], Q[b]]
      // If A & B intersect, update inflag. 
      CGAL::Object isect_object = CGAL::intersection(Segment(P[a1], P[a]),
						     Segment(Q[b1], Q[b]));
      Point p;
      if (CGAL::assign(p, isect_object))
	{
	  if (inflag == Unknown && FirstPoint)
	    {
	      FirstPoint = false;
	      aa = ba = 0;
	    }
	  r.push_back(p);
	  
	  // Update inflag. 
	  if (aHB == CGAL::POSITIVE)
	    inflag = Pin;
	  else if (bHA == CGAL::POSITIVE)
	    inflag = Qin;
	}
      
      //-----Advance rules-----
      // Special case: A & B overlap and oppositely oriented. */
      Segment s;
      if (CGAL::assign(s, isect_object))
	{
	  // the interior of the intersection is empty
	  if ( (P[a1] - P[a]) * (Q[b1] - Q[b]) < 0)
	    return;
	  
	  // FIXME: not sure about the other cases...
	}
      
      if ( (cross == CGAL::COLLINEAR) &&
	   (aHB == CGAL::NEGATIVE) &&
	   (bHA == CGAL::NEGATIVE) )
	{
	  // Special case: A & B parallel and separated -> P&Q disjoint. 
	  return;
	}
      else if ( (cross == CGAL::COLLINEAR) &&
		(aHB == CGAL::COLLINEAR) &&
		(bHA == CGAL::COLLINEAR) )
	{
	  // Special case: A & B collinear. 
	  // Advance but do not output point.
	  if (inflag == Pin)
	    b = advance<K>(b, &ba, m, inflag == Qin, r, Q[b]);
	  else
	    a = advance<K>(a, &aa, n, inflag == Pin, r, P[a]);
	}      
      else if ( (cross == CGAL::COLLINEAR) ||
		(cross == CGAL::POSITIVE))
	{
	  // Generic cases. 
	  if (bHA == CGAL::POSITIVE)
	    a = advance<K> (a, &aa, n, inflag == Pin, r, P[a]);
	  else
	    b = advance<K> (b, &ba, m, inflag == Qin, r, Q[b]);
	}
      else // if ( cross < 0 ) 
	{
	  if ( aHB == CGAL::POSITIVE)
	    b = advance<K> (b, &ba, m, inflag == Qin, r, Q[b]);
	  else
	    a = advance<K> (a, &aa, n, inflag == Pin, r, P[a]);
	}
    }
  while ( ((aa < n) || (ba < m)) && (aa < 2*n) && (ba < 2*m) );
  
  if (inflag == Unknown) 
    {
      if (P.has_on_bounded_side(Q[0]))
	r = Q;
      else if (Q.has_on_bounded_side(P[0]))
	r = P;
    }
  
  return;
}

}

#endif
