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

#ifndef MA_PREDICATES_HPP
#define MA_PREDICATES_HPP

//#include <CGAL/Regular_triangulation_euclidean_traits_2.h>

namespace MA
{
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

  class Voronoi_intersection_traits_base
  {
  public:    
    class Construct_dual
    {
    public:
      template <class K>
      typename CGAL::Line_2<K>
      operator()(const typename CGAL::Point_2<K> &p,
		 const typename CGAL::Point_2<K> &q) const
      {
	return CGAL::bisector(p,q);
      }
      
      template <class K>
      typename CGAL::Line_2<K>
      operator ()(const typename CGAL::Weighted_point_2<K> &p,
		  const typename CGAL::Weighted_point_2<K> &q) const
      {
	return CGAL::radical_axis(p,q);
      }

      template <class K>
      typename CGAL::Point_2<K>
      operator ()(const typename CGAL::Weighted_point_2<K> &p,
		  const typename CGAL::Weighted_point_2<K> &q,
		  const typename CGAL::Weighted_point_2<K> &r) const
      {
	return CGAL::weighted_circumcenter(p,q,r);
      }

      template <class K>
      typename CGAL::Point_2<K>
      operator()(const typename CGAL::Point_2<K> &p,
		 const typename CGAL::Point_2<K> &q,
		 const typename CGAL::Point_2<K> &r) const
      {
	return CGAL::circumcenter(p,q,r);
      }
    };
    
    class Side1
    {
    public:
      template <class K>
      bool
      operator()(const typename CGAL::Weighted_point_2<K> &p,
		 const typename CGAL::Weighted_point_2<K> &q,
		 const typename CGAL::Point_2<K> &E) const
      {
	typename K::Compare_power_distance_2 compare_power_distance =
	  K().compare_power_distance_2_object();
	//CGAL::Comparison_result c = compare_power_distance(p, q, E);
	CGAL::Comparison_result c = compare_power_distance(E, p, q);
	return (c == CGAL::SMALLER);
      }
      
      template <class K>
      bool
      operator()(const typename CGAL::Point_2<K> &p,
		 const typename CGAL::Point_2<K> &q,
		 const typename CGAL::Point_2<K> &E) const
      {
	CGAL::Comparison_result c = 
	  CGAL::compare_distance_to_point(E, p, q);
	return (c == CGAL::SMALLER);
      }
    };

    class Side2
    {
    public:
      typedef bool result_type;

      template <class Point>
      bool
      operator ()(const Point &v, const Point &u2, 
		  const Point &u3, const Point &w) const
      { 
	Construct_dual dual;
	Side1 side1;
	return side1(v,w,dual(v,u2,u3));
      }
    };

    class Side3
    {
    public:
      typedef bool result_type;

      template <class Point, class SegPoint>
      bool
      operator()(const Point &v, const Point &u,
		 const SegPoint &Sa, const SegPoint &Sb,
		 const Point &w) const
      {
	typedef typename CGAL::Kernel_traits<Point>::Kernel K;
	Construct_dual dual;
	Side1 side1;
	const CGAL::Line_2<K> L = dual(v,u);
	const CGAL::Line_2<K> M (Sa,Sb);
	return side1(v,w,line_line_intersection(M, L));
      }  
    };

  };

  template <class K>
  class Voronoi_intersection_traits
  { 
  public:
    typedef typename K::Exact_kernel_rt EK;
    typedef typename K::Approximate_kernel Approximate_kernel;

    
    typedef Voronoi_intersection_traits_base Exact_traits;
    typedef Voronoi_intersection_traits_base Filtering_traits;

    typedef typename K::C2E C2E;
    typedef typename K::C2F C2F;
    
    typedef Voronoi_intersection_traits_base::Construct_dual Construct_dual;

    // Side1 (= power_test_2 or side_of_oriented_circle_2) is already
    // filtered by CGAL
    typedef Voronoi_intersection_traits_base::Side1 Side1;

    typedef CGAL::Filtered_predicate<
      typename Exact_traits::Side2,
      typename Filtering_traits::Side2,
      C2E, C2F> Side2;

    typedef CGAL::Filtered_predicate<
      typename Exact_traits::Side3,
      typename Filtering_traits::Side3,
      C2E, C2F> Side3;
  };

}

#endif

