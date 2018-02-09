#ifndef COMMON_RT_HPP
#define COMMON_RT_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
//#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
//#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Polygon_2.h>

namespace MA
{
  namespace details
  {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Polygon_2<K> Polygon;
    typedef K::FT FT;
    //typedef CGAL::Regular_triangulation_filtered_traits_2<K> RT_Traits;
    typedef CGAL::Regular_triangulation_vertex_base_2<K> Vbase;
    typedef CGAL::Triangulation_vertex_base_with_info_2 <size_t, K, Vbase> Vb;
    typedef CGAL::Regular_triangulation_face_base_2<K> Cb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Cb> Tds;
    typedef CGAL::Regular_triangulation_2<K, Tds> RT;

    typedef RT::Vertex_handle Vertex_handle_RT;
    typedef RT::Weighted_point Weighted_point;
    typedef CGAL::Point_2<K> Point;
    
    // Helper function to insert points with indices into a regular triangulation
    template <class Matrix,
	      class Vector>
    RT
    make_regular_triangulation(const Matrix &X,
			       const Vector &weights)
    {
      size_t N = X.rows();
      assert(weights.rows() == N);
      assert(weights.cols() == 1);
      assert(X.cols() == 2);

      // insert points with indices in the regular triangulation
      std::vector<std::pair<Weighted_point,size_t> > Xw(N);
      for (size_t i = 0; i < N; ++i)
	{
	  Xw[i] = std::make_pair(Weighted_point(Point(X(i,0), X(i,1)),
						weights(i)), i);
	}

      RT dt (Xw.begin(), Xw.end());
      dt.infinite_vertex()->info() = -1;
      return dt;
    }
  }
}

#endif
