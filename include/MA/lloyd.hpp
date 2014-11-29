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

#ifndef MA_LLOYD_HPP
#define MA_LLOYD_HPP

#include <Eigen/Dense>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <MA/voronoi_triangulation_intersection.hpp>
#include <MA/voronoi_polygon_intersection.hpp>
#include <MA/quadrature.hpp>
#include <MA/misc.hpp>
#include <MA/functions.hpp>

namespace MA
{
  template <class T, class Functions, class Matrix, class Vector>
  void lloyd (const T &densityT,
	      const Functions &densityF,
	      const Matrix &X,
	      const Vector &weights,
	      Matrix &centroids,
	      Vector &masses)
  {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Polygon_2<K> Polygon;
    typedef K::FT FT;
    typedef CGAL::Regular_triangulation_filtered_traits_2<K> RT_Traits;
    typedef CGAL::Regular_triangulation_vertex_base_2<RT_Traits> Vbase;
    typedef CGAL::Triangulation_vertex_base_with_info_2
      <size_t, RT_Traits, Vbase> Vb;
    typedef CGAL::Regular_triangulation_face_base_2<RT_Traits> Cb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Cb> Tds;
    typedef CGAL::Regular_triangulation_2<RT_Traits, Tds> RT;

    typedef RT::Vertex_handle Vertex_handle_RT;
    typedef RT::Weighted_point Weighted_point;
    typedef typename CGAL::Point_2<K> Point;
    typedef typename CGAL::Vector_2<K> Vector_2;
    
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
    
    typedef MA::Voronoi_intersection_traits<K> Traits;
    typedef typename MA::Tri_intersector<T,RT,Traits> Tri_isector;  
    typedef typename Tri_isector::Pgon Pgon;

    masses = Vector::Zero(N);
    centroids = Matrix::Zero(N,2);

    MA::voronoi_triangulation_intersection
      (densityT,dt,
       [&] (const Polygon &poly,
	    typename T::Face_handle f,
	    Vertex_handle_RT v)
       {
	 size_t idv = v->info();
	 auto fit = densityF.find(f);
	 assert(fit != densityF.end());
	 auto fv = fit->second; // function to integrate 
	 
	 FT area = MA::integrate_1(poly, fv);
	 Vector_2 bary = MA::integrate_3(poly, [&](Point p) 
	 				 {
	 				   return (fv(p) * 
	 					   Vector_2(p-CGAL::ORIGIN));
	 				 });
	 masses[idv] = masses[idv] + area;
	 centroids(idv,0) += bary.x();
	 centroids(idv,1) += bary.y();
       });
    for (size_t i = 0; i < N; ++i)
      {
	centroids(i,0) /= masses[i];
	centroids(i,1) /= masses[i];
      }
  }
}

#endif

