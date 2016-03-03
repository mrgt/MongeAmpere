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

#include <MA/common_rt.hpp>
#include <MA/voronoi_triangulation_intersection.hpp>
#include <MA/voronoi_polygon_intersection.hpp>
#include <MA/quadrature.hpp>
#include <MA/misc.hpp>
#include <MA/functions.hpp>
#include <Eigen/Dense>

namespace MA
{

  template <class T, class Functions, class Matrix, class Vector>
  void first_moment (const T &densityT,
		     const Functions &densityF,
		     const Matrix &X,
		     const Vector &weights,
		     Vector &masses,
		     Matrix &centroids)
  {
    using namespace details;
    typedef Eigen::Vector3d Vector3d;

    size_t N = X.rows();
    masses = Vector::Zero(N);
    centroids = Matrix::Zero(N,2);

    MA::voronoi_triangulation_intersection
      (densityT,
       details::make_regular_triangulation(X,weights),
       [&] (const details::Polygon &poly,
	    typename T::Face_handle f,
	    details::Vertex_handle_RT v)
       {
	 size_t idv = v->info();
	 auto fit = densityF.find(f);
	 assert(fit != densityF.end());
	 auto fv = fit->second; // function to integrate 

	 Vector3d intg = MA::integrate_3<Vector3d>(poly, Vector3d::Zero(),
						   [&](Point p) 
						   {
						     FT fp = fv(p);
						     return Vector3d(fp,         // area
								     fp * p.x(), // first moment
								     fp * p.y());
						   });
	 masses[idv] += intg[0];
	 centroids(idv,0) += intg[1];
	 centroids(idv,1) += intg[2];
       });
  }

  template <class T, class Functions, class Matrix, class Vector>
  void second_moment (const T &densityT,
		      const Functions &densityF,
		      const Matrix &X,
		      const Vector &weights,
		      Vector &masses,
		      Matrix &centroids,
		      Matrix &inertia)
  {
    using namespace details;
    typedef Eigen::Matrix<double, 6, 1> Vector6d;

    size_t N = X.rows();    
    masses = Vector::Zero(N);
    centroids = Matrix::Zero(N,2);
    inertia = Matrix::Zero(N,3);

    MA::voronoi_triangulation_intersection
      (densityT,
       details::make_regular_triangulation(X,weights),
       [&] (const details::Polygon &poly,
	    typename T::Face_handle f,
	    details::Vertex_handle_RT v)
       {
	 size_t idv = v->info();
	 auto fit = densityF.find(f);
	 assert(fit != densityF.end());
	 auto fv = fit->second; // function to integrate 

	 Vector6d intg = MA::integrate_3<Vector6d>(poly, Vector6d::Zero(),
						   [&](Point p) 
						   {
						     FT fp = fv(p);
						     Vector6d r;
						     // weighted area
						     r[0] = fp; 
						     // first moments
						     r[1] =  fp * p.x();  
						     r[2] = fp * p.y();
						     // second moments
						     r[3] = fp * p.x() * p.x();
						     r[4] = fp * p.y() * p.y(); 
						     r[5] = fp * p.x() * p.y();
						     return r;
						   });
	 masses[idv] += intg[0];
	 centroids(idv,0) += intg[1];
	 centroids(idv,1) += intg[2];
	 inertia(idv,0) += intg[3];
	 inertia(idv,1) += intg[4];
	 inertia(idv,2) += intg[5];
       });
  }


  template <class T, class Functions, class Matrix, class Vector>
  void lloyd (const T &densityT,
	      const Functions &densityF,
	      const Matrix &X,
	      const Vector &weights,
	      Vector &masses,
	      Matrix &centroids)
  {
    // compute first moments (integral of coordinates) and rescale
    // them so as to get centroids of Voronoi cells.
    first_moment(densityT, densityF, X, weights, masses, centroids);

    size_t N = X.rows();
    for (size_t i = 0; i < N; ++i)
      {
	centroids(i,0) /= masses[i];
	centroids(i,1) /= masses[i];
      }
  }
}

#endif

