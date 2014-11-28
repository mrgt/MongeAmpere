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

#ifndef MA_KANTOROVICH_HPP
#define MA_KANTOROVICH_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

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
  template <class T, class Functions, class Matrix,
	    class Vector, class SparseMatrix>
  double kantorovich (const T &densityT,
		      const Functions &densityF,
		      const Matrix &X,
		      const Vector &weights,
		      Vector &g,
		      SparseMatrix &h)
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
    
    // compute the quadratic part
    typedef MA::Voronoi_intersection_traits<K> Traits;
    typedef typename MA::Tri_intersector<T,RT,Traits> Tri_isector;  
    typedef typename Tri_isector::Pgon Pgon;
    
    typedef Eigen::Triplet<FT> Triplet;
    std::vector<Triplet> htri;
    
    FT total(0), fval(0), total_area(0);
    g = Vector::Zero(N);

    MA::voronoi_triangulation_intersection_raw
      (densityT,dt,
       [&] (const Pgon &pgon,
	    typename T::Face_handle f,
	    Vertex_handle_RT v)
       {
	 Tri_isector isector;

	 Polygon p;
	 std::vector<Vertex_handle_RT> adj;
	 for (size_t i = 0; i < pgon.size(); ++i)
	   {
             size_t ii = (i==0)?(pgon.size()-1):(i-1);
	     //size_t ii = (i+1)%pgon.size();
	     p.push_back(isector.vertex_to_point(pgon[i], pgon[ii]));
	     adj.push_back((pgon[i].type == Tri_isector::EDGE_DT) ?
			   pgon[i].edge_dt.second : 0);
	   }

	 size_t idv = v->info();
	 auto fit = densityF.find(f);
	 assert(fit != densityF.end());
	 auto fv = fit->second; // function to integrate 
	 
	 // compute hessian
	 size_t num_adj = 0;
	 for (size_t i = 0; i < p.size(); ++i)
	   {
	     if (adj[i] == 0)
	       continue;
	     Vertex_handle_RT w = adj[i];
	     size_t idw = w->info();
	     
	     FT r = MA::integrate_1(p.edge(i), fv);
	     FT d = 2*sqrt(CGAL::squared_distance(v->point(),
						  w->point()));
	     htri.push_back(Triplet(idv, idw, -r/d));
	     htri.push_back(Triplet(idv, idv, +r/d));
	   }
	 
	 // compute value and gradient
	 FT warea = MA::integrate_1(p, fv);
	 FT intg = MA::integrate_3(p, [&](Point p) 
				   {
				     return fv(p) * 
				     CGAL::squared_distance(p,
							    v->point());
				   });
	 fval = fval + warea * weights[idv] - intg; 
	 g[idv] = g[idv] + warea;
	 total += warea;
         total_area += p.area();
       });
    h = SparseMatrix(N,N);
    h.setFromTriplets(htri.begin(), htri.end());
    h.makeCompressed();
    return fval;
  }
}

#endif

