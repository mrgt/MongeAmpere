#ifndef MA_OPTIMAL_TRANSPORT_HPP
#define MA_OPTIMAL_TRANSPORT_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <MA/voronoi_triangulation_intersection.hpp>
#include <MA/voronoi_polygon_intersection.hpp>
#include <MA/quadrature.hpp>
#include <MA/misc.hpp>
#include <MA/functions.hpp>
#include <MA/wolfe.hpp>

namespace MA
{
  template <class T, class Functions, class Matrix,
	    class Vector, class SparseMatrix>
  double ot_eval (const T &densityT,
		  Functions &densityF,
		  Matrix &X,
		  Vector &masses,
		  const Vector &weights,
		  Vector &g,
		  SparseMatrix &h)
  {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Polygon_2<K> Polygon;
    typedef K::FT FT;
    typedef CGAL::Regular_triangulation_filtered_traits_2<K> RT_Traits;
    typedef CGAL::Regular_triangulation_2<RT_Traits> RT;
    typedef RT::Vertex_handle Vertex_handle_RT;
    typedef RT::Weighted_point Weighted_point;
    typedef typename CGAL::Point_2<K> Point;
    
    size_t N = X.rows();
    assert(weights.rows() == N);
    assert(weight.cols() == 1);
    assert(X.cols() == 2);
    
    std::vector<Weighted_point> Xw(N);
    std::map<Point,size_t> indices;
    for (size_t i = 0; i < N; ++i)
    {
      Point p(X(i,0), X(i,1));
      indices[p] = i;
      Xw[i] = Weighted_point(p, weights(i));
    }
    RT dt (Xw.begin(), Xw.end());
    
    // compute the linear part of the function
    g = - masses;
    FT fval = - masses.dot(weights);
    
    // compute the quadratic part
    typedef MA::Voronoi_intersection_traits<K> Traits;
    typedef typename MA::Tri_intersector<T,RT,Traits> Tri_isector;  
    typedef typename Tri_isector::Pgon Pgon;
    
    typedef Eigen::Triplet<FT> Triplet;
    std::vector<Triplet> htri;
    
    FT total(0);
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
	     p.push_back(isector.vertex_to_point(pgon[i]));
	     typename Tri_isector::Pgon_edge e =
	       MA::common(pgon[i], pgon[(i+1)%pgon.size()]);
	     adj.push_back((e.type == Tri_isector::EDGE_DT) ?
			   e.edge_dt.second : 0);
	   }
	 
	 size_t idv = indices[v->point()];
	 
	 // compute hessian
	 for (size_t i = 0; i < p.size(); ++i)
	   {
	     if (adj[i] == 0)
	       continue;
	     Vertex_handle_RT w = adj[i];
	     size_t idw = indices[w->point()];
	     
	     FT r = MA::integrate_1(p.edge(i), densityF[f]);
	     FT d = 2*sqrt(CGAL::squared_distance(v->point(),
						  w->point()));
	     htri.push_back(Triplet(idv, idw, -r/d));
	     htri.push_back(Triplet(idv, idv, +r/d));
	   }
	 
	 // compute value and gradient
	 FT area = MA::integrate_1(p, densityF[f]);
	 FT intg = MA::integrate_3(p, [&](Point p) 
				   {
				     return densityF[f](p) * 
				     CGAL::squared_distance(p,
							    v->point());
				   });
	 fval = fval + area * weights[idv] - intg; 
	 g[idv] = g[idv] + area;
	 total += area;
       });
    h = SparseMatrix(N,N);
    h.setFromTriplets(htri.begin(), htri.end());
    h.makeCompressed();
    // std::cerr << "total = " << total  << "\n";
    //Eigen::VectorXd gg = fval.derivatives();
    //  std::cerr << "g[0] = " << g[0] << "/ " << gg(0) << "\n";
    return fval;
  }
}

#endif



