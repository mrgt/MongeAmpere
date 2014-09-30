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
#include <boost/timer/timer.hpp>
#include <lbfgs.hpp>
#include <cstdlib>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Vector_2<K> Vector;
typedef CGAL::Line_2<K> Line;
typedef CGAL::Polygon_2<K> Polygon;

//typedef CGAL::Simple_cartesian<AD> K_ad;

typedef CGAL::Regular_triangulation_filtered_traits_2<K> Traits;
typedef CGAL::Regular_triangulation_2<Traits> RT;
typedef RT::Vertex_handle Vertex_handle_RT;
typedef RT::Weighted_point Weighted_point;

typedef CGAL::Delaunay_triangulation_2<K> T;

double rr() 
{ 
  return 2*double(rand() / (RAND_MAX + 1.0))-1;
}

typedef Eigen::Triplet<double> Triplet;
typedef Eigen::SparseMatrix<double> SparseMatrix;
typedef Eigen::SparseVector<double> SparseVector;
typedef Eigen::VectorXd VectorXd;

template <class Functions>
double ot_eval (const T &densityT,
		Functions &densityF,
		const std::vector<Point> &X,
		VectorXd &masses,
		std::map<Point,size_t> indices,
		const VectorXd &weights,
		VectorXd &g)
{
  size_t N = X.size();
  std::vector<Weighted_point> Xw(N);
  for (size_t i = 0; i < N; ++i)
    Xw[i] = Weighted_point(X[i], weights[i]);
  RT dt (Xw.begin(), Xw.end());

  // compute the linear part of the function
  g = - masses;
  FT fval = - masses.dot(weights);

  // compute the quadratic part
  typedef MA::Voronoi_intersection_traits<K> Traits;
  typedef MA::Tri_intersector<T,RT,Traits> Tri_isector;  
  typedef typename Tri_isector::Pgon Pgon;

  FT total(0);
  MA::voronoi_triangulation_intersection_raw
  (densityT,dt,
   [&] (const Pgon &pgon,
	T::Face_handle f,
	RT::Vertex_handle v)
   {
     Tri_isector isector;

     Polygon p;
     std::vector<Vertex_handle_RT> adj;
     for (size_t i = 0; i < pgon.size(); ++i)
       {
	 p.push_back(isector.vertex_to_point(pgon[i]));
	 Tri_isector::Pgon_edge e = MA::common(pgon[i], 
					       pgon[(i+1)%pgon.size()]);
	 adj.push_back((e.type == Tri_isector::EDGE_DT) ?
		       e.edge_dt.second : 0);
       }

     size_t idv = indices[v->point()];
     
     // compute hessian
     for (size_t i = 0; i < p.size(); ++i)
       {
	 if (adj[i] == 0)
	   continue;
	 FT r = MA::integrate_1(p.edge(i), densityF[f]);
	 Vertex_handle_RT w = adj[i];
	 size_t idw = indices[w->point()];
       }

     // compute value and gradient
     FT area = MA::integrate_1(p, densityF[f]);
     FT intg = MA::integrate_3(p, [&](Point p) 
			       {
                                 return
			         (densityF[f](p) * 
			          CGAL::squared_distance(p, v->point()));
			       });
     fval = fval + area * weights[idv] - intg; 
     g[idv] = g[idv] + area;
     total += area;
   });

  std::cerr << "total = " << total  << "\n";
  //Eigen::VectorXd gg = fval.derivatives();
  //  std::cerr << "g[0] = " << g[0] << "/ " << gg(0) << "\n";
  return fval;
}

int main(int argc, const char **argv)
{
  if (argc < 2)
    return -1;

  std::map<T::Face_handle,
	   MA::Linear_function<K>> functions;
  T t;
  cimg_library::CImg<double> image(argv[1]);
  double total_mass = MA::image_to_pl_function(image, t, functions);
  boost::timer::auto_cpu_timer tm(std::cerr);

  // generate points
  size_t N = 50;
  std::vector<Point> X(N);
  VectorXd masses(N);
  std::map<Point, size_t> indices;
  for (size_t i = 0; i < N; ++i)
    {
      Point p(rr(), rr());
      X[i] = p;
      masses[i] = total_mass/N;
      indices[p] = i;
    }

  std::cerr << total_mass << "\n";

  auto eval = [&](const VectorXd &weights,
		  VectorXd &g)
    {
      return ot_eval(t, functions, X, masses, indices, weights, g);
    };

  Lbfgs lb (N);
  //lb._pgtol = 1e-6 / N;
  //lb._factr = 1e6;

  VectorXd g = VectorXd::Zero(N), weights = g;
  double f = eval(weights, g);
  size_t k = 0;

  while (1)
    {
      int r = lb.iterate((double *) &weights[0], f,
			 (double *) &g[0]);
      if (r == LBFGS_FG)
	{
	  f = eval(weights, g);
#if 0
	  auto wp = weights, wm = weights;
	  auto gg = g;
	  double eps = 1e-5;
	  wp[0] = weights[0]+eps;
	  wm[0] = weights[0]-eps;
	  double fp = eval(wp, gg), fm = eval(wm, gg);
	  std::cerr << "Gfd=" << ((fp - fm)/(2*eps)) << " / G=";
	  std::cerr << g[0] << "\n";
#endif
	}
      else if (r == LBFGS_NEW_X)
        {
          std::cerr << (++k) << ": f="
                    << f << " |df| = " << g.maxCoeff() << "\n";
        }
      else
        break;
    }
  return 0;
}

