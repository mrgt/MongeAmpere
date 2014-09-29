#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <MA/Autodiff_nt.hpp>
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
#include <Eigen/Dense>

typedef AD FT;
typedef CGAL::Filtered_kernel<CGAL::Simple_cartesian<AD>> K;

//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Vector_2<K> Vector;
typedef CGAL::Line_2<K> Line;
typedef CGAL::Polygon_2<K> Polygon;
//typedef K::FT FT;

typedef CGAL::Regular_triangulation_filtered_traits_2<K> Traits;
typedef CGAL::Regular_triangulation_2<Traits> Regular_triangulation;
typedef Regular_triangulation::Weighted_point Weighted_point;

typedef CGAL::Delaunay_triangulation_2<K> Triangulation;

double rr() 
{ 
  return 2*double(rand() / (RAND_MAX + 1.0))-1;
}

template <class Functions>
double ot_eval (const Triangulation &densityT,
		Functions &densityF,
		const std::vector<Point> &X,
		const std::vector<double> &masses,
		std::map<Point,size_t> indices,
		const std::vector<double> &w,
		std::vector<double> &g)
{
  size_t N = X.size();
  std::vector<Weighted_point> Xw(N);
  std::vector<FT> weights;
  for (size_t i = 0; i < N; ++i)
    {
      weights.push_back(AD(w[i], N, i));
      Xw[i] = Weighted_point(X[i], weights[i]);
    }
  Regular_triangulation dt (Xw.begin(), Xw.end());

  // compute the linear part of the function
  FT fval(0), total(0), tmass(0);
  for (size_t i = 0; i < N; ++i)
    {
      fval = fval - masses[i]*weights[i];
      g[i] = - masses[i];
      tmass += masses[i];
    }

  // compute the quadratic part
  MA::voronoi_triangulation_intersection
  (densityT,dt,
   [&] (const Polygon &p,
	Triangulation::Face_handle f,
	Regular_triangulation::Vertex_handle v)
   {
     FT area = MA::integrate_1(p, densityF[f]);
     FT intg = MA::integrate_3(p, [&](Point p) 
			       {
                                 return
			         (densityF[f](p) * 
			          CGAL::squared_distance(p, v->point()));
			       });
     size_t idx = indices[v->point()];
     fval = fval + area * weights[idx] - intg; 
     g[idx] = g[idx] + area.value();
     total += area;
   });

  std::cerr << "total = " << total << "/ " << tmass << "\n";
  Eigen::VectorXd gg = fval.derivatives();
  std::cerr << "g[0] = " << g[0] << "/ " << gg(0) << "\n";
  return fval.value();
}

template <class FT>
FT norm(const std::vector<FT> &v)
{
  FT r (0);
  for (auto p:v)
    r = std::max(r, fabs(p));
  return r;
}

int main(int argc, const char **argv)
{
  if (argc < 2)
    return -1;

  std::map<Triangulation::Face_handle, MA::Linear_function<K>> functions;
  Triangulation t;
  cimg_library::CImg<double> image(argv[1]);
  double total_mass = MA::image_to_pl_function(image, t, functions);
  boost::timer::auto_cpu_timer tm(std::cerr);

  // generate points
  size_t N = 50;
  std::vector<Point> X(N);
  std::vector<double> masses(N);
  std::map<Point, size_t> indices;
  for (size_t i = 0; i < N; ++i)
    {
      Point p(rr(), rr());
      X[i] = p;
      masses[i] = total_mass/N;
      indices[p] = i;
    }

  std::cerr << total_mass << "\n";

  auto eval = [&](const std::vector<double> &weights,
		  std::vector<double> &g)
    {
      return ot_eval(t, functions, X, masses, indices, weights, g);
    };

  Lbfgs lb (N);
  //lb._pgtol = 1e-6 / N;
  //lb._factr = 1e6;

  std::vector<double> g (N,0), weights(N,0);
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
                    << f << " |df| = " << norm(g) << "\n";
        }
      else
        break;
    }
  return 0;
}

