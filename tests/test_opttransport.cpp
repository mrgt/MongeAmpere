#include <MA/optimal_transport.hpp>
#include <boost/timer/timer.hpp>
#include <lbfgs.hpp>
#include <cstdlib>

#undef Success

// Suitesparse 4.3.1 does not define UF_long, which is expected by the
// Eigen wrapper classes
#include <cs.h>
#ifndef UF_long
#define  UF_long cs_long_t
#endif
#include <Eigen/SPQRSupport>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Vector_2<K> Vector;
typedef CGAL::Line_2<K> Line;
typedef CGAL::Polygon_2<K> Polygon;

typedef Eigen::SparseMatrix<FT> SparseMatrix;
typedef Eigen::SparseVector<FT> SparseVector;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::MatrixXd MatrixXd;

typedef CGAL::Delaunay_triangulation_2<K> T;

double rr() 
{ 
  return 2*double(rand() / (RAND_MAX + 1.0))-1;
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
  size_t N = 20000;
  MatrixXd X(N,2);
  VectorXd masses(N);
  for (size_t i = 0; i < N; ++i)
    {
      X(i,0) = rr();
      X(i,1) = rr();
      masses(i) = total_mass/N;
    }

  std::cerr << total_mass << "\n";

  auto f = [&](const VectorXd &x, VectorXd &g, SparseMatrix &h)
    {
      return ot_eval(t, functions, X, masses, x, g, h);
    };

  VectorXd x = VectorXd::Zero(N);
  double eps_g = 1e-6;
  size_t iteration = 1;
  do
    {
      VectorXd g;
      SparseMatrix h;
      double fx = f(x, g, h);

      if (g.norm() < eps_g)
	break;

      Eigen::SPQR<SparseMatrix> solver(h);
      VectorXd mg = -g;
      VectorXd d = solver.solve(mg);

      // choose the step length by a simple backtracking, ensuring the
      // invertibility (up to the invariance under the addition of a
      // constance) of the hessian at the next point
      double alpha = 1;
      while(1)
	{
	  VectorXd xx = x + alpha * d;
	  VectorXd gg; SparseMatrix hh;
	  f(xx,gg,hh);
	  gg = gg + masses;
	  if (gg.minCoeff() > 1e-7)
	    break;
	  alpha = alpha / 2;
	  std::cerr << alpha << "\n";
	}

      x = x + alpha * d;
      std::cerr << "it " << iteration++ << ":"
		<< " f=" << fx 
		<< " |df|=" << g.norm()
		<< " tau = " << alpha << "\n";
    }
  while (1);
  return 0;
}

