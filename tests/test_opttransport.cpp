#define EIGEN_CHOLMOD_SUPPORT
#include <MA/optimal_transport.hpp>
#include <boost/timer/timer.hpp>
#include <lbfgs.hpp>
#include <cstdlib>

#undef Success

//#include <Eigen/SparseLU>
// fix for suitesparse
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

//typedef CGAL::Simple_cartesian<AD> K_ad;

typedef Eigen::Triplet<FT> Triplet;
typedef Eigen::SparseMatrix<FT> SparseMatrix;
typedef Eigen::SparseVector<FT> SparseVector;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::MatrixXd MatrixXd;

typedef CGAL::Delaunay_triangulation_2<K> T;

double rr() 
{ 
  return 2*double(rand() / (RAND_MAX + 1.0))-1;
}

template <class F>
void eval_gradient_fd(F f, const VectorXd &x, VectorXd &g,
		      FT eps = 1e-5)
{
  size_t N = x.size();
  g = VectorXd::Zero(N);

  for (size_t i = 0; i < N; ++i)
    {
      VectorXd xp = x, xm = x;
      VectorXd gg = g;
      SparseMatrix hh;
      xp[i] += eps;
      xm[i] -= eps;
      FT fp = f(xp, gg, hh), fm = f(xm, gg, hh);
      g[i] = (fp - fm)/(2*eps);
    }
}

template <class F>
void eval_hessian_fd(F f, const VectorXd &x, MatrixXd &h,
		      FT eps = 1e-5)
{
  size_t N = x.size();
  h = MatrixXd::Zero(N,N);
  for (size_t i = 0; i < N; ++i)
    {
      VectorXd xp = x, xm = x;
      VectorXd gp = VectorXd::Zero(N), gm = VectorXd::Zero(N);
      SparseMatrix hh;
      xp[i] += eps;
      xm[i] -= eps;
      FT fp = f(xp, gp, hh), fm = f(xm, gm, hh);
      h.row(i) = (gp - gm)/(2*eps);
    }
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

  auto eval = [&](const VectorXd &weights,
		  VectorXd &g,
		  SparseMatrix &h)
    {
      return ot_eval(t, functions, X, masses, weights, g, h);
    };
  auto evalg = [&](const VectorXd &x, VectorXd &g)
    {
      SparseMatrix h;
      return eval(x,g,h);
    };

  VectorXd x = VectorXd::Zero(N);
  double eps_g = 1e-6;
  size_t iteration = 1;
  do
    {
      VectorXd g;
      SparseMatrix h;
      double fx = eval(x, g, h);

      if (g.norm() < eps_g)
	break;

#if 0
      // Eigen::SimplicialLDLT<SparseMatrix> solver;      
      Eigen::SPQR<SparseMatrix> solver;
      Eigen::VectorXd dH = h.diagonal();
      FT mindiag = dH.minCoeff();

      // Attempt repeated Cholesky factorization until the Hessian
      // becomes positive semidefinite.
      FT beta = 1e-5;
      FT tau = (mindiag > 0) ? 0 : (-mindiag + beta);
      size_t factorizations = 0;
      while (true)
	{
	  // Add tau*I to the Hessian.
	  if (tau > 0)
	    {
	      std::cerr << tau << "\n";
	      for (size_t i = 0; i < N; ++i)
		{
		  int ii = static_cast<int>(i);
		  h.coeffRef(ii, ii) = dH(i) + tau;
		}
	    }
	  // Attempt Cholesky factorization.
	  solver.compute(h);
	  bool success = (solver.info() == Eigen::Success);
	  factorizations++;
	  // Check for success.
	  if (success)
	    break;
	  tau = std::max(2*tau, beta);
	  assert(factorizations <= 100);
	}
      VectorXd d = solver.solve(-g);
#endif
#if 1
      Eigen::SPQR<SparseMatrix> solver(h);
      VectorXd mg = -g;
      VectorXd d = solver.solve(mg);
#else
      // SparseMatrix id(N,N);
      // id.setIdentity();
      // h = h + 1e-9 * id;

      Eigen::CholmodDecomposition<SparseMatrix> solver;
      solver.setShift(1e-8);
      solver.compute(h);
      VectorXd mg = -g;
      VectorXd d = solver.solve(mg);
#endif

#if 0
      double alpha = MA::perform_Wolfe_linesearch(evalg, x, fx, g, d, 1);
      if (alpha < 1e-15)
	{
	  d = -g;
	  alpha = MA::perform_Wolfe_linesearch(evalg, x, fx, g, d, 1);
	}
#endif

      double alpha = 1;
      while(1)
	{
	  VectorXd xx = x + alpha * d;
	  VectorXd gg;
	  evalg(xx,gg);
	  gg = gg + masses;
	  if (gg.minCoeff() > 1e-7)
	    break;
	  alpha = alpha / 2;
	  std::cerr << alpha << "\n";
	}

      x = x + alpha * d;
      std::cerr << "iteration " << iteration++
		<< " f=" << fx 
		<< " |df|=" << g.norm()
		<< " tau = " << alpha << "\n";
    }
  while (1);
  return 0;

#if 1
  {
    srand(time(NULL));
    VectorXd x = VectorXd::Random(N), gfd, gx;
    SparseMatrix hx;
    MatrixXd hfd;
    eval(x, gx, hx);
    eval_gradient_fd(eval, x, gfd, 1e-5);
    eval_hessian_fd(eval, x, hfd, 1e-5);
    
    std::cerr << (gx-gfd).norm() << "\n";
    std::cerr << (hfd-MatrixXd(hx)).norm() << "\n";
    std::cerr << hx << "\n";
  }
#endif

  Lbfgs lb (N);
  //lb._pgtol = 1e-6 / N;
  //lb._factr = 1e6;

  VectorXd g = VectorXd::Zero(N), weights = g;
  SparseMatrix h;
  double f = eval(weights, g, h);
  size_t k = 0;

  while (1)
    {
      int r = lb.iterate((double *) &weights[0], f,
			 (double *) &g[0]);
      if (r == LBFGS_FG)
	{
	  f = eval(weights, g, h);
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

