#ifndef MA_OPTIMAL_TRANSPORT_HPP
#define MA_OPTIMAL_TRANSPORT_HPP

#include <MA/kantorovich.hpp>
#include <Eigen/IterativeLinearSolvers>

#undef Success

// Suitesparse 4.3.1 does not define UF_long, which is expected by the
// Eigen wrapper classes
#include <cs.h>
#ifndef UF_long
#define  UF_long cs_long_t
#endif
#include <Eigen/SPQRSupport>


namespace MA
{
  struct Statistics
  {
    size_t niter;
    size_t neval;
  };

  template <class T, class Functions, class Matrix,
	    class Vector>
  void ot_solve(const T &t,
		const Functions &functions,
		const Matrix &X,
		const Vector &masses,
		Vector &x, // result and initial guess
		double eps_g = 1e-7,
		size_t maxiter = 100,
		bool verbose = true,
		struct Statistics *stats = 0)

  {
    typedef Eigen::SparseMatrix<double> SparseMatrix;

    size_t neval = 0, niter = 0;
    auto f = [&](const Vector &x, Vector &g, SparseMatrix &h)
      {
	++neval;
	return kantorovich(t, functions, X, x, g, h) - masses.dot(x);
      };
    
    size_t N = X.rows();
    assert(X.cols() == 2);
    assert(masses.rows() == N);
    assert(masses.cols() == 1);

    // if no initial guess is provided, start with zero, and compute
    // function value, gradient and hessian
    x = (x.size() == N) ? x : Vector::Zero(N);
    Vector g;
    SparseMatrix h;
    double fx = f(x, g, h); 

    // we impose a minimum weighted area for Laguerre cells during the
    // execution of the algorithm:
    // eps = min(minimum of cells areas at beginning,
    //           minimum of target areas).
    double eps0 = std::min(g.minCoeff(),
			   masses.minCoeff())/2;
    assert(eps0 > 0);

    while ((g-masses).norm() >= eps_g && 
	   niter++ <= maxiter)
      {
        Eigen::SPQR<SparseMatrix> solver(h);
        //Eigen::ConjugateGradient<SparseMatrix> solver(h);
	Vector d = solver.solve(Vector(masses-g));
	double mean = d.mean();
	for (size_t i = 0; i < d.size(); ++i)
	  d(i) -= mean;
	double err = (h*d - (masses-g)).norm();
	assert(err < 1e-7);
	std::cerr << "err=" << err << "\n";
	std::cerr << "|d|=" << d.norm() << "\n";
	std::cerr << "max(d)=" << d.maxCoeff() << "\n";
	std::cerr << "min(d)=" << d.minCoeff() << "\n";
	double kappa = .5;
	
	// choose the step length by a simple backtracking, ensuring the
	// invertibility (up to the invariance under the addition of a
	// constant) of the hessian at the next point
	double alpha = 1;
	Vector x0 = x;
	double n0 = (g-masses).norm();
	size_t nlinesearch = 0;

	while(1)
	  {
	    x = x0 + alpha * d;
	    fx = f(x,g,h); 
	    if (g.minCoeff() >= eps0 &&
		(g-masses).norm() <= (1-alpha/2)*n0)
	      break;
	    alpha *= kappa;
	    if (verbose)
	      {
		std::cerr << "subit " << niter << "." << (nlinesearch++)
			  << ": min(|g|)=" << g.minCoeff() << "\n";
	      }
	  }
	if (verbose)
	  {
	    std::cerr << "it " << niter << ":"
		      << " f=" << fx 
		      << " |df|=" << (g-masses).norm()
		      << " tau = " << alpha 
		      << " eval = " << neval << "\n";
	  }
      }
    
    if (stats)
      {
	stats->niter = niter;
	stats->neval = neval;
      }
  }
}

#endif



