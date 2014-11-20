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
#include <Eigen/CholmodSupport>

namespace MA
{
  struct Statistics
  {
    size_t niter;
    size_t neval;
  };

  template <class SparseMatrix>
  void check_laplacian_matrix(const SparseMatrix &h)
  {
    size_t N = h.rows(); 
    assert(N == h.cols());
    auto v = h.diagonal();
    size_t i;
    std::cerr << "diag = " << v.head(10)
	      << " ... in [" << v.minCoeff()  << "," 
	      << v.maxCoeff(&i) << "]\n";;
    std::cerr << "maxCoeff => " << i << "\n";
  }

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
    size_t N = X.rows();
    assert(X.cols() == 2);
    assert(masses.rows() == N);
    assert(masses.cols() == 1);

    auto f = [&](const Vector &x,
		 Vector &m, // masses of the power cells
		 Vector &g, // gradient
		 SparseMatrix &h)
      {
	++neval;
	double r = kantorovich(t, functions, X, x, g, h);
	m = g;
	g = g - masses;
	return r - masses.dot(x);
      };
    

    // if no initial guess is provided, start with zero, and compute
    // function value, gradient and hessian
    if (x.size() != N)
      {
	x = Vector::Zero(N);
      }
    Vector g,m;
    SparseMatrix h;
    double fx = f(x, m, g, h); 

    // we impose a minimum weighted area for Laguerre cells during the
    // execution of the algorithm:
    // eps = min(minimum of cells areas at beginning,
    //           minimum of target areas).
    double eps0 = std::min(m.minCoeff(),
			   masses.minCoeff())/2;
    if (eps0 <= 0)
      {
	std::cerr << "Error: computed minimum mass is non-positive\n";
	size_t i;
	m.minCoeff(&i);
	std::cerr << "This is because the Laguerre cell for the point ["
		  << X(i,0) << ", " << X(i,1) << "], i=" << i << " is empty.\n";
	return;
      }

    while (g.norm() >= eps_g && 
	   niter++ <= maxiter)
      {
	check_laplacian_matrix(h);
#if 1
	// remove last row so that the linear system is invertible
	Vector gs = g.head(N-1);
	SparseMatrix hs = h.block(0,0,N-1,N-1); // top-left submatrix
	Eigen::SimplicialLDLT<SparseMatrix> solver2(h);
        //Eigen::SparseLU<SparseMatrix> solver(h);
        Eigen::SPQR<SparseMatrix> solver(hs);
        //Eigen::ConjugateGradient<SparseMatrix> solver(h);
	Vector ds = -solver.solve(Vector(gs));
	Vector d(N);
	d.head(N-1) = ds;
	d(N-1) = 0;

	Vector ds2 = -solver2.solve(Vector(gs));
	Vector d2(N);
	d2.head(N-1) = ds2;
	d2(N-1) = 0;

	std::cerr << "rank(h) = " << solver.rank() << "\n";
	// determine multiplicative constant:
	// double k = g.norm() / (h*d).norm();
	// d = k*d;
#else
        //Eigen::SimplicialLDLT<SparseMatrix> solver(h);
        //Eigen::SparseLU<SparseMatrix> solver(h);
        Eigen::SPQR<SparseMatrix> solver(h);
        //Eigen::ConjugateGradient<SparseMatrix> solver(h);
	Vector d = -solver.solve(Vector(g));
#endif
	double err = (h*d + g).norm();
	double err2 = (h*d2 + g).norm();
	assert(err < 1e-7);


	std::cerr << "err=" << err << "\n";
	std::cerr << "err2=" << err2 << "\n";
	std::cerr << "|d|=" << d.norm() << "\n";
	std::cerr << "max(d)=" << d.maxCoeff() << "\n";
	std::cerr << "min(d)=" << d.minCoeff() << "\n";
	
	// choose the step length by a simple backtracking, ensuring the
	// invertibility (up to the invariance under the addition of a
	// constant) of the hessian at the next point
	double alpha = 1;
	Vector x0 = x;
	double n0 = g.norm();
	size_t nlinesearch = 0;

	while(1)
	  {
	    x = x0 + alpha * d;
	    fx = f(x,m,g,h); 
	    if (m.minCoeff() >= eps0 &&
		g.norm() <= (1-alpha/2)*n0)
	      break;
	    alpha *= .5;
	    if (verbose)
	      {
		std::cerr << "subit " << niter << "." << (nlinesearch++)
			  << ": min(masses)=" << m.minCoeff() << "\n";
	      }
	  }
	if (verbose)
	  {
	    std::cerr << "it " << niter << ":"
		      << " f=" << fx 
		      << " |df|=" << (g-masses).norm()
		      << " min(m)=" << masses.minCoeff()
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



