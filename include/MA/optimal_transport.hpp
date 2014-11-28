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

#ifndef MA_OPTIMAL_TRANSPORT_HPP
#define MA_OPTIMAL_TRANSPORT_HPP

#include <MA/kantorovich.hpp>
#include <Eigen/SparseCholesky>

#ifdef MA_USE_SUITESPARSE_QR
#undef Success

// Suitesparse 4.3.1 does not define UF_long, which is expected by the
// Eigen wrapper classes
#include <cs.h>
#ifndef UF_long
#define  UF_long cs_long_t
#endif
#include <Eigen/SPQRSupport>
#include <Eigen/CholmodSupport>
#endif

namespace MA
{
  struct Statistics
  {
    size_t niter;
    size_t neval;
  };

  template <class SparseMatrix, class Vector>
  Vector solve_laplacian_matrix(const SparseMatrix &h, 
				const Vector &g,
				bool verbose = false)
  {
    size_t N = h.rows(); 
    assert(N == h.cols());
    auto v = h.diagonal();
    if (v.minCoeff() == 0)
      {
	size_t i;
	std::cerr << "Error: hessian of Kantorovich's functional "
		  << "is not invertible:\n";
	std::cerr << "diag = " << v.head(10)
		  << " ... in [" << v.minCoeff()  << "," 
		  << v.minCoeff(&i) << "]\n";;
	std::cerr << "minCoeff => " << i << "\n";
      }

    // remove last row and column so that the linear system is
    // invertible
    Vector gs = g.head(N-1);
    SparseMatrix hs = h.block(0,0,N-1,N-1); // top-left submatrix
    Eigen::SimplicialLLT<SparseMatrix> solver(hs);
    Vector ds = solver.solve(gs);

    // if cannot solve with Cholesky, use QR decomposition
    double err = (hs*ds - gs).norm();
    if (err > 1e-7) // FIXME: threshold
      {
	std::cerr << "WARNING: in solve_laplacian_matrix: err=" << err << "\n";
#ifdef MA_USE_SUITESPARSE_QR
	std::cerr << "Resorting to QR decomposition";
	Eigen::SPQR<SparseMatrix> solver(hs);
	ds = solver.solve(gs);
	if (verbose)
	  std::cerr << "rank(h) = " << solver.rank() << "\n";
#endif
      }

    // assemble result
    Vector d(N);
    d.head(N-1) = ds;
    d(N-1) = 0;

    return d;
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
	std::cerr << "This is because the Laguerre cell for the "
		  << "point [" << X(i,0) << ", " << X(i,1) << "],"
		  << " i=" << i << " is empty.\n";
	return;
      }

    while (g.norm() >= eps_g && 
	   niter++ <= maxiter)
      {
	Vector d = -solve_laplacian_matrix(h,g);
	
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
		      << " |df|=" << g.norm()
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



