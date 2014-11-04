#include <MA/optimal_transport.hpp>
#include <boost/timer/timer.hpp>
#include <lbfgs.hpp>
#include <cstdlib>
#include <fstream>

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

  // generate points
  //size_t ns[] = {10, 30, 50, 100, 200, 300, 400, 500};
  size_t ns[] = {300, 400, 500};
  std::ofstream os("bench_opttransport.res");
  
  for (size_t n : ns)
    {
      size_t N = n*n;

#if 0
      double dx = 2/double(n-1), x0=-1.0;
      double dy = 2/double(n-1), y0=-1.0;

      MatrixXd X(N,2);
      VectorXd masses(N);
      size_t cur = 0;
      for (size_t i = 0; i < n; ++i)
	{
	  for (size_t j = 0; j < n; ++j, ++cur)
	    {
	      X(cur,0) = x0 + i * dx + rr()/(2.0*n);
	      X(cur,1) = y0 + j * dy + rr()/(2.0*n);
	      masses(cur) = total_mass/N;
	    }
	}
#else
  MatrixXd X(N,2);
  VectorXd masses(N);
  for (size_t i = 0; i < N; ++i)
    {
      X(i,0) = rr();
      X(i,1) = rr();
      masses(i) = total_mass/N;
    }
#endif

      VectorXd res;
      boost::timer::cpu_timer tm;
      MA::Statistics stats;
      std::cerr << "CASE n=" << n << " / N=" << N << "\n";
      MA::ot_solve(t, functions, X, masses, res, 
		   1e-7, 100, true, &stats);
      auto ut = tm.format(boost::timer::default_places,
			  "%u");
      std::cerr << "TIME=" << ut <<  "s\n\n";
      os << "N=" << N
	 << " niter=" << stats.niter
	 << " neval=" << stats.neval
	 << " cpu=" << ut << "\n";
      std::flush(os);
    }
  return 0;
}

