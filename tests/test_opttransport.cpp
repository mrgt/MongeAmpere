#include <MA/optimal_transport.hpp>
#include <boost/timer/timer.hpp>
#include <lbfgs.hpp>
#include <cstdlib>

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
  if (argc < 3)
    return -1;
  size_t N = atoi(argv[1]);

  std::map<T::Face_handle,
	   MA::Linear_function<K>> functions;

  // generate points
  MatrixXd X(N,2);
  VectorXd masses(N);
  for (size_t i = 0; i < N; ++i)
    {
      X(i,0) = rr();
      X(i,1) = rr();
      masses(i) = 1.0;
    }

  VectorXd res;
  for (size_t i = 2; i < argc; ++i)
  {
      T t;
      cimg_library::CImg<double> image(argv[i]);
      double total_mass = MA::image_to_pl_function(image, t, functions);
      for (size_t i = 0; i < N; ++i)
	masses(i) = total_mass/N;
      boost::timer::auto_cpu_timer tm(std::cerr);
      MA::ot_solve(t, functions, X, masses, res);
  }
  return 0;
}

