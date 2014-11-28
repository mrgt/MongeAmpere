#include <MA/lloyd.hpp>
#include <boost/timer/timer.hpp>
#include <cstdlib>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Vector_2<K> Vector;
typedef CGAL::Line_2<K> Line;
typedef CGAL::Polygon_2<K> Polygon;
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
  size_t N = atoi(argv[2]);

  size_t niter = 100;
  if (argc > 3)
    niter = atoi(argv[3]);

  std::map<T::Face_handle,
	   MA::Linear_function<K>> functions;

  // generate points
  MatrixXd X(N,2);
  for (size_t i = 0; i < N; ++i)
    {
      X(i,0) = rr();
      X(i,1) = rr();
    }

  T t;
  cimg_library::CImg<double> image(argv[1]);
  MA::image_to_pl_function(image, t, functions);
  boost::timer::auto_cpu_timer tm(std::cerr);
  
  VectorXd masses = VectorXd::Zero(N);
  VectorXd weights = VectorXd::Zero(N);
  MatrixXd centroids = MatrixXd::Zero(N,2);

  for (size_t it = 0; it < niter; ++it)
    {
      MA::lloyd(t, functions, X, weights, centroids, masses);
      X = centroids;
    }

  // output to stdio
  for (size_t i = 0; i < N; ++i)
    std::cout << X(i,0) << " " << X(i,1) << "\n";
}
