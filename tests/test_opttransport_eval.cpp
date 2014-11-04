#include <MA/optimal_transport.hpp>
#include <boost/timer/timer.hpp>
#include <lbfgs.hpp>
#include <cstdlib>

#undef Success
#include <Eigen/Sparse>

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

int main()
{
  srand(time(NULL));

  std::map<T::Face_handle,
	   MA::Linear_function<K>> functions;
  T t;
  cimg_library::CImg<double> image(argv[1]);
  double total_mass = MA::image_to_pl_function(image, t, functions);
  boost::timer::auto_cpu_timer tm(std::cerr);

  // generate points
  size_t N = 50;
  MatrixXd X(N,2);
  for (size_t i = 0; i < N; ++i)
    {
      X(i,0) = rr();
      X(i,1) = rr();
    }

  auto eval = [&](const VectorXd &weights,
		  VectorXd &g,
		  SparseMatrix &h)
    {
      return MA::kantorovich(t, functions, X, weights, g, h);
    };


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
