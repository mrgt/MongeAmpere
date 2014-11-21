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

void barycenters(const T &t,
		 const MatrixXd &X,
		 const VectorXd &weights,
		 std::vector<Vector> &barys)
{
  size_t N = weights.size();
  typedef CGAL::Regular_triangulation_filtered_traits_2<K> RT_Traits;
  typedef CGAL::Regular_triangulation_2<RT_Traits> RT;
  typedef RT::Vertex_handle Vertex_handle_RT;
  typedef RT::Weighted_point Weighted_point;
  typedef typename CGAL::Point_2<K> Point;
  
  std::vector<Weighted_point> Xw(N);
  std::map<Point,size_t> indices;
  for (size_t i = 0; i < N; ++i)
    {
      Point p(X(i,0), X(i,1));
      indices[p] = i;
      Xw[i] = Weighted_point(p, weights(i));
    }
  RT dt (Xw.begin(), Xw.end());
  
  // compute the barycenters 
  typedef MA::Voronoi_intersection_traits<K> Traits;
  typedef typename MA::Tri_intersector<T,RT,Traits> Tri_isector;  
  typedef typename Tri_isector::Pgon Pgon;
  
  typedef Eigen::Triplet<FT> Triplet;
  std::vector<Triplet> htri;
  
  barys.clear();
  barys.resize(N,Vector(0,0));
  std::vector<FT> areas(N,0);
  
  MA::voronoi_triangulation_intersection
    (t,dt,
     [&] (const Polygon &P,
	  typename T::Face_handle f,
	  Vertex_handle_RT v)
     {
       size_t idv = indices[v->point()];
       Vector bary = MA::integrate_1(P, [&](Point p)
				     {
				       return Vector(p-CGAL::ORIGIN);
				     });
       FT area = P.area();
       areas[idv] = areas[idv] + area;
       barys[idv] = barys[idv] + bary;
     });
  for (size_t i = 0; i < N; ++i)
    barys[i] = barys[i] /areas[i];
}


int main(int argc, const char **argv)
{
  std::map<T::Face_handle,
	   MA::Linear_function<K>> functions;
  T t;
  cimg_library::CImg<double> image(2,2);
  image.fill(255);
  double total_mass = MA::image_to_pl_function(image, t, functions);
  boost::timer::auto_cpu_timer tm(std::cerr);
  srand(time(NULL));

  // generate points
  size_t N = 10000;
  MatrixXd X(N,2);
  VectorXd masses(N);
  for (size_t i = 0; i < N; ++i)
    {
      X(i,0) = rr()/1.1;
      X(i,1) = rr()/1.1;
      masses(i) = total_mass/N;
    }

  double eps = 0.03;
  for (size_t i = 0; i < 100; ++i)
    {
      VectorXd weights;
      std::vector<Vector> barys;
      std::cerr << "ITERATION " << i << "\n";
      MA::ot_solve(t, functions, X, masses, weights);
      std::cerr << "\n\n";
      barycenters(t, X,  weights, barys);

      char fname[256];
      snprintf(fname, 255, "X%d", int(i));
      std::ofstream ofs(fname);
      
      for (size_t k = 0; k < N; ++k)
	{
	  ofs << X(k,0) << " " << X(k,1) << "\n";
	  X(k,0) = X(k,0) + eps * (X(k,0) - barys[k].x());
	  X(k,1) = X(k,1) + eps * (X(k,1) - barys[k].y());
	}
    }
  
      return 0;
}

