#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <MA/voronoi_triangulation_intersection.hpp>
#include <MA/voronoi_polygon_intersection.hpp>
#include <MA/quadrature.hpp>
#include <MA/misc.hpp>
#include <MA/functions.hpp>
#include <boost/timer/timer.hpp>

#include <cstdlib>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Vector_2<K> Vector;
typedef CGAL::Line_2<K> Line;
typedef CGAL::Polygon_2<K> Polygon;
typedef K::FT FT;

typedef CGAL::Delaunay_triangulation_2<K> Delaunay_triangulation;
typedef Delaunay_triangulation Triangulation;

double rr() 
{ 
  return 2*double(rand() / (RAND_MAX + 1.0))-1;
}

int main(int argc, const char **argv)
{
  if (argc < 2)
    return -1;

  std::map<Triangulation::Face_handle, MA::Linear_function<K>> functions;
  Triangulation t;
  cimg_library::CImg<double> image(argv[1]);
  FT tot_orig = MA::image_to_pl_function(image, t, functions);
  boost::timer::auto_cpu_timer tm(std::cerr);

  // generate points and their Delaunay triangulation
  std::vector<Point> pts;
  for (size_t i = 0; i < 5000; ++i)
    pts.push_back(Point(rr(), rr()));
  Delaunay_triangulation dt (pts.begin(), pts.end());

  // compute the integrals over Voronoi cells
  std::map<Delaunay_triangulation::Vertex_handle, FT> integrals;
  FT total(0);
  
  MA::voronoi_triangulation_intersection
  (t,dt,
   [&] (const CGAL::Polygon_2<K> &p,
	Triangulation::Face_handle f,
	Delaunay_triangulation::Vertex_handle v)
   {
     FT z = MA::integrate_1<FT>(p, FT(0), functions[f]);
     integrals[v] += z; 
     total += z;
   });
  
  Polygon P;
  P.push_back(Point(-1,-1));
  P.push_back(Point(1,-1));
  P.push_back(Point(1,1));
  P.push_back(Point(-1,1));

  MA::ps_begin(std::cout);
  for (auto v = dt.finite_vertices_begin();
       v != dt.finite_vertices_end(); ++v)
    {
      Polygon R = MA::voronoi_polygon_intersection(P, dt, v);
      double a = CGAL::to_double(R.area());
      double ig = integrals[v]/a;
      MA::ps_polygon(std::cout, R,
		     0.001, 
		     ig, ig, ig, true);
    }

  std::cerr << total << " vs "
	    << tot_orig << "\n";
  MA::ps_end(std::cout);
  return 0;
}

