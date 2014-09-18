#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <MA/voronoi_polygon_intersection.hpp>
#include <MA/misc.hpp>

#include <cstdlib>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Vector_2<K> Vector;
typedef CGAL::Line_2<K> Line;
typedef CGAL::Polygon_2<K> Polygon;

typedef CGAL::Regular_triangulation_filtered_traits_2<K> Traits;
typedef CGAL::Regular_triangulation_2<Traits> RT;

double rr() 
{ 
  return 2*double(rand() / (RAND_MAX + 1.0))-1;
}

#include <boost/timer/timer.hpp>

int main()
{
  Polygon P;
  P.push_back(Point(-1,-1));
  P.push_back(Point(1,-1));
  P.push_back(Point(1,1));
  P.push_back(Point(-1,1));

  std::vector<Point> pts;
  for (size_t i = 0; i < 100000; ++i)
    pts.push_back(RT::Weighted_point(Point(rr(), rr()), 0));
  RT rt (pts.begin(), pts.end());

  boost::timer::auto_cpu_timer t(std::cerr);
  MA::ps_begin(std::cout);
  double A (0);
  for (auto v = rt.finite_vertices_begin();
       v != rt.finite_vertices_end(); ++v)
  {
    Polygon R = MA::voronoi_polygon_intersection(P, rt, v);
    MA::ps_polygon(std::cout, R, 0.001);
    A += CGAL::to_double(R.area());
  }
  MA::ps_end(std::cout);
  std::cerr << A << "\n";
}

