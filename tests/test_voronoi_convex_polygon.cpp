#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <MA/voronoi_convex_polygon_intersection.hpp>
#include <MA/misc.hpp>

#include <cstdlib>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Vector_2<K> Vector;
typedef CGAL::Line_2<K> Line;
typedef CGAL::Polygon_2<K> Polygon;


typedef CGAL::Delaunay_triangulation_2<K> DT;

double rr() 
{ 
  return 2*double(rand() / (RAND_MAX + 1.0))-1;
}

#include <boost/timer/timer.hpp>

int main()
{
  std::vector<Line> square;
  square.push_back(Line(0,1,1));
  square.push_back(Line(0,-1,1));
  square.push_back(Line(1,0,1));
  square.push_back(Line(-1,0,1));
  
  // std::vector<size_t> ch;
  // MA::halfplanes_intersection(square, Point(0,0), ch);
  // for (auto i:ch)
  //   std::cerr << i << "\n";
  std::vector<Point> pts;
  for (size_t i = 0; i < 100000; ++i)
    pts.push_back(Point(rr(), rr()));
  DT dt (pts.begin(), pts.end());

  boost::timer::auto_cpu_timer t (std::cerr);
  MA::ps_begin(std::cout);
  double A (0);
  for (auto v = dt.finite_vertices_begin();
       v != dt.finite_vertices_end(); ++v)
  {
    Polygon R;
    MA::voronoi_convex_polygon_intersection(square, dt, v, R);
    MA::ps_polygon(std::cout, R, 0.01);
    A += CGAL::to_double(R.area());
  }
  MA::ps_end(std::cout);
  std::cerr << A << "\n";
}

