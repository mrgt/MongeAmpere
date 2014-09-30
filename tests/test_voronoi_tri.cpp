#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <MA/voronoi_triangulation_intersection.hpp>
#include <MA/misc.hpp>

#include <cstdlib>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Vector_2<K> Vector;
typedef CGAL::Line_2<K> Line;
typedef CGAL::Polygon_2<K> Polygon;
typedef K::FT FT;

typedef CGAL::Delaunay_triangulation_2<K> DT;
typedef DT T;

double rr() 
{ 
  return 2*double(rand() / (RAND_MAX + 1.0))-1;
}

#include <boost/timer/timer.hpp>

struct F
{
  double a;
  F(): a(0) {}

  template <class K, class FH, class VH>
  void operator() (const CGAL::Polygon_2<K> &p,
		   FH tri, VH v)
  {
    MA::ps_polygon(std::cout, p, 0.001);
    a = a + p.area();
  }
};


int main()
{
  std::vector<Point> pts;
  for (size_t i = 0; i < 100; ++i)
    pts.push_back(Point(rr(), rr()));
  DT dt (pts.begin(), pts.end());

  // 
  std::vector<Point> grid;
  int n = 1;
  for (int i = -n; i <= n; ++i)
    {
      for (int j = -n; j <= n; ++j)
	{
	  Point p(double(i)/double(n),
		  double(j)/double(n));
	  grid.push_back(p);
	}
    }
  T t (grid.begin(), grid.end());

  F f;
  boost::timer::auto_cpu_timer tm(std::cerr);
  MA::ps_begin(std::cout);
  MA::voronoi_triangulation_intersection(t,dt,f);
  MA::ps_end(std::cout);
  std::cerr << "area = " << f.a << "\n";
}

