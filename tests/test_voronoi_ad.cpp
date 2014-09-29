#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>


#include <MA/Autodiff_nt.hpp>
#include <MA/voronoi_polygon_intersection.hpp>
#include <MA/misc.hpp>

#include <cstdlib>

typedef AD FT;
typedef CGAL::Filtered_kernel<CGAL::Simple_cartesian<AD>> K;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Vector_2<K> Vector;
typedef CGAL::Line_2<K> Line;
typedef CGAL::Polygon_2<K> Polygon;


typedef CGAL::Delaunay_triangulation_2<K> DT;

double rr() 
{ 
  return 2*double(rand() / (RAND_MAX + 1.0))-1;
}

int main()
{
  Polygon cross;
  
  cross.push_back(Point(+1,-1));
  cross.push_back(Point(+2,-1));
  cross.push_back(Point(+2,+1));
  cross.push_back(Point(+1,+1));
  cross.push_back(Point(+1,+2));
  cross.push_back(Point(-1,+2));
  cross.push_back(Point(-1,+1));
  cross.push_back(Point(-2,+1));
  cross.push_back(Point(-2,-1));
  cross.push_back(Point(-1,-1));
  cross.push_back(Point(-1,-2));
  cross.push_back(Point(+1,-2));

  DT dt;
  dt.insert(Point(0,0));
  for (size_t i = 0; i < 10000; ++i)
    dt.insert(Point(rr(), rr()));

  MA::ps_begin(std::cout);
  double A (0);
  for (auto v = dt.finite_vertices_begin();
       v != dt.finite_vertices_end(); ++v)
  {
    Polygon R = MA::voronoi_polygon_intersection(cross, dt, v);
    MA::ps_polygon(std::cout, R, 0.01);
    A += CGAL::to_double(R.area());
  }
  MA::ps_end(std::cout);
  std::cerr << A << "\n";
}

