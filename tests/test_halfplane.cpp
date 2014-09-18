#include <MA/polygon_intersection.hpp>
#include <MA/misc.hpp>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Vector_2<K> Vector;
typedef CGAL::Line_2<K> Line;
typedef CGAL::Polygon_2<K> Polygon;

#include <cstdlib>
double rr() 
{ 
  return double(rand() / (RAND_MAX + 1.0));
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

  srand(time(NULL));
  Polygon isect;
  for (size_t i = 0; i < 100000; ++i)
    {
      Line L(Point(rr(),rr()), Vector(rr(),rr()));
      isect = Polygon();
      MA::polygon_halfplane_intersection (cross, L, isect);
    }
  MA::ps_begin(std::cout);
  MA::ps_polygon(std::cout, isect);
  MA::ps_end(std::cout);
}

