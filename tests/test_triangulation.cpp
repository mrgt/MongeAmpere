#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_incremental_builder_2.h>
#include <CGAL/Triangulation_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_2<K> T;
typedef CGAL::Point_2<K> Point;

int main()
{
  std::vector<Point> points;
  points.push_back(Point(0,0));
  points.push_back(Point(1,0));
  points.push_back(Point(1,1));
  points.push_back(Point(0,1));

  typedef std::array<size_t, 3> Triangle;
  std::vector<Triangle> triangles;
  triangles.push_back(Triangle({0,1,2}));
  triangles.push_back(Triangle({0,2,3}));

  T t;
  CGAL::Triangulation_incremental_builder_2<T> builder(t);

  std::vector<T::Vertex_handle > vvh;
  vvh.reserve(points.size());

  builder.begin_triangulation();
  for (auto p:points)
    vvh.push_back(builder.add_vertex(p));
  for (auto t:triangles)
    builder.add_face(vvh[t[0]], vvh[t[1]], vvh[t[2]]);
  builder.end_triangulation();
}
