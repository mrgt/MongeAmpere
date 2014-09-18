#ifndef CONVEX_HULL_CONSTRUCTED
#define CONVEX_HULL_CONSTRUCTED
#include <CGAL/convex_hull_2.h>
#include <map>

namespace MA
{
   template <class K>
   typename CGAL::Point_2<K>
   line_to_dual_point(const typename CGAL::Line_2<K> &line,
		      const typename CGAL::Point_2<K> &origin)
   {
      typename K::FT newl = line.c() + (origin.x() * line.a() +
					origin.y() * line.b());
      return typename CGAL::Point_2<K>(-line.a() / newl, 
				       -line.b() / newl);
   }

  template <class K>
  void
  halfplanes_intersection
  (const std::vector<typename CGAL::Line_2<K>> &lines,
   const typename CGAL::Point_2<K> &origin,
   std::vector<size_t> &out)
  {
     typedef typename CGAL::Point_2<K> Point;
     std::vector<Point> points, ch;
     std::map<Point, size_t> pid;
     for (size_t i = 0; i < lines.size(); ++i)
     {
	auto p = line_to_dual_point(lines[i], origin);
	pid[p] = i;
	points.push_back(p);
     }
     convex_hull_2(points.begin(), points.end(),
		   std::back_inserter(ch));
     out.clear();
     for (auto p:ch)
	out.push_back(pid[p]);
  }
}


#endif
