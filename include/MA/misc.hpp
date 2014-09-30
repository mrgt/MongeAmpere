#ifndef MA_MISC_HPP
#define MA_MISC_HPP

#include <CGAL/Polygon_2.h>

namespace MA
{
  void ps_begin(std::ostream &os)
  {
    os << "%!PS\n";
    os << " matrix currentmatrix /originmat exch def\n"
       << "/umatrix {originmat matrix concatmatrix setmatrix} def\n"
       << "[28.3465 0 0 28.3465 10.5 100.0] umatrix\n\n";
  }

  void ps_end(std::ostream &os)
  {
    os << "showpage\n";
  }

  template <class K>
  void ps_polygon(std::ostream &os,
		  const typename CGAL::Polygon_2<K> &P,
		  double linewidth = 0.1,
		  double r=0, double g=0, double b=0,
		  bool filled = false)
  {
    size_t n = P.size();
    if (n == 0)
      return;
    
    os << linewidth << " setlinewidth\n";
    os << r << " " << g << " " << b << " setcolor\n";
    double tx = 10.5, ty = 14.15;
    os << (10*P[n-1].x()+tx) << " "
       << (10*P[n-1].y()+ty) << " newpath moveto\n";
    for (auto E = P.vertices_begin(); E != P.vertices_end(); ++E)
      os << (10*E->x()+tx) << " " << (10*E->y()+ty) << " lineto\n";
    os << (filled ? "fill" : "stroke") << "\n\n";
  }
}

#endif

