#ifndef MA_MISC_HPP
#define MA_MISC_HPP

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
		  double linewidth = 0.1)
  {
    size_t n = P.size();
    if (n == 0)
      return;
    
    os << linewidth << " setlinewidth\n";
    double tx = 10.5, ty = 14.15;
    os << (P[n-1].x()+tx) << " "
       << (P[n-1].y()+ty) << " newpath moveto\n";
    for (auto E = P.vertices_begin(); E != P.vertices_end(); ++E)
      os << (E->x()+tx) << " " << (E->y()+ty) << " lineto\n";
    os << "0 setgray\n";
    os << "stroke\n\n";
  }
}

#endif

