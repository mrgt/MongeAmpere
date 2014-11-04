#ifndef MA_FUNCTIONS_HPP
#define MA_FUNCTIONS_HPP

#include <CImg.h>
#include <MA/quadrature.hpp>

namespace MA
{
// Compute barycentric coordinates (u, v, w) for
// point p with respect to triangle (a, b, c)
template <class Point, class FT>
void barycentric(const Point &p,
		 const Point &a,
		 const Point &b,
		 const Point &c,
		 FT &u, FT &v, FT &w)
{
  auto v0 = b - a, v1 = c - a, v2 = p - a;
  FT d00 = v0*v0;
  FT d01 = v0*v1;
  FT d11 = v1*v1;
  FT d20 = v2*v0;
  FT d21 = v2*v1;
  FT denom = d00 * d11 - d01 * d01;
  v = (d11 * d20 - d01 * d21) / denom;
  w = (d00 * d21 - d01 * d20) / denom;
  u = FT(1) - v - w;
}

template <class Point, class FT>
FT extrapolate(const Point &p,
	       const Point &a, FT fa,
	       const Point &b, FT fb,
	       const Point &c, FT fc)
{
  FT u, v, w;
  barycentric(p, a, b, c, u, v, w);
  return u * fa + v * fb + w * fc;
}

template<class K>
class Linear_function
{
  typedef typename K::Point_2 Point;
  typedef typename K::FT FT;
  FT _a, _b, _c;

public:
  Linear_function(): _a(0), _b(0), _c(0) {}
  Linear_function(const Point &p, FT fp, 
		  const Point &q, FT fq,
		  const Point &r, FT fr) 
  {
    _c = extrapolate(Point(0,0), p, fp, q, fq, r, fr); // const term
    _a = extrapolate(Point(1,0), p, fp, q, fq, r, fr) - _c;
    _b = extrapolate(Point(0,1), p, fp, q, fq, r, fr) - _c;
  }

  FT
  operator () (const CGAL::Point_2<K> &p) const
  {
    return _a * p.x() + _b * p.y() + _c;
  }

  typedef FT result_type;
};

template <class T, class Function>
double
image_to_pl_function(const cimg_library::CImg<double> &image,
		     T &t,
		     std::map<typename T::Face_handle, Function> &fs)
{
  typedef typename T::Point Point;
  std::vector<Point> grid;
  std::map<Point, double> fgrid;
  size_t n = image.width();
  size_t m = image.height();
  double dx = 2/double(n-1), x0=-1.0;
  double dy = 2/double(m-1), y0=-1.0;
  for (size_t i = 0; i < n; ++i)
    {
      for (size_t j = 0; j < m; ++j)
	{
	  Point p(x0 + i * dx,
		  y0 + j * dy);
	  grid.push_back(p);
	  fgrid[p] = image(i,m-j-1)/double(255) + 1e-3;
	}
    }
  t = T(grid.begin(), grid.end());

  double tot_orig(0);
  for (auto f = t.finite_faces_begin(); 
       f != t.finite_faces_end(); ++f)
    {
      Point p = f->vertex(0)->point(), 
	q = f->vertex(1)->point(), 
	r = f->vertex(2)->point();
      fs[f] = Function(p, fgrid[p], 
		       q, fgrid[q],
		       r, fgrid[r]);
      tot_orig += CGAL::to_double(MA::integrate_centroid(p,q,r, fs[f]));
    }
  return tot_orig;
}
}

#endif
