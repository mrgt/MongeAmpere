#ifndef MA_QUADRATURE_HPP
#define MA_QUADRATURE_HPP

#include <CGAL/Simple_cartesian.h>

namespace MA
{

// order 3 formula with 6 points
template <class K, class F>
auto
integrate_albrecht_collatz(const typename CGAL::Point_2<K> &a, 
			   const typename CGAL::Point_2<K> &b, 
			   const typename CGAL::Point_2<K> &c,
			   const F &f) -> decltype(f(a))
{
  typedef typename K::FT FT;
  const FT _1_2 = FT(1)/FT(2), _1_6 = FT(1)/FT(6), _2_3 = FT(2)/FT(3);
  const FT _1_60 = FT(1)/FT(60), _9_60 = FT(9)/FT(60);
  auto  u = b-a, v = c-a;
  auto r1 = _1_60*f(a + _1_2*u + _1_2*v);
  auto r2 = _1_60*f(a + _1_2*u);
  auto r3 = _1_60*f(a + _1_2*v);
  auto r4 = _9_60*f(a + _1_6*u + _2_3*v);
  auto r5 = _9_60*f(a + _1_6*v + _2_3*u);
  auto r6 = _9_60*f(a + _1_6*u + _1_6*v);
  return CGAL::area(a,b,c)*(r1+r2+r3+r4+r5+r6);
}

template <class K, class F>
auto
integrate_midedge(const typename CGAL::Point_2<K> &a, 
		  const typename CGAL::Point_2<K> &b, 
		  const typename CGAL::Point_2<K> &c,
		  const F &f) -> decltype(f(a))
{
  typedef typename K::FT FT;
  const FT _1_2 = FT(1)/FT(2), _1_6 = FT(1)/FT(6);
  auto u = b-a, v = c-a;
  auto r1 = _1_6*f(a + _1_2*u + _1_2*v);
  auto r2 = _1_6*f(a + _1_2*u);
  auto r3 = _1_6*f(a + _1_2*v);
  return CGAL::area(a,b,c)*(r1+r2+r3);
}

double r01() 
{ 
  return double(rand() / (RAND_MAX + 1.0));
}

template <class K>
typename CGAL::Point_2<K>
rand_in_triangle(const typename CGAL::Point_2<K> &a, 
		 const typename CGAL::Point_2<K> &b, 
		 const typename CGAL::Point_2<K> &c)
{
  double r1 = sqrt(r01()), r2 = r01();
  auto A = a - CGAL::ORIGIN, B = b - CGAL::ORIGIN, C = c - CGAL::ORIGIN;
  return CGAL::ORIGIN + ((1 - r1) * A +
			 (r1 * (1 - r2)) * B +
			 (r1 * r2) * C);
}

template <class K, class F>
auto
integrate_monte_carlo(const typename CGAL::Point_2<K> &a, 
		      const typename CGAL::Point_2<K> &b, 
		      const typename CGAL::Point_2<K> &c,
		      const F &f,
		      size_t N) -> decltype(f(a))
{
  typedef typename K::FT FT;
  decltype(f(a)) r;
  for (size_t i = 0; i < N; ++i)
    {
      auto p = rand_in_triangle(a,b,c);
      r += f(p);
    }
  return CGAL::area(a,b,c)*(r/N);
}

}

#endif

