#ifndef CONVEX_HULL_CONSTRUCTED
#define CONVEX_HULL_CONSTRUCTED

#include <vector>
#include <algorithm>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_2 Point;
typedef K::Line_2 Line;
typedef K::Vector_2 Vector;

template <class AK, class EK>
class Constructed_point_2
{
  typedef typename AK::Point_2 Approximate_point_2;
  typedef typename EK::Point_2 Exact_point_2;
 public:
  virtual Approximate_point_2 approx();
};

template <class NT>
int orientation(NT x1, NT y1,
		NT x2, NT y2,
		NT x3, NT y3)
{
  const NT d1 = x1*y2 - x2*y1;
  const NT d2 = x3*y1 - x1*y3;
  const NT d3 = x2*y3 - x3*y2;
  return sign(d1+d2+d3);
}

template <class Constructed_point>
int orientation(const Constructed_point &p1, 
		const Constructed_point &p2,
		const Constructed_point &p3)
{
  int sign;
  try
    {
      sign = orientation(p1.x_interval(), p1.y_interval(),
			 p2.x_interval(), p2.y_interval(),
			 p3.x_interval(), p3.y_interval());
    } 
  catch(...)
    {
      std::cerr << "here!\n";
      sign = orientation(p1.x_exact(), p1.y_exact(),
			 p2.x_exact(), p2.y_exact(),
			 p3.x_exact(), p3.y_exact());
    }
  return sign;
}

template <class NT>
bool lexicographically_less(NT x1, NT y1,
			    NT x2, NT y2)
{
  CGAL::Sign xx = CGAL::sign(x1 - x2);
  //std::cerr << (x1 - x2) << " -> " << xx  <<"\n";
  if (xx == CGAL::NEGATIVE)
    return true;
  if (xx == CGAL::POSITIVE)
    return false;
  CGAL::Sign yy = CGAL::sign(y1 - y2);
  if (yy == CGAL::NEGATIVE || yy == CGAL::ZERO)
    return true;
  return false;
}

template  <class Constructed_point>
bool lexicographically_less(const Constructed_point &p1, 
			    const Constructed_point &p2)
{
  try
    {
      return lexicographically_less(p1.x_interval(), p1.y_interval(),
				    p2.x_interval(), p2.y_interval());
    } 
  catch(...)
    {
      std::cerr << "here!\n";
      return lexicographically_less(p1.x_exact(), p1.y_exact(),
				    p2.x_exact(), p2.y_exact());
    }
}

class Lexicographically_less
{
public:
  template <class Constructed_point>
  bool operator ()(const Constructed_point &p1, 
		   const Constructed_point &p2)
  {
    return lexicographically_less(p1,p2);
  }
};


template <class Point>
void simple_convex_hull_2(std::vector<Point> &P,
			  std::vector<size_t> &extreme_points)
{
  int n = P.size(), k = 0;
  H.resize(2*n);
  sort(P.begin(), P.end(), Lexicographically_less());
  
  // Build lower hull
  for (int i = 0; i < n; i++) {
    while (k >= 2 && orientation(H[k-2], H[k-1], P[i]) <= 0) k--;
    H[k++] = P[i];
  }
  
  // Build upper hull
  for (int i = n-2, t = k+1; i >= 0; i--) {
    while (k >= t && orientation(H[k-2], H[k-1], P[i]) <= 0) k--;
    H[k++] = P[i];
  }
  
  H.pop_back();
  H.resize(k-1);
}


#endif
