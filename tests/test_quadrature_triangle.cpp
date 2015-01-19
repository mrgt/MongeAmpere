#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <boost/type_traits.hpp>
#include <MA/quadrature.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_2 Point;


int main()
{
  srand(time(NULL));
  Point a(MA::r01(), MA::r01());
  Point b(MA::r01(), MA::r01());
  Point c(MA::r01(), MA::r01());

  auto f = [](Point p) { return p.y()*p.y()*p.x(); };
  std::cerr << MA::integrate_albrecht_collatz<FT>(a,b,c,f) << "\n";
  std::cerr << MA::integrate_monte_carlo<FT>(a,b,c,f,100000) << "\n";
  std::cerr << MA::integrate_midedge<FT>(a,b,c,f) << "\n";
}
