MongeAmpere++
====================

This library is designed to solve large instance of semi-discrete optimal transport problems, also known as
the constrained least-square assignment problem. More precisely, MongeAmpere++ can be used to solve the quadratic
optimal transport problem between a piecewise linear density on a 2D triangulation and a finite sum of Dirac masses.

Other features that might be useful for other numerical applications include :

* the robust computation of the intersections between a Voronoi diagram (or more generally a power/Laguerre diagram)
  and a triangulation ;
* the computation of the discrete Monge-Ampère operator, and its first and second derivatives with respect to the values
  of the function ; 
* an implementation of Lloyd's algorithm for a piecewise linear density with exact integrals.

The documentation is currently (very) scarse, but you can look at the test/ directory to have a taste of what can be 
done with this MA++.

Dependencies
------------

This software depends on recent versions of CGAL, Boost, SparseSuite (optional), CImg and CMake. It also requires a C++14
compatible compiler.

Tests have been performed with these versions of the following librairies:

* **Gcc** 4.9.1 (C++14 compatibility is required)
+ **CMake** 2.8.12.2
+ **CImg** 1.5.7
+ **Eigen** 3.2.1
+ **Suitesparse** 4.4.1 (optionnal)
+ **Boost** 1.55
+ **Cgal** 4.11
