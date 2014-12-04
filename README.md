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

This software depends on recent versions of CGAL, Boost, SparseSuite (optional), CImg and CMake. It also requires a C++11
compatible compiler.

Tests have been lead with these versions of the following librairies:

* **Gcc** 4.9.1 (C++11 compatibility)
* **CMake** 2.8.12.2
+ **Eigen** 3.2.1
+ **Suitesparse** 4.4.1 (optionnal)
+ **Boost** 1.55
+ **Cgal** 4.4

``` sh
sudo apt-get install cmake-curses-gui
sudo apt-get install libeigen3-dev
sudo apt-get install cimg-dev
sudo apt-get install libboost-all-dev
sudo apt-get install libcgal-dev
```

Installation
------------

Create a temporary folder

``` sh
cd path/to/MongeAmpere/..
mkdir buildMA
cd buildMA
```

Build with ccmake:

``` sh
ccmake ../MongeAmpere
```

Then type successively c, c, g in ccmake to configure MongeAmpere. Then build it:

``` sh
make
```

# For old systems

For older linux versions or systems where you are not root, the sources may not be up to date. You can install the dependences localy via linuxbrew, which has is own sources.

Dependences of linuxbrew
------------------------

``` sh
sudo apt-get install build-essential curl git m4 ruby texinfo libbz2-dev libcurl4-openssl-dev libexpat-dev libncurses-dev zlib1g-dev
```

Install of linuxbrew
--------------------

``` sh
ruby -e "$(wget -O- https://raw.github.com/Homebrew/linuxbrew/go/install)"
```

Environment variables
---------------------

Add these lines at the end of your ~/.bashrc, to adapt your ENV variables

``` sh
export PATH="$HOME/.linuxbrew/bin:$PATH"
export MANPATH="$HOME/.linuxbrew/share/man:$MANPATH"
export INFOPATH="$HOME/.linuxbrew/share/info:$INFOPATH"
```

Errors and updates
------------------

``` sh
brew doctor
brew update
brew install cgal
```

Library installation
--------------------

``` sh
brew install <libraryname>
```

Finally install the program as described higher

# For old systems with gcc < 4.8.1

To use the linuxbrew compiler chain instead of the old one installed on your machine, check the instructions on the [Standalone-Installation](https://github.com/Homebrew/linuxbrew/wiki/Standalone-Installation) website.
For more information about how to use linuxbrew and how to correct possible errors, check [linuxbrew](https://github.com/Homebrew/linuxbrew/) website.
