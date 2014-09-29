
// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Modified by Q. Merigot: removed expression templates for
// compatibility with CGAL

#ifndef AUTODIFF_SCALAR_H
#define AUTODIFF_SCALAR_H
#include <Eigen/Sparse>

  class AD
  {
    typedef double Scalar;
    typedef Eigen::SparseVector<double> Vector;
    Scalar m_value;
    Vector m_derivatives;
    
  public:
    /** Default constructor without any initialization. */
    AD() {}
    
    /** Constructs an active scalar from its \a value,
	and initializes the \a nbDer derivatives such that it corresponds to the \a derNumber -th variable */
    AD(const Scalar& value, int nbDer, int derNumber)
      : m_value(value), m_derivatives(Vector(nbDer))
    {
      m_derivatives.coeffRef(derNumber) = Scalar(1);
    }
    
    /** Conversion from a scalar constant to an active scalar.
     * The derivatives are set to zero. */
    /* explicit */ AD(const Scalar& value)
      : m_value(value)
    {
      if(m_derivatives.size()>0)
	m_derivatives.setZero();
    }
    
    /** Constructs an active scalar from its \a value and derivatives \a der */
    AD(const Scalar& value, const Vector& der)
      : m_value(value), m_derivatives(der)
    {}
    
    friend inline
    std::ostream & operator << (std::ostream & s, const AD& a)
    {
      return s << a.value();
    }
    
    AD(const AD& other)
      : m_value(other.value()), m_derivatives(other.derivatives())
    {}
    
    inline AD& operator=(const AD& other)
    {
      m_value = other.value();
      m_derivatives = other.derivatives();
      return *this;
    }
    
    inline const Scalar& value() const { return m_value; }
    inline Scalar& value() { return m_value; }
    
    inline const Vector& derivatives() const { return m_derivatives; }
    inline Vector& derivatives() { return m_derivatives; }
    
    inline bool operator< (const Scalar& other) const 
    { return m_value <  other; }
    inline bool operator< (const AD& other) const  
    { return m_value <  other.value(); }
    inline bool operator<=(const Scalar& other) const  
    { return m_value <= other; }
    inline bool operator<= (const AD& other) const  
    { return m_value <= other.value(); }
    inline bool operator> (const Scalar& other) const  
    { return m_value >  other; }
    inline bool operator> (const AD& other) const  
    { return m_value >  other.value(); }
    inline bool operator>=(const Scalar& other) const  
    { return m_value >= other; }
    inline bool operator>= (const AD& other) const  
    { return m_value >= other.value(); }
    inline bool operator==(const Scalar& other) const  
    { return m_value == other; }
    inline bool operator==(const AD& other) const  
    { return m_value == other; }
    inline bool operator!=(const Scalar& other) const  
    { return m_value != other; }
    inline bool operator!=(const AD& other) const  
    { return m_value == other; }
    
    friend inline bool operator< (const Scalar& a, const AD& b) 
    { return a <  b.value(); }
    friend inline bool operator<=(const Scalar& a, const AD& b) 
    { return a <= b.value(); }
    friend inline bool operator> (const Scalar& a, const AD& b) 
    { return a >  b.value(); }
    friend inline bool operator>=(const Scalar& a, const AD& b) 
    { return a >= b.value(); }
    friend inline bool operator==(const Scalar& a, const AD& b) 
    { return a == b.value(); }
    friend inline bool operator!=(const Scalar& a, const AD& b) 
    { return a != b.value(); }
    
    inline AD operator+(const Scalar& other) const
    {
      return AD(m_value + other, m_derivatives);
    }
    
    friend inline AD operator+(const Scalar& a, const AD& b)
    {
      return AD(a + b.value(), b.derivatives());
    }

    inline AD& operator+=(const Scalar& other)
    {
      value() += other;
      return *this;
    }

    AD operator+(const AD& other) const
    {
      return AD(m_value + other.value(),
		m_derivatives + other.derivatives());
    }
    
    AD& operator+=(const AD& other)
    {
      (*this) = (*this) + other;
      return *this;
    }
    
    AD operator-(const Scalar& b) const
    {
      return AD(m_value - b, m_derivatives);
    }
    
    friend inline AD operator-(const Scalar& a, const AD& b)
    {
      return AD(a - b.value(), -b.derivatives());
    }
    
    AD& operator-=(const Scalar& other)
    {
      value() -= other;
      return *this;
    }
    
    AD operator-(const AD& other) const
    {
      return AD(m_value - other.value(),
		m_derivatives - other.derivatives());
    }
    
    AD& operator-=(const AD& other)
    {
      *this = *this - other;
      return *this;
    }
    
    AD operator-() const
    {
      return AD(-m_value, -m_derivatives);
    }
    
    AD operator*(const Scalar& other) const
    {
      return AD(m_value * other, (m_derivatives * other));
    }
    
    friend inline AD operator*(const Scalar& other, const AD& a)
    {
      return AD(a.value() * other,
		a.derivatives() * other);
    }
    
    AD operator/(const Scalar& other) const
    {
      return AD(m_value / other, (m_derivatives * (Scalar(1)/other)));
    }
    
    friend inline AD operator/(const Scalar& other, const AD& a)
    {
      return AD(other / a.value(),
		a.derivatives() *
		(Scalar(-other) / (a.value()*a.value())));
    }
    
    AD operator/(const AD& other) const
    {
      return AD(m_value / other.value(),
		((m_derivatives * other.value()) - 
	       (m_value * other.derivatives())) * 
		(Scalar(1)/(other.value()*other.value())));
    }
    
    AD operator*(const AD& other) const
    {
      return AD(m_value * other.value(),
		(m_derivatives * other.value()) +
		(m_value * other.derivatives()));
    }
    
    inline AD& operator*=(const Scalar& other)
    {
      *this = *this * other;
      return *this;
    }
    
    inline AD& operator*=(const AD& other)
    {
      *this = *this * other;
      return *this;
    }
    
    inline AD& operator/=(const Scalar& other)
    {
      *this = *this / other;
      return *this;
    }
    
    inline AD& operator/=(const AD& other)
    {
      *this = *this / other;
      return *this;
    }
  };
  
  AD sqrt(const AD &x)
  {
    using std::sqrt;
    double sqrtx = sqrt(x.value());
    return AD(sqrtx,x.derivatives() * (double(0.5) / sqrtx));
  }

namespace CGAL
{
  template <> class Algebraic_structure_traits<AD>
    : public Algebraic_structure_traits_base<AD, Field_tag >  
  {
  public:
    typedef Tag_false            Is_exact;
    typedef Tag_true             Is_numerical_sensitive;
  };

  template <>
  class Real_embeddable_traits< AD >
    : public INTERN_RET::Real_embeddable_traits_base< AD , CGAL::Tag_true>
  {
  public:
    class To_double : public std::unary_function<AD, double >
    {
    public:
      double operator()( const AD& x ) const {
    	return x.value();
      }
    };
    class To_interval 
      : public std::unary_function< AD, std::pair<double,double> > {     
    public:
      std::pair<double,double> operator()( const AD& x ) const {
	double dx(x.value());
	return std::make_pair(dx,dx);
      }
    };

  };

  template <>
  struct NT_converter <AD, Gmpq>
    : public std::unary_function<AD, Gmpq>
  {
    Gmpq
    operator()(const AD &a) const
    {
      return Gmpq(a.value());
    }
  };

}

  
#endif // EIGEN_AUTODIFF_SCALAR_H
