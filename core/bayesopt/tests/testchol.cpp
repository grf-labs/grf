/** -*- c++ -*- \file cholesky_test.cpp \brief test cholesky decomposition */
/*
 -   begin                : 2005-08-24
 -   copyright            : (C) 2005 by Gunter Winkler, Konstantin Kutzkow
 -   email                : guwi17@gmx.de

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/

#include <cassert>
#include <limits>

#include <boost/timer.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/numeric/ublas/io.hpp>

#include "ublas_cholesky.hpp"

namespace ublas = boost::numeric::ublas;
namespace bayesopt 
{
  namespace utils
  {
/** \brief make a immutable triangular adaptor from a matrix
 *
 * \usage: 
<code>
 A = triangular< lower >(B);
 A = triangular(B, lower());
</code>
 */
template < class TYPE, class MATRIX >
ublas::triangular_adaptor<const MATRIX, TYPE>
triangular(const MATRIX & A, const TYPE& uplo = TYPE())
{
  return ublas::triangular_adaptor<const MATRIX, TYPE>(A);
}

/** \brief make a immutable banded adaptor from a matrix
 *
 * \usage: 
<code>
 A = banded(B, lower, upper);
</code>
 */
template < class MATRIX >
ublas::banded_adaptor<const MATRIX>
banded(const MATRIX & A, const size_t lower, const size_t upper)
{
  return ublas::banded_adaptor<const MATRIX>(A, lower, upper);
}

/** \brief make a immutable symmetric adaptor from a matrix
 *
 * \usage: 
<code>
 A = symmetric< lower >(B);
 A = symmetric(B, lower());
</code>
 */
template < class TYPE, class MATRIX >
ublas::symmetric_adaptor<const MATRIX, TYPE>
symmetric(const MATRIX & A, const TYPE& uplo = TYPE())
{
  return ublas::symmetric_adaptor<const MATRIX, TYPE>(A);
}


/** \brief fill lower triangular matrix L 
 */
template < class MATRIX >
void fill_symm(MATRIX & L, const size_t bands = std::numeric_limits<size_t>::max() )
{
  typedef typename MATRIX::size_type size_type;
  
  assert(L.size1() == L.size2());

  size_type size = L.size1();
  for (size_type i=0; i<size; i++) {
    for (size_type j = ((i>bands)?(i-bands):0); j<i; j++) {
      L(i,j) = 1 + (1.0 + i)/(1.0 + j) + 1.0 / (0.5 + i - j);
    }
    L(i,i) = 1+i+size;
  }

  return;
}

int main_(int argc, char * argv[] )
{
  size_t size = 10;
  if (argc > 1)
    size = ::atoi (argv [1]);
  boost::timer  t1;
  double pr, de, sv;

  typedef double DBL;
  typedef ublas::row_major  ORI;
  {
    // use dense matrix
    ublas::matrix<DBL, ORI> A (size, size);
    ublas::matrix<DBL, ORI> T (size, size);
    ublas::matrix<DBL, ORI> L (size, size);

    A = ublas::zero_matrix<DBL>(size, size);

    ublas::vector<DBL> b (size);
    ublas::vector<DBL> x (size);
    ublas::vector<DBL> y (size);

    std::fill(y.begin(), y.end(), 1.0);

    fill_symm(T);
    t1.restart();
    A = ublas::prod(T, trans(T));
    pr = t1.elapsed();

    b = prod( A, y);

    t1.restart();
    size_t res = cholesky_decompose(A, L);
    de = t1.elapsed();

    t1.restart();
    x = b;
    cholesky_solve(L, x, ublas::lower());
    sv = t1.elapsed();

    std::cout << res << ": " 
              << ublas::norm_inf(L-T) << " "
              << ublas::norm_2(x-y) << " "
              << " (deco: " << de << " sec)"
              << " (prod: " << pr << " sec)"
              << " (solve: " << sv << " sec)"
              << " " << size
              << std::endl;
  }

  {
    // use dense triangular matrices
    ublas::triangular_matrix<DBL, ublas::lower, ORI> A (size, size);
    ublas::triangular_matrix<DBL, ublas::lower, ORI> T (size, size);
    ublas::triangular_matrix<DBL, ublas::lower, ORI> L (size, size);
    
    A = ublas::zero_matrix<DBL> (size, size) ;
    A = triangular<ublas::lower>( ublas::zero_matrix<DBL> (size, size) );
    A = triangular( ublas::zero_matrix<DBL> (size, size), ublas::lower() );

    ublas::vector<DBL> b (size);
    ublas::vector<DBL> x (size);
    ublas::vector<DBL> y (size);

    std::fill(y.begin(), y.end(), 1.0);

    fill_symm(T);
    t1.restart();
    A = triangular<ublas::lower>( ublas::prod(T, trans(T)) );
    pr = t1.elapsed();

    b = prod( symmetric<ublas::lower>(A), y);

    t1.restart();
    size_t res = cholesky_decompose(A, L);
    de = t1.elapsed();

    t1.restart();
    x = b;
    cholesky_solve(L, x, ublas::lower());
    sv = t1.elapsed();

    std::cout << res << ": " 
              << ublas::norm_inf(L-T) << " "
              << ublas::norm_2(x-y) << " "
              << " (deco: " << de << " sec)"
              << " (prod: " << pr << " sec)"
              << " (solve: " << sv << " sec)"
              << " " << size
              << std::endl;
//     std::cout << L << std::endl;
//     std::cout << b << std::endl;
//     std::cout << x << std::endl;
//     std::cout << y << std::endl;
//     std::cout << (b - prod(symmetric<ublas::lower>(A), x)) << std::endl;
  }

  {
    // use banded matrices matrix
    typedef ublas::banded_matrix<DBL, ORI> MAT;

    size_t bands = std::min<size_t>( size, 50 );
    MAT A (size, size, bands, 0);
    MAT T (size, size, bands, 0);
    MAT L (size, size, bands, 0);
    
    A = ublas::zero_matrix<DBL> (size, size) ;
    A = banded( ublas::zero_matrix<DBL> (size, size), bands, 0 );

    ublas::vector<DBL> b (size);
    ublas::vector<DBL> x (size);
    ublas::vector<DBL> y (size);

    std::fill(y.begin(), y.end(), 1.0);

    fill_symm(T, bands);
    t1.restart();
    A = banded( ublas::prod(T, trans(T)), bands, 0 );
    pr = t1.elapsed();

    b = prod( symmetric<ublas::lower>(A), y);

    t1.restart();
    size_t res = cholesky_decompose(A, L);
    de = t1.elapsed();

    t1.restart();
    x = b;
    cholesky_solve(L, x, ublas::lower());
    sv = t1.elapsed();

    std::cout << res << ": " 
              << ublas::norm_inf(L-T) << " "
              << ublas::norm_2(x-y) << " "
              << " (deco: " << de << " sec)"
              << " (prod: " << pr << " sec)"
              << " (solve: " << sv << " sec)"
              << " " << size
              << std::endl;
  }

  return EXIT_SUCCESS;
}

  }
}

int main(int argc, char * argv[] )
{ return bayesopt::utils::main_(argc,argv); }

/****************

(note: dense + tria + 50 banded)

$ g++-4.0 -I $HOME/include -o cholesky_test cholesky_test.cpp -Wall -g -O2 -DNDEBUG

(column major, double)

0: 3.05533e-13 1.78207e-14  (deco: 0 sec) (prod: 0 sec) (solve: 0 sec) 100
0: 3.05533e-13 1.78207e-14  (deco: 0 sec) (prod: 0.01 sec) (solve: 0 sec) 100
0: 1.97176e-13 9.0801e-15  (deco: 0 sec) (prod: 0 sec) (solve: 0 sec) 100
0: 1.18172e-12 5.82623e-14  (deco: 0 sec) (prod: 0.06 sec) (solve: 0 sec) 200
0: 1.18172e-12 5.82623e-14  (deco: 0.03 sec) (prod: 0.03 sec) (solve: 0 sec) 200
0: 1.15463e-13 1.01481e-14  (deco: 0.01 sec) (prod: 0.01 sec) (solve: 0 sec) 200
0: 4.80105e-12 2.13833e-13  (deco: 0.1 sec) (prod: 0.58 sec) (solve: 0 sec) 400
0: 4.80105e-12 2.13833e-13  (deco: 0.22 sec) (prod: 0.22 sec) (solve: 0.01 sec) 400
0: 6.39488e-14 1.14428e-14  (deco: 0.02 sec) (prod: 0.02 sec) (solve: 0.01 sec) 400
0: 1.80402e-11 3.72721e-13  (deco: 1.27 sec) (prod: 9.87 sec) (solve: 0 sec) 800
0: 1.80402e-11 3.72721e-13  (deco: 2.16 sec) (prod: 2.13 sec) (solve: 0.04 sec) 800
0: 4.04121e-14 1.52388e-14  (deco: 0.05 sec) (prod: 0.05 sec) (solve: 0.01 sec) 800
0: 7.85829e-11 1.51858e-12  (deco: 33.5 sec) (prod: 279.15 sec) (solve: 0.05 sec) 1600
0: 7.85829e-11 1.51858e-12  (deco: 21.8 sec) (prod: 22.45 sec) (solve: 0.18 sec) 1600
0: 2.26485e-14 1.98783e-14  (deco: 0.19 sec) (prod: 0.09 sec) (solve: 0.01 sec) 1600

(row major, double)

0: 3.05533e-13 1.78207e-14  (deco: 0 sec) (prod: 0 sec) (solve: 0 sec) 100
0: 3.05533e-13 1.78207e-14  (deco: 0 sec) (prod: 0 sec) (solve: 0 sec) 100
0: 1.97176e-13 9.0801e-15  (deco: 0.01 sec) (prod: 0 sec) (solve: 0 sec) 100
0: 1.18172e-12 5.82623e-14  (deco: 0.01 sec) (prod: 0.02 sec) (solve: 0 sec) 200
0: 1.18172e-12 5.82623e-14  (deco: 0.02 sec) (prod: 0.02 sec) (solve: 0 sec) 200
0: 1.15463e-13 1.01481e-14  (deco: 0.01 sec) (prod: 0.01 sec) (solve: 0 sec) 200
0: 4.80105e-12 2.13833e-13  (deco: 0.03 sec) (prod: 0.21 sec) (solve: 0.01 sec) 400
0: 4.80105e-12 2.13833e-13  (deco: 0.08 sec) (prod: 0.09 sec) (solve: 0.01 sec) 400
0: 6.39488e-14 1.14428e-14  (deco: 0.02 sec) (prod: 0.02 sec) (solve: 0 sec) 400
0: 1.80402e-11 3.72721e-13  (deco: 0.37 sec) (prod: 1.51 sec) (solve: 0 sec) 800
0: 1.80402e-11 3.72721e-13  (deco: 0.69 sec) (prod: 0.57 sec) (solve: 0.03 sec) 800
0: 4.04121e-14 1.52388e-14  (deco: 0.06 sec) (prod: 0.04 sec) (solve: 0 sec) 800
0: 7.85829e-11 1.51858e-12  (deco: 2.57 sec) (prod: 11.26 sec) (solve: 0.05 sec) 1600
0: 7.85829e-11 1.51858e-12  (deco: 4.98 sec) (prod: 4.43 sec) (solve: 0.15 sec) 1600
0: 2.26485e-14 1.98783e-14  (deco: 0.21 sec) (prod: 0.06 sec) (solve: 0.01 sec) 1600

$ g++-3.3.5 -I $HOME/include -o cholesky_test cholesky_test.cpp -Wall -g -O2 -DNDEBUG

0: 6.72351e-13 7.10356e-14  (deco: 0.01 sec) (prod: 0.02 sec) (solve: 0 sec) 100
0: 6.72351e-13 7.10356e-14  (deco: 0.01 sec) (prod: 0.01 sec) (solve: 0 sec) 100
0: 1.12355e-13 1.96079e-14  (deco: 0.01 sec) (prod: 0 sec) (solve: 0 sec) 100
0: 2.61746e-12 2.75266e-13  (deco: 0.02 sec) (prod: 0.15 sec) (solve: 0 sec) 200
0: 2.61746e-12 2.75266e-13  (deco: 0.05 sec) (prod: 0.07 sec) (solve: 0 sec) 200
0: 7.32747e-14 1.38831e-14  (deco: 0.02 sec) (prod: 0.02 sec) (solve: 0 sec) 200
0: 1.02713e-11 1.33167e-12  (deco: 0.18 sec) (prod: 1.04 sec) (solve: 0 sec) 400
0: 1.02713e-11 1.33167e-12  (deco: 0.39 sec) (prod: 0.4 sec) (solve: 0.01 sec) 400
0: 3.33067e-14 1.4327e-14  (deco: 0.03 sec) (prod: 0.04 sec) (solve: 0 sec) 400
0: 4.21783e-11 3.23897e-12  (deco: 1.55 sec) (prod: 8.1 sec) (solve: 0.01 sec) 800
0: 4.21783e-11 3.23897e-12  (deco: 3.17 sec) (prod: 3.07 sec) (solve: 0.02 sec) 800
0: 1.42109e-14 1.87039e-14  (deco: 0.09 sec) (prod: 0.06 sec) (solve: 0 sec) 800
0: 1.6946e-10 9.93422e-12  (deco: 12.09 sec) (prod: 64.96 sec) (solve: 0.07 sec) 1600
0: 1.6946e-10 9.93422e-12  (deco: 24.72 sec) (prod: 24.28 sec) (solve: 0.1 sec) 1600
0: 1.15463e-14 2.46404e-14  (deco: 0.28 sec) (prod: 0.1 sec) (solve: 0 sec) 1600


 ****************/
