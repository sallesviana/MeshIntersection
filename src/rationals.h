/*
Copyright 2016 Salles V. G. Magalhaes, W. R. Franklin, Marcus Andrade

This file is part of PinMesh.

PinMesh is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PinMesh is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PinMesh.  If not, see <http://www.gnu.org/licenses/>.
*/


//Some wrappers for rational and big numbers
//

#ifndef RATIONALS_H
#define RATIONALS_H

#define USE_GMPXX
#undef USE_BOOST

//#define USE_BOOST
//#undef USE_GMPXX

#ifdef USE_GMPXX
#include <gmpxx.h>
typedef  mpq_class rational;
typedef  mpz_class big_int;
#endif

#ifdef USE_BOOST
#include <boost/multiprecision/cpp_int.hpp>
using namespace boost::multiprecision;
typedef  cpp_rational rational;
typedef  cpp_int big_int;
#endif

inline const double castDoubleWrap(const rational &q) {
#ifdef USE_GMPXX
	return q.get_d();
#endif
#ifdef USE_BOOST
	return (double)q;
#endif
}

/*
inline const Area castAreaWrap(const rational &q) {
#ifdef USE_GMPXX
	return q.get_d();
#endif
#ifdef USE_BOOST
	return q.convert_to<Area>();
#endif
}
*/

inline const int castIntWrap(const big_int &z) {
#ifdef USE_GMPXX
	return (int)z.get_si();
#endif
#ifdef USE_BOOST
	return int(z);
#endif
}


inline const big_int numeratorWrap(const rational &q) {
#ifdef USE_GMPXX
	return q.get_num();
#endif
#ifdef USE_BOOST
	return numerator(q);
#endif
}

inline const big_int denominatorWrap(const rational &q) {
#ifdef USE_GMPXX
	return q.get_den();
#endif
#ifdef USE_BOOST
	return denominator(q);
#endif
}

inline const int convertToInt(const rational &q, big_int tempVars[]) {
	tempVars[0] = q.get_num();
	tempVars[1] = q.get_den();
	tempVars[0] /= tempVars[1];
	return castIntWrap(tempVars[0]);
	//return q.get_d();
}

inline void divide_qrWrap(const big_int &n,const big_int &d, big_int &i,big_int &r) {
#ifdef USE_GMPXX
	i = n/d;
	r = n%d;
#endif
#ifdef USE_BOOST
	divide_qr(n,d,i,r);
#endif
}


#endif