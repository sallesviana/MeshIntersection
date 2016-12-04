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

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <time.h>
#include <omp.h>
using namespace std;


//===============================================================
// General templates (from WRF source codes)
//===============================================================

template <class T>
T max(const T a, const T b, const T c, const T d) {
  return max(max(a,b), max(c,d));
}

template <class T>
T min(const T a, const T b, const T c, const T d) {
  return min(min(a,b), min(c,d));
}

template <class T, class U>
void accum_max(T &maxsofar, const U x) {
  if (x>maxsofar) maxsofar = x;
}

template <class T>
void accum_min(T &minsofar, const T x) {
  if (x<minsofar) minsofar = x;
}

template <class T>
T square(const T x) { return x*x; }

template <class T>    // Squared length from dx and dy
T sqlen(const T dx, const T dy) { return square(dx)+square(dy); }


//From Boost, for using pairs in unordered_sets:
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
  template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };
}


// Utils for measuring time...

double convertTimeMsecs(const timespec &t);

//source: http://www.guyrutenberg.com/2007/09/22/profiling-code-using-clock_gettime/
timespec diff(timespec start, timespec end);



#endif