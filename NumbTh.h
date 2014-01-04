/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
 
#ifndef _NumbTh
#define _NumbTh

#include <vector>
#include <cmath>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZX.h>
#include <NTL/GF2X.h>
#include <NTL/vec_ZZ.h>
#include <NTL/xdouble.h>

NTL_CLIENT

// MinGW hack
#ifndef lrand48
#define drand48() (((double)rand()) / RAND_MAX)
#define lrand48() rand()
#define srand48(X) srand(X)
#endif

inline double log2(const xdouble& x){ return log(x) * 1.442695040889; }
inline double log2(const double x){ return log(x) * 1.442695040889; }

// Factoring by trial division, only works for N<2^{60}, only the primes
// are recorded, not their multiplicity. class zz can be long or ZZ
void factorize(vector<long> &factors, long N);
void factorize(vector<ZZ> &factors, const ZZ& N);

/* Compute Phi(N) and also factorize N */
void phiN(long &phiN, vector<long> &facts, long N);
void phiN(ZZ &phiN, vector<ZZ> &facts, const ZZ &N);

/* Compute Phi(N) */
int phi_N(int N);

// Find e-th root of unity modulo the current modulus
void FindPrimitiveRoot(zz_p &r, unsigned e);
void FindPrimitiveRoot(ZZ_p &r, unsigned e);

/* Compute mobius function (naive method as n is small) */
int mobius(int n);

/* Compute cyclotomic polynomial */
ZZX Cyclotomic(int N);

/* Find a primitive root modulo N */
int primroot(int N,int phiN);

int ord(int N,int p);


// Rand mod p poly of degree < n
ZZX RandPoly(int n,const ZZ& p);


/* When abs=false reduce to interval (-q/2,...,q/2), when abs=true reduce
 * to [0,q). When abs=false and q=2, maintains the same sign as the input.
 */
void PolyRed(ZZX& out, const ZZX& in,       int q, bool abs=false);
void PolyRed(ZZX& out, const ZZX& in, const ZZ& q, bool abs=false);
inline void PolyRed(ZZX& F, int q, bool abs=false) { PolyRed(F,F,q,abs); }
inline void PolyRed(ZZX& F, const ZZ& q, bool abs=false)
{ PolyRed(F,F,q,abs); }

GF2X to_GF2X(const ZZX& a);
ZZX to_ZZX(const GF2X& a);



/* Finds whether x is an element of the set X
   of size sz, returns -1 it not and the location if true
*/
int is_in(int x,int* X,int sz);

/* Incremental integer CRT for vectors. Expects co-primes p,q with q odd,
 * and such that all the entries in v1 are in [-p/2,p/2). Returns in v1 the
 * CRT of vp mod p and vq mod q, as integers in [-pq/2, pq/2). Uses the
 * formula:
 *                   CRT(vp,p,vq,q) = vp + [(vq-vp)*p^{-1}]_q * p,
 *
 * where [...]_q means reduction to the interval [-q/2,q/2). Notice that if q
 * is odd then this is the same as reducing to [-(q-1)/2,(q-1)/2], which means
 * that [...]_q * p is in [-p(q-1)/2, p(q-1)/2], and since vp is in [-p/2,p/2)
 * then the sum is indeed in [-pq/2,pq/2).
 *
 * return true is both vectors are of the same length, false otherwise
 */
template <class zzvec>        // zzvec can be vec_ZZ or vec_long
bool intVecCRT(vec_ZZ& vp, const ZZ& p, const zzvec& vq, long q);

// argmax(v)/argmin(v) finds the index of the (first) largest/smallest
//   element in the vector v. They are roughly just simpler variants of
//   std::max_element and std::mim_element.
// argmin/argmax are implemented as a template, so the code must be placed
// in the header file for the comiler to find it. The class T must have an
// implementation of operator> and operator< for this template to work.

template <class T, bool maxFlag>
long argminmax(vector<T>& v)
{
  if (v.size()<1) return -1; // error: this is an empty array
  unsigned idx = 0;
  T target = v[0];
  for (unsigned i=1; i<v.size(); i++)
    if (maxFlag) { if (v[i] > target) { target = v[i]; idx = i;} }
    else         { if (v[i] < target) { target = v[i]; idx = i;} }
  return (long) idx;
}

template <class T> long argmax(vector<T>& v)
{  return argminmax<T,true>(v); }

template <class T> long argmin(vector<T>& v)
{  return argminmax<T,false>(v); }


// Sample polynomials with entries {-1,0,1}. These functions are similar to
// the SampleSmall class from v1, but without a class around it.

// In sampleSmall, each coeff is 0 w/ prob. 1/2 and +-1 w/ prob. 1/4. In
// sampleHWt, min(Hwt,n) random coefficients are chosen at random in {-1,+1}
// and the others are set to zero. If n=0 then n=poly.deg()+1 is used. 

void sampleSmall(ZZX &poly, long n=0);
void sampleHWt(ZZX &poly, long Hwt, long n=0);

// Sample polynomials with Gaussian coefficients. This function is similar
// to the MultivariateGauss class from v1, but without a class around it.

void sampleGaussian(ZZX &poly, long n=0, double stdev=1.0);


// NTL's random number generation faciliity is pretty limited,
// and does not provide a way to save/restore the state of a
// pseudo-random stream.  The following gives us that ability.

class RandomState {
private:
  ZZ state;
  bool restored;

public:
  RandomState() {
    RandomBits(state, 512);
    restored = false;
  }

  void restore() {
    if (!restored) {
      SetSeed(state);
      restored = true;
    }
  }

  ~RandomState() {
    restore();
  }

private:
  RandomState(const RandomState&); // disable copy constructor
  RandomState& operator=(const RandomState&); // disable assignment
};

// usage: 
//  { RandomState state; ... }
//
// The destructor will restore state at end of scope,
// but this can also be done explicitly via the restore method.


// An experimental facility...it is annoying that vector::size()
// is an unsigned quantity...this leads to all kinds of annoying
// warning messages...

template <typename T>
inline long lsize(const vector<T>& v) {
  return (long) v.size();
}


// modular composition: res = g(h) mod f
void ModComp(ZZX& res, const ZZX& g, const ZZX& h, const ZZX& f);

// returns the largest absolute value coefficient
ZZ largestCoeff(const ZZX& f);

#endif
