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

/* DoubleCRT.cpp - This class holds an integer polynomial in double-CRT form
 *
 * Double-CRT form is a matrix of L rows and phi(m) columns. The i'th row
 * contains the FFT of the element wrt the ith prime, i.e. the evaluations of
 * the polynomial at the primitive mth roots of unity mod the ith prime. The
 * polynomial thus represented is defined modulo the product of all the primes
 * in use. The list of primes is defined by the data member modChain, which is
 * a vector of Cmodulus objects. 
 */
#include <NTL/ZZX.h>
NTL_CLIENT

#include "NumbTh.h"
#include "PAlgebra.h"
#include "CModulus.h"

#include "DoubleCRT.h"
#include "SingleCRT.h"

// NTL implementation of mat_long

//NTL_matrix_impl(long,vec_long,vec_vec_long,mat_long)
//NTL_io_matrix_impl(long,vec_long,vec_vec_long,mat_long)
//NTL_eq_matrix_impl(long,vec_long,vec_vec_long,mat_long)


// representing an integer polynomial as DoubleCRT. If the number of moduli
// to use is not specified, the resulting object uses all the moduli in
// the context. If the coefficients of poly are larger than the product of
// the used moduli, they are effectively reduced modulo that product


// a "sanity check" function, verifies consistency of matrix with current
// moduli chain an error is raised if they are not consistent
void DoubleCRT::verify()
{
  const IndexSet& s = map.getIndexSet();
  if (s.last() >= context.numPrimes())
    Error("DoubleCRT object has too many rows");

  long phim = context.zMstar.phiM();

  // check that the content of i'th row is in [0,pi) for all i
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];

    if (row.length() != phim) 
      Error("DoubleCRT object has bad row length");

    long pi = context.ithPrime(i); // the i'th modulus
    for (long j=0; j<phim; j++)
      if (row[j]<0 || row[j]>= pi) 
	Error("DoubleCRT object has inconsistent data");
  }
}



// Arithmetic operations. Only the "destructive" versions are used,
// i.e., a += b is implemented but not a + b.

// Generic operation, Fnc is AddMod, SubMod, or MulMod (from NTL's ZZ module)
DoubleCRT& DoubleCRT::Op(const DoubleCRT &other, long (*Fnc)(long, long, long),
			 bool matchIndexSets)
{
  if (&context != &other.context)
    Error("DoubleCRT::Op: incompatible objects");

  // Match the index sets, if needed
  if (matchIndexSets && !(map.getIndexSet() >= other.map.getIndexSet()))
    addPrimes(other.map.getIndexSet() / map.getIndexSet());  // This is expensive

  // INVARIANT: map.getIndexSet() >= other.map.getIndexSet()) 

  // If you need to mod-up the other, do it on a temporary scratch copy
  DoubleCRT tmp(context); 
  const IndexMap<vec_long>* other_map = &other.map;
  if (map.getIndexSet() > other.map.getIndexSet()) { // Even more expensive
    tmp = other;
    tmp.addPrimes(map.getIndexSet() / other.map.getIndexSet());
    other_map = &tmp.map;
  }

  const IndexSet& s = map.getIndexSet();
  long phim = context.zMstar.phiM();

  // add/sub/mul the data, element by element, modulo the respective primes
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    vec_long& row = map[i];
    const vec_long& other_row = (*other_map)[i];
    
    for (long j = 0; j < phim; j++)
      row[j] = Fnc(row[j], other_row[j], pi);
  }
  return *this;
}

DoubleCRT& DoubleCRT::Op(const ZZ &num, long (*Fnc)(long, long, long))
{

  const IndexSet& s = map.getIndexSet();
  long phim = context.zMstar.phiM();
  
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    long n = rem(num, pi);  // n = num % pi
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++)
      row[j] = Fnc(row[j], n, pi);
  }
  return *this;
}

DoubleCRT& DoubleCRT::Op(const ZZX &poly, long (*Fnc)(long, long, long))
{
  const IndexSet& s = map.getIndexSet();
  DoubleCRT other(poly, context, s); // other defined wrt same primes as *this

  return Op(other, Fnc);
}


// expand index set by s1.
// it is assumed that s1 is disjoint from the current index set.
void DoubleCRT::addPrimes(const IndexSet& s1)
{
  if (card(s1) == 0) return; // nothing to do
  assert( disjoint(s1,map.getIndexSet()) ); // s1 is disjoint from *this

  ZZX poly;
  toPoly(poly); // recover in coefficient representation

  map.insert(s1);  // add new rows to the map

  // fill in new rows
  for (long i = s1.first(); i <= s1.last(); i = s1.next(i)) {
    context.ithModulus(i).FFT(map[i], poly); // reduce mod p_i and store FFT image
  }
}


// expand index set by s1, and for every q in s1 scale by q*(q^{-1} mod p).
// it is assumed that s1 is disjoint from the current index set.
// Returns the logarithm of the product of the added factors.
double DoubleCRT::addPrimesAndScale(const IndexSet& s1)
{
  const ZZ p = context.ModulusP();
  if (card(s1) == 0) return 0.0; // nothing to do
  assert(p>=2);
  assert(card(s1 & map.getIndexSet()) == 0); // s1 is disjoint from *this

  // compute factor to scale existing rows
  ZZ factor = to_ZZ(1);
  double logFactor = 0.0;
  for (long i = s1.first(); i <= s1.last(); i = s1.next(i)) {
    long qi = context.ithPrime(i);
    factor *= qi;
    logFactor += log((double)qi);
  }

  //  multiply factor by factor^{-1} mod p
  ZZ tmp;
  rem(tmp, factor, p);

  ZZ prodInv = InvMod(tmp,p);
  factor *= prodInv;
  logFactor += log(prodInv);

  long phim = context.zMstar.phiM();

  // scale existing rows
  const IndexSet& iSet = map.getIndexSet();
  for (long i = iSet.first(); i <= iSet.last(); i = iSet.next(i)) {
    long qi = context.ithPrime(i);
    long f = rem(factor, qi);     // f = factor % qi
    vec_long& row = map[i];
    // scale row by a factor of f modulo qi
    mulmod_precon_t bninv = PrepMulModPrecon(f, qi);
    for (long j=0; j<phim; j++) 
      row[j] = MulModPrecon(row[j], f, qi, bninv);
  }

  // insert new rows and fill them with zeros
  map.insert(s1);  // add new rows to the map
  for (long i = s1.first(); i <= s1.last(); i = s1.next(i)) {
    vec_long& row = map[i];
    for (long j=0; j<phim; j++) row[j] = 0;
  }

  return logFactor;
}



DoubleCRT::DoubleCRT(const ZZX& poly, const FHEcontext &_context, const IndexSet& s)
: context(_context), map(new DoubleCRTHelper(_context))
{
  assert(s.last() < context.numPrimes());

  map.insert(s);

  // convert the integer polynomial to FFT representation modulo the primes
  for (long i = 0; i <= s.last(); i = s.next(i)) {
    const Cmodulus &pi = context.ithModulus(i);
    pi.FFT(map[i], poly); // reduce mod pi and store FFT image
  }
}


DoubleCRT::DoubleCRT(const ZZX& poly, const FHEcontext &_context)
: context(_context), map(new DoubleCRTHelper(_context))
{
  IndexSet s = context.ctxtPrimes;//IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);

  // convert the integer polynomial to FFT representation modulo the primes
  for (long i = 0; i <= s.last(); i = s.next(i)) {
    const Cmodulus &pi = context.ithModulus(i);
    pi.FFT(map[i], poly); // reduce mod pi and store FFT image
  }
}



DoubleCRT::DoubleCRT(const ZZX& poly)
: context(*activeContext), map(new DoubleCRTHelper(*activeContext))
{
  IndexSet s = context.ctxtPrimes;//IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);

  // convert the integer polynomial to FFT representation modulo the primes
  for (long i = 0; i <= s.last(); i = s.next(i)) {
    const Cmodulus &pi = context.ithModulus(i);
    pi.FFT(map[i], poly); // reduce mod pi and store FFT image
  }
}



DoubleCRT::DoubleCRT(const FHEcontext &_context, const IndexSet& s)
: context(_context), map(new DoubleCRTHelper(_context))
{
  assert(s.last() < context.numPrimes());

  map.insert(s);
  long phim = context.zMstar.phiM();

  for (long i = 0; i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++) row[j] = 0;
  }
}


DoubleCRT::DoubleCRT(const FHEcontext &_context)
: context(_context), map(new DoubleCRTHelper(_context))
{
  IndexSet s = context.ctxtPrimes;//IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);

  long phim = context.zMstar.phiM();

  for (long i = 0; i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++) row[j] = 0;
  }
}

DoubleCRT::DoubleCRT(const DoubleCRT& other)
: context(other.context), map(other.map) { }
// NOTE: this is just the standard copy constructor...maybe
// remove the declaration and definition

DoubleCRT::DoubleCRT()
: context(*activeContext), map(new DoubleCRTHelper(*activeContext))
{
  IndexSet s = context.ctxtPrimes;//IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);

  long phim = context.zMstar.phiM();

  for (long i = 0; i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++) row[j] = 0;
  }
}

DoubleCRT& DoubleCRT::operator=(const DoubleCRT& other)
{
   if (&context != &other.context) 
      Error("DoubleCRT assigment: incompatible contexts");

   map = other.map;
   return *this;
}


DoubleCRT& DoubleCRT::operator=(const ZZX&poly)
{
  const IndexSet& s = map.getIndexSet();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) 
    context.ithModulus(i).FFT(map[i], poly); // reduce mod pi and store FFT image

  return *this;
}

DoubleCRT& DoubleCRT::operator=(const ZZ& num)
{
  const IndexSet& s = map.getIndexSet();
  long phim = context.zMstar.phiM();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];
    long pi = context.ithPrime(i);
    long n = rem(num, pi);

    for (long j = 0; j < phim; j++) row[j] = n;
  }

  return *this;
}

void DoubleCRT::toPoly(ZZX& poly, const IndexSet& s,
		       bool positive) const
{
  IndexSet s1 = map.getIndexSet() & s;

  if (card(s1) == 0) {
    clear(poly);
    return;
  }

  ZZ p = to_ZZ(context.ithPrime(s1.first()));  // the first modulus

  // Get poly modulo the first prime in coefficent form
  long i = s1.first();
  const Cmodulus& mod = context.ithModulus(i);
  mod.iFFT(poly, map[i]);

  vec_ZZ& vp = poly.rep;

  // ensure that vp is of size phi(m) with entries in [-p/2,p/2]
  long phim = context.zMstar.phiM();
  long vpLength = vp.length();
  if (vpLength < phim) { // just in case of leading zeros in poly
    vp.SetLength(phim);
    for (long j = vpLength; j < phim; j++) vp[j]=0;
  }
  ZZ p_over_2 = p/2;
  for (long j = 0; j < phim; j++) if (vp[j] > p_over_2) vp[j] -= p;

  // do incremental integer CRT for other levels

  ZZX current;
  for (i = s1.next(i); i <= s1.last(); i = s1.next(i)) {
    long q = context.ithPrime(i);          // the next modulus
    context.ithModulus(i).iFFT(current, map[i]); // Poly mod q in coeff form

    // CRT the coefficient vectors of poly and current
    intVecCRT(vp, p, current.rep, q);    // defined in the module NumbTh
    p *= q;     // update the modulus
  }

  // The above yeilds polynomial with coefficients in [-p/2,p/2]
  // If we need positive, just add p to all the negative coefficients
  if (positive) 
    for (long j=0; j<poly.rep.length(); j++) {
      if (poly.rep[j] < 0) poly.rep[j] += p;
    }

  poly.normalize(); // need to call this after we work on the coeffs
}

void DoubleCRT::toPoly(ZZX& p, bool positive) const
{
  const IndexSet& s = map.getIndexSet();
  toPoly(p, s, positive);
}

// Division by constant
DoubleCRT& DoubleCRT::operator/=(const ZZ &num)
{
  const IndexSet& s = map.getIndexSet();
  long phim = context.zMstar.phiM();
  
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    long n = InvMod(rem(num, pi),pi);  // n = num^{-1} mod pi
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++)
      row[j] = MulMod(row[j], n, pi);
  }
  return *this;
}

// Small-exponent polynomial exponentiation
void DoubleCRT::Exp(long e)
{
  const IndexSet& s = map.getIndexSet();
  long phim = context.zMstar.phiM();
  
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++)
      row[j] = PowerMod(row[j], e, pi);
  }
}



// Apply the automorphism F(X) --> F(X^k)  (with gcd(k,m)=1)
void DoubleCRT::automorph(long k)
{
  const PAlgebra& zmStar = context.zMstar;
  if (!zmStar.inZmStar(k))
    Error("DoubleCRT::automorph: k not in Zm*");

  long m = zmStar.M();
  vector<long> tmp(m);  // temporary array of size m

  const IndexSet& s = map.getIndexSet();

  // go over the rows, permute them one at a time
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];

    for (long j=1; j<m; j++) { // 1st pass: copy to temporary array
      long idx = zmStar.indexInZmstar(j); // returns -1 if j \notin (Z/mZ)*
      if (idx>=0) tmp[j] = row[idx];
    }

    for (long j=1; j<m; j++) { // 2nd pass: copy back from temporary array
      long idx = zmStar.indexInZmstar(j); // returns -1 if j \notin (Z/mZ)*
      if (idx>=0) row[idx] = tmp[MulMod(j,k,m)];
                                           // new[j] = old[j*k mod m]
    }    
  }
}

// fills ith row with random integers mod pi
void DoubleCRT::randomize(const ZZ* seed) 
{
  if (seed != NULL) SetSeed(*seed);

  const IndexSet& s = map.getIndexSet();
  long phim = context.zMstar.phiM();
  
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++)
      row[j] = RandomBnd(pi);   // RandomBnd is defined in NTL's module ZZ
  }
}


DoubleCRT& DoubleCRT::operator=(const SingleCRT& scrt)
{
  if (&context != &scrt.getContext())
    Error("DoubleCRT=SingleCRT -- incompatible contexts");

  map.clear();  // empty the map
  const IndexSet& s = scrt.getMap().getIndexSet();
  map.insert(s);
  
  for (long i = s.first(); i <= s.last(); i = s.next(i)) 
    context.ithModulus(i).FFT(map[i],scrt.getMap()[i]); // compute FFT image
  return *this;
}

void DoubleCRT::toSingleCRT(SingleCRT& scrt, const IndexSet& s) const 
{
  if (&context != &scrt.getContext())
    Error("DoubleCRT::toSingleCRT -- incompatible contexts");

  IndexSet s1 = s & map.getIndexSet();
  scrt.map.clear();
  scrt.map.insert(s1);

  for (long i = s1.first(); i <= s1.last(); i = s1.next(i)) 
    context.ithModulus(i).iFFT(scrt.map[i], map[i]); // inverse FFT
}

void DoubleCRT::toSingleCRT(SingleCRT& scrt) const 
{
  const IndexSet& s = map.getIndexSet();
  toSingleCRT(scrt, s);
}


void DoubleCRT::scaleDownToSet(const IndexSet& s)
{
  const IndexSet& indexSet = getIndexSet();

  IndexSet intersect = s & indexSet;
  IndexSet diff = indexSet / s;
  
  assert(card(intersect) > 0);
  assert(card(diff) > 0);

  ZZ diffProd = context.productOfPrimes(diff);
  *this *= (diffProd % context.ModulusP());
  
  ZZX delta;
  toPoly(delta, diff);

  //ZZX delta2 = -delta;
  //ReduceCoefficientsSlow(delta2, diffProd);
  
  long delta_len = delta.rep.length();
  ZZ factor = diffProd * InvMod(to_ZZ(diffProd % context.ModulusP()), to_ZZ(context.ModulusP()));
  for (long i = 0; i < delta_len; i++) {
    ZZ c = delta.rep[i];
    delta.rep[i] *= factor;
    delta.rep[i] -= c;
  }
  delta.normalize();
  ReduceCoefficientsSlow(delta, diffProd * context.ModulusP());
  
/*  
  cout << "pt: " << diffProd << endl;
  cout << "-c bar: " << delta2 << endl;
  ZZX delta3 = delta;
  cout << "delta: " << delta << endl;
  ReduceCoefficientsSlow(delta3, diffProd);
  cout << "delta3: " << delta3 << endl << endl;
*/
  removePrimes(diff);
  *this += delta;
  *this /= diffProd;
}

