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
 
#ifndef _DoubleCRT_H_
#define _DoubleCRT_H_
/* DoubleCRT.h - This class holds an integer polynomial in double-CRT form
 *
 * Double-CRT form is a matrix of L rows and phi(m) columns. The i'th row
 * contains the FFT of the element wrt the ith prime, i.e. the evaluations
 * of the polynomial at the primitive mth roots of unity mod the ith prime.
 * The polynomial thus represented is defined modulo the product of all the
 * primes in use.
 *
 * The list of primes is defined by the data member indexMap.
 * indexMap.getIndexSet() defines the set of indices of primes
 * associated with this DoubleCRT object: they index the
 * primes stored in the associated FHEContext.
 *
 * Arithmetic operations are computed modulo the product of the primes in use
 * and also modulo Phi_m(X). Arithmetic operations can only be applied to
 * DoubleCRT objects relative to the same context, trying to add/multiply
 * objects that have different FHEContext objects will raise an error.
 */


#include <vector>
#include <NTL/ZZX.h>
#include <NTL/vec_vec_long.h>
#include <NTL/matrix.h>
#include "NumbTh.h"
#include "IndexMap.h"

NTL_CLIENT

// NTL decleration of a matrix of long integers

//NTL_matrix_decl(long,vec_long,vec_vec_long,mat_long)
//NTL_eq_matrix_decl(long,vec_long,vec_vec_long,mat_long)
//NTL_io_matrix_decl(long,vec_long,vec_vec_long,mat_long)

#include "FHEContext.h"

class SingleCRT;


class DoubleCRTHelper : public IndexMapInit<vec_long> {
private: 
  long val;

public:

  DoubleCRTHelper(const FHEcontext& context) { 
    val = context.zMstar.phiM(); 
  }

  virtual void init(vec_long& v) { 
    v.FixLength(val); 
  }

  virtual IndexMapInit<vec_long> * clone() const { 
    return new DoubleCRTHelper(*this); 
  }

private:
  DoubleCRTHelper(); // disable default constructor
};



class DoubleCRT {
  const FHEcontext& context; // the context
  IndexMap<vec_long> map; // the data itself: if the i'th prime is in use then
                          // map[i] is the vector of evaluation wrt this prime

  // a "sanity check" function, verifies consistency of the map with current
  // moduli chain, an error is raised if they are not consistent
  void verify();


  // Generic operators, Fnc is either AddMod, SubMod, or MulMod. (This should
  // have been a template, but gcc refuses to cooperate.) 
  // The behavior when *this and other use different primes depends on the flag
  // matchIndexSets. When it is set to true then the effective modulus is
  // determined by the union of the two index sets; otherwise, the index set
  // of *this.

  DoubleCRT& Op(const DoubleCRT &other, long (*Fnc)(long, long, long),
		bool matchIndexSets=true);

  DoubleCRT& Op(const ZZ &num, long (*Fnc)(long, long, long));

  DoubleCRT& Op(const ZZX &poly, long (*Fnc)(long, long, long));

public:

  // Constructors and assignment operators

  // representing an integer polynomial as DoubleCRT. If the set of primes
  // to use is not specified, the resulting object uses all the primes in
  // the context. If the coefficients of poly are larger than the product of
  // the used primes, they are effectively reduced modulo that product

  DoubleCRT(const DoubleCRT& other); // copy constructor

  DoubleCRT(const ZZX&poly, const FHEcontext& _context, const IndexSet& indexSet);
  DoubleCRT(const ZZX&poly, const FHEcontext& _context);

  explicit DoubleCRT(const ZZX&poly); 
  // uses the "active context", run-time error if it is NULL
  // declare "explicit" to avoid implicit type conversion

 // Without specifying a ZZX, we get the zero polynomial

  DoubleCRT(const FHEcontext &_context, const IndexSet& indexSet);

  explicit DoubleCRT(const FHEcontext &_context);
  // declare "explicit" to avoid implicit type conversion

  DoubleCRT(); 
  // uses the "active context", run-time error if it is NULL


  // Assignment operator, the following two lines are equivalent:
  //    DoubleCRT dCRT(poly, context, indexSet);
  // or
  //    DoubleCRT dCRT(context, indexSet); dCRT = poly;

  DoubleCRT& operator=(const DoubleCRT& other);

  DoubleCRT& operator=(const SingleCRT& other);
  DoubleCRT& operator=(const ZZX& poly);
  DoubleCRT& operator=(const ZZ& num);
  DoubleCRT& operator=(const long num) { *this = to_ZZ(num); return *this; }

  // Recovering the polynomial in coefficient representation. This yields an
  // integer polynomial with coefficients in [-P/2,P/2], unless the positive
  // flag is set to true, in which case we get coefficients in [0,P-1] (P is
  // the product of all moduli used). Using the optional IndexSet param
  // we compute the polynomial reduced modulo the product of only the ptimes
  // in that set.

  void toPoly(ZZX& p, const IndexSet& s, bool positive=false) const;
  void toPoly(ZZX& p, bool positive=false) const;

  // The variant toPolyMod has another argument, which is a modulus Q, and it
  // computes toPoly() mod Q. This is offerred as a separate function in the
  // hope that one day we will figure out a more efficient method of computing
  // this. Right now it is not implemented
  // 
  // void toPolyMod(ZZX& p, const ZZ &Q, const IndexSet& s) const;


  // FIXME: comparison between objects with different contexts should be an error
  bool operator==(const DoubleCRT& other) const {
    return &context == &other.context && map == other.map;
  }

  bool operator!=(const DoubleCRT& other) const { 
    return !(*this==other);
  }

  // Set to zero, one
  DoubleCRT& SetZero() { 
    *this = ZZ::zero(); 
    return *this; 
  }

  DoubleCRT& SetOne()  { 
    *this = 1; 
    return *this; 
  }

  // expand the index set by s1.
  // it is assumed that s1 is disjoint from the current index set.
  void addPrimes(const IndexSet& s1);

  // expand index set by s1, and scale by f*(f^{-1} mod p),
  // where f is the product of all the primes in s1.
  // it is assumed that s1 is disjoint from the current index set.
  // Returns the logarithm of the product of the added factors.
  double addPrimesAndScale(const IndexSet& s1);

  // remove s1 from the index set
  void removePrimes(const IndexSet& s1) {
    map.remove(s1);
  }


  // Arithmetic operations. Only the "destructive" versions are used,
  // i.e., a += b is implemented but not a + b.

  // Addition, negation, subtraction
  DoubleCRT& operator+=(const DoubleCRT &other) {
    return Op(other, NTL::AddMod);
  }

  DoubleCRT& operator+=(const ZZX &poly) {
    return Op(poly, NTL::AddMod);
  }

  DoubleCRT& operator+=(const ZZ &num) { 
    return Op(num, NTL::AddMod);
  }

  DoubleCRT& operator+=(long num) { 
    return Op(to_ZZ(num), NTL::AddMod);
  }

  DoubleCRT& operator-=(const DoubleCRT &other) {
    return Op(other,NTL::SubMod);
  }

  DoubleCRT& operator-=(const ZZX &poly) {
    return Op(poly,NTL::SubMod);
  }
  
  DoubleCRT& operator-=(const ZZ &num) { 
    return Op(num, NTL::SubMod);
  }

  DoubleCRT& operator-=(long num) { 
    return Op(to_ZZ(num), NTL::SubMod);
  }

  // These are the prefix versions, ++dcrt and --dcrt. 
  DoubleCRT& operator++() { return (*this += 1); };
  DoubleCRT& operator--() { return (*this -= 1); };

  // These are the postfix versions -- return type is void,
  // so it is offered just for style...
  void operator++(int) { *this += 1; };
  void operator--(int) { *this -= 1; };


  // Multiplication
  DoubleCRT& operator*=(const DoubleCRT &other) {
    return Op(other,NTL::MulMod);
  }

  DoubleCRT& operator*=(const ZZX &poly) {
    return Op(poly,NTL::MulMod);
  }

  DoubleCRT& operator*=(const ZZ &num) { 
    return Op(num,NTL::MulMod);
  }

  DoubleCRT& operator*=(long num) { 
    return Op(to_ZZ(num),NTL::MulMod);
  }


  // Procedural equivalents, supporting also the matchIndexSets flag
  void Add(const DoubleCRT &other, bool matchIndexSets=true) {
    Op(other, NTL::AddMod, matchIndexSets); 
  }

  void Sub(const DoubleCRT &other, bool matchIndexSets=true) {
    Op(other, NTL::SubMod, matchIndexSets); 
  }

  void Mul(const DoubleCRT &other, bool matchIndexSets=true) {
    Op(other, NTL::MulMod, matchIndexSets); 
  }

  // Division by constant
  DoubleCRT& operator/=(const ZZ &num);
  DoubleCRT& operator/=(long num) { return (*this /= to_ZZ(num)); }


  // Small-exponent polynomial exponentiation
  void Exp(long k);


  // Apply the automorphism F(X) --> F(X^k)  (with gcd(k,m)=1)
  void automorph(long k);
  DoubleCRT& operator>>=(long k) { automorph(k); return *this; }


  // Utilities

  const FHEcontext& getContext() const { return context; }
  const IndexMap<vec_long>& getMap() const { return map; }
  const IndexSet& getIndexSet() const { return map.getIndexSet(); }

  void setMap(const IndexMap<vec_long> &newMap) { map = newMap; }
  
  // Choose random DoubleCRT's, either at random or with small/Gaussian
  // coefficients. 

  // fills each row i w/ random ints mod pi, uses NTL's PRG
  void randomize(const ZZ* seed=NULL);

  // Coefficients are -1/0/1, Prob[0]=1/2
  void sampleSmall() {
    ZZX poly; 
    ::sampleSmall(poly,context.zMstar.phiM()); // degree-(phi(m)-1) polynomial
    *this = poly; // convert to DoubleCRT
  }

  // Coefficients are -1/0/1 with pre-specified number of nonzeros
  void sampleHWt(long Hwt) {
    ZZX poly; 
    ::sampleHWt(poly,Hwt,context.zMstar.phiM());
    *this = poly; // convert to DoubleCRT
  }

  // Coefficients are Gaussians
  void sampleGaussian(double stdev=0.0) {
    if (stdev==0.0) stdev=to_double(context.stdev); 
    ZZX poly; 
    ::sampleGaussian(poly, context.zMstar.phiM(), stdev);
    *this = poly; // convert to DoubleCRT
  }


  // makes a corresponding SingleCRT object, restricted to
  // the given index set, if specified
  void toSingleCRT(SingleCRT& scrt, const IndexSet& s) const;
  void toSingleCRT(SingleCRT& scrt) const;

  // used to implement modulus switching
  void scaleDownToSet(const IndexSet& s);

  ZZ getCoefficientModulus() const {
    return context.productOfPrimes();
  }


  friend ostream &operator<<(ostream &os, const DoubleCRT &p) {
    ZZX poly;
    p.toPoly(poly);
    
    return (os << poly);
  }

#if 0
  // This makes a SingleCRT objects with the same levels as *this, converts
  // the numLvls starting levels fromLvl to ZZX format, and SETS TO ZERO all
  // the other levels in the SingleCRT object.

  // FIXME: not sure about the semantics here... VJS...


  // I/O: ONLY the matrix is outputted/recovered, not the moduli chain!! 
  // An error is raised on input if this is not consistent with the current chain

  friend istream& operator>> (istream &s, DoubleCRT &d)
  { s >> d.data; d.fixIndexSet(); d.verify(); return s; }
#endif

};




inline void conv(DoubleCRT &d, const ZZX &p) { d=p; }

inline DoubleCRT to_DoubleCRT(const ZZX& p) {
  return DoubleCRT(p);
}

inline void conv(ZZX &p, const DoubleCRT &d) { d.toPoly(p); }

inline ZZX to_ZZX(const DoubleCRT &d)  { ZZX p; d.toPoly(p); return p; }

inline void conv(DoubleCRT& d, const SingleCRT& s) { d=s; }

#endif // #ifndef _DoubleCRT_H_
