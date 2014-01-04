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

#ifndef _CModulus_H_
#define _CModulus_H_
/* CModulus.h - supports forward and backward length-m FFT transformations
 *
 * This is a wrapper around the bluesteinFFT routines, for one modulus q.
 * Two classes are defined here, Cmodulus for a small moduli (long) and
 * CModulus for a large ones (ZZ). These classes are otherwise identical
 * hence they are implemented using a class template.
 *
 * On initialization, it initizlies NTL's zz_pContext/ZZ_pContext for this q
 * and computes a 2m-th root of unity r mod q and also r^{-1} mod q.
 * Thereafter this class provides FFT and iFFT routines that converts between
 * time & frequency domains. Some tables are computed the first time that
 * each dierctions is called, which are then used in subsequent computations.
 * 
 * The "time domain" polynomials are represented as ZZX, whic are reduced
 * modulo Phi_m(X). The "frequency domain" are jusr vectors of integers
 * (vec_long or vec_ZZ), that store only the evaluation in primitive m-th
 * roots of unity.
 */
#include "PAlgebra.h"
#include "bluestein.h"

//NTL_CLIENT

template <class zz,class zp,class zpx,class zzv,class fftrep,class zpContext>
class Cmod {
  zz          q;       // the modulus
  zpContext   context; // NTL's tables for this modulus

  // points to the Zm* structure, m is FFT size
  const PAlgebra*   zmStar;
  // That should have been a reference-counted pointer,
  // using regular pointer is just asking for trouble

  zz          root;    // 2m-th root of unity modulo q
  zz          rInv;    // root^{-1} mod q

  zpx*        powers;  // tables for forward FFT
  fftrep*     Rb;

  zpx*        ipowers; // tables for backward FFT
  fftrep*     iRb;

  void privateInit(const PAlgebra&, const zz& rt);// Allocate memory and compute roots
  void freeSpace() {
    if (powers)  { delete powers; powers=NULL; }
    if (Rb)      { delete Rb; Rb=NULL; }
    if (ipowers) { delete ipowers; ipowers=NULL; }
    if (iRb)     { delete iRb; iRb = NULL; }
  }

 public:

  // Destructor and constructors

  ~Cmod() { freeSpace(); } // destructor

  // Default constructor
  Cmod():zmStar(NULL),powers(NULL),Rb(NULL),ipowers(NULL),iRb(NULL){}

  Cmod(const Cmod &other):
    zmStar(NULL),powers(NULL),Rb(NULL),ipowers(NULL),iRb(NULL)
  { *this = other; }

  // specify m, q, and optional root
  Cmod(const PAlgebra &zms, const zz &qq, const zz &rt):
    zmStar(NULL),powers(NULL),Rb(NULL),ipowers(NULL),iRb(NULL)
  { init(zms, qq, rt); }

  // specify only m and optional root, q is the NTL current modulus
  Cmod(const PAlgebra &zms, const zz &rt):
    zmStar(NULL),powers(NULL),Rb(NULL),ipowers(NULL),iRb(NULL)
  { init(zms, rt); }

  // a "twisted" copy constructor with a different value of m
  Cmod(const PAlgebra &zms, const Cmod &other):
    zmStar(NULL),powers(NULL),Rb(NULL),ipowers(NULL),iRb(NULL)
  { init(zms,other); }

  Cmod& operator=(const Cmod &other) {
    q       = other.q;
    context = other.context;
    zmStar  =  other.zmStar; // Yes, really copy this point

    root = other.root;
    rInv = other.rInv;

    // copy data, not pointers in these fields
    freeSpace(); // just in case
    if (other.powers)  powers  = new zpx(*(other.powers));
    if (other.Rb)      Rb      = new fftrep(*(other.Rb));
    if (other.ipowers) ipowers = new zpx(*(other.ipowers));
    if (other.iRb)     iRb     = new fftrep(*(other.iRb));

    return *this;
  }

  // Initializing an empty Cmod object

  // specify m, q, and optional root
  void init(const PAlgebra &zms, zz qq, const zz &rt) {
    if (zmStar) return; // do not overwrite an initialized object
    context = zpContext(qq);
    privateInit(zms,rt);
  }

  // specify only m and optional root, q is the NTL current modulus
  void init(const PAlgebra &zms, const zz &rt) {
    if (zmStar) return; // do not overwrite an initialized object
    context.save();
    privateInit(zms,rt); 
  }

  // "twisted" copy with a different value of m
  void init(const PAlgebra &zms, const Cmod &other) {
    if (zmStar) return; // do not overwrite an initialized object

    if (other.getM() == zms.M()) {
      *this = other; return; // just copy the whole thing
    }

    zz rt;
    context = other.context;
    if (other.getM()%zms.M()==0) { // use the info in other
      context.restore();  // setup zp::modulus()
      zp rtp;
      conv(rtp,other.root);
      power(rtp, rtp, other.getM()/zms.M());
      rt = rep(rtp);
    }
    else conv(rt,0);

    privateInit(zms, rt);
  }

  // utility methods

  const PAlgebra &ZmStar() const { return *zmStar; }
  unsigned getM() const    { return zmStar->M(); }
  unsigned getPhiM() const { return zmStar->phiM(); }
  const zz& getQ() const          { return q; }
  const zz& getRoot() const       { return root; }

  void restoreModulus() {context.restore();} // restore NTL's current modulus

  // FFT routines

  void FFT(zzv &y, const ZZX& x) const;  // y = FFT(x)
  void iFFT(ZZX &x, const zzv& y) const; // x = FFT^{-1}(y)
};

typedef Cmod<long,zz_p,zz_pX,vec_long,fftRep,zz_pContext> Cmodulus;
typedef Cmod<ZZ,  ZZ_p,ZZ_pX,vec_ZZ,  FFTRep,ZZ_pContext> CModulus;

#endif // ifdef _CModulus_H_
