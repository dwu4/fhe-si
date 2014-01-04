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

/* 
 * CModulus.cpp - supports forward and backward length-m FFT transformations
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
#include "NumbTh.h"
#include "CModulus.h"

// Some simple functions that should have been provided by NTL but are not
inline bool IsZero(long i) { return (i==0); }
inline void conv(NTL::vec_zz_p& to, NTL::vec_long& from)
{
  to.SetLength(from.length());
  for (long i=0; i<from.length(); i++) conv(to[i], from[i]);
}
inline void conv(NTL::vec_long& to, NTL::vec_zz_p& from)
{
  to.SetLength(from.length());
  for (long i=0; i<from.length(); i++) to[i]=rep(from[i]);
}

inline void SetCoeff(NTL::ZZ_pX& poly, unsigned int idx, const NTL::ZZ& val)
{ SetCoeff(poly, idx, to_ZZ_p(val)); }

// It is assumed that m,q,context, and root are already set. If root is set
// to zero, it will be computed by the compRoots() method. Then rInv is
// computed as the inverse of root.

template <class zz, class zp, class zpx, class zzv, class fftrep, class zpContext>
void Cmod<zz,zp,zpx,zzv,fftrep,zpContext>::privateInit(const PAlgebra& zms,
						       const zz& rt)
{
  context.restore(); // set NTL's current modulus
  zmStar = &zms;
  q = zp::modulus();
  root = rt;

  // First find a 2m-th root of unity modulo q, if not given
  if (IsZero(root)) {
    context.restore(); // Set the current modulus to q
    zp rtp;
    unsigned e = 2*getM();
    FindPrimitiveRoot(rtp,e);
    if (IsZero(rtp)) // sanity check
      Error("Cmod::compRoots(): no 2m'th roots of unity mod q");
    root = rep(rtp);
  }
  rInv = InvMod(root,q); // set rInv = root^{-1} mod q

  // allocate memory (current modulus was defined above)
  freeSpace();  // just in case
  powers  = new zpx();
  Rb      = new fftrep();
  ipowers = new zpx();
  iRb     = new fftrep();
}


template <class zz,class zp,class zpx,class zzv,class fftrep,class zpContext>
void  Cmod<zz,zp,zpx,zzv,fftrep,zpContext>::FFT(zzv &y, const ZZX& x) const
{
  context.restore();
  zp rt;
  zpx in, out;

  conv(in,x);      // convert input to zpx format
  conv(rt, root);  // convert root to zp format

  BluesteinFFT(out, in, getM(), rt, *powers, *Rb); // call the FFT routine

  // copy the result to the output vector y, keeping only the
  // entries corresponding to primitive roots of unity
  y.SetLength(zmStar->phiM());
  unsigned i,j;
  for (i=j=0; i<getM(); i++)
    if (zmStar->inZmStar(i)) y[j++] = rep(coeff(out,i));
}

template <class zz,class zp,class zpx,class zzv,class fftrep,class zpContext>
void  Cmod<zz,zp,zpx,zzv,fftrep,zpContext>::iFFT(ZZX &x, const zzv& y) const
{
  context.restore();
  zp rt;
  zpx in, out, pwrs;

  // convert input to zpx format, initializing only the coeffs i s.t. (i,m)=1
  in.SetMaxLength(getM());
  unsigned i,j;
  for (i=j=0; i<getM(); i++)
    if (zmStar->inZmStar(i)) SetCoeff(in, i, y[j++]);
  in.normalize();
  conv(rt, rInv);  // convert rInv to zp format

  BluesteinFFT(out, in, getM(), rt, *ipowers, *iRb); // call the FFT routine
  out /= getM();                                // normalization

  // reduce the result mod (Phi_m(X),q) and copy to the output polynomial x
  conv(in, zmStar->PhimX());  // convert Phi_m(X) to zpx format
  rem(out, out, in); // out %= (Phi_m(X),q)

  conv(x,out); // convert output to ZZX format
}


// instantiating the template classes
template class Cmod<long,zz_p,zz_pX,vec_long,fftRep,zz_pContext>; // small q
template class Cmod<ZZ,  ZZ_p,ZZ_pX,vec_ZZ,  FFTRep,ZZ_pContext>; // large q


