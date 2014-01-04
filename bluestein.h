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
 
#ifndef _Bluestein
#define _Bluestein
/* bluestein.h - declaration of BluesteinFFT(x, a, n, root, powers, Rb):
 * ---------------------------------------------------------------------
 * Compute length-n FFT of the coefficient-vector of a, put the result in x
 * If the degree of a is less than n then it treats the top coefficients
 * as 0, if the degree of a is more than n then the extra coefficients are
 * ignored. Similarly, if the top entries in x are zeros then x will have
 * degree smaller than n. The argument root is a 2n-th root of unity, namely
 * BluesteinFFT(...,root,...)=DFT(...,root^2,...).
 *
 * The inverse-FFT is obtained just by calling BluesteinFFT(... root^{-1}),
 * but this procedure is *NOT SCALED*, so BluesteinFFT(x,a,n,root,...) and
 * then BluesteinFFT(a2,x,n,root^{-1},...) will result in a2=a*n.
 *
 * In addition to the size-n FFT of a which is returned in x, this procedure
 * also returns the powers of root in the powers argument:
 *    powers = [1, root, root^4, root^9, ..., root^{(n-1)^2}]
 * and in Rb it returns the size-N FFT representation of the negative
 * powers (with N>=2n-1, N a power of two):
 *    b = [0,...,0, root^{-(n-1)^2},...,root^{-4}, root^{-1}, 1, 
 *                  root^{-1},root^{-4},...,root^{-(n-1)^2}, 0...,0]
 * On subsequent calls with these 'powers' and 'Rb', these arrays are
 * not computed again but taken from these pre-comuted variables.
 *
 * If the powers and Rb arguments are initialized, then it is assumed that
 * they were computed correctly from root. The bahavior is undefined when
 * calling with initialized powers and Rb but a different root. In particular,
 * to compute the inverse-FFT (using root^{-1}), one must provide different
 * powers and Rb arguments than those that were given when computing in the
 * forward direction using root. To reset these arguments between calls
 * with different root values, use clear(powers); Rb.SetSize(0);
 *
 * This procedure assume that the global zp::modulus() was defined before
 * calling it. Also, this procedure *cannot* be used for in-place FFT,
 * calling BluesteinFFT(x,x,...) will just zero-out the polynomial x.
 *
 * This module builds on Shoup's NTL, and contains both a bigint version
 * with types ZZ_p and ZZ_pX and a smallint version with types zz_p and zz_pX.
 */
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>

NTL_CLIENT

void BluesteinFFT(ZZ_pX& x, const ZZ_pX& a, long n,
                  const ZZ_p& root, ZZ_pX& powers, FFTRep& Rb);

void BluesteinFFT(zz_pX& x, const zz_pX& a, long n,
                  const zz_p& root, zz_pX& powers, fftRep& Rb);

#endif

