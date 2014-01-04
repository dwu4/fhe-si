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
 
#include "NumbTh.h"

#include <fstream>
#include <cassert>

using namespace std;

// Factoring by trial division, only works for N<2^{60}.
// Only the primes are recorded, not their multiplicity
template<class zz> static void factorT(vector<zz> &factors, const zz &N)
{
  factors.resize(0); // reset the factors

  if (N<2) return;   // sanity check

  PrimeSeq s;
  zz n = N;
  while (true) {
    if (ProbPrime(n)) { // we are left with just a single prime
      factors.push_back(n);
      return;
    }
    // if n is a composite, check if the next prime divides it
    long p = s.next();
    if ((n%p)==0) {
      zz pp;
      conv(pp,p);
      factors.push_back(pp);
      do { n /= p; } while ((n%p)==0);
    }
    if (n==1) return;
  }
}
void factorize(vector<long> &factors, long N) { factorT<long>(factors, N);}
void factorize(vector<ZZ> &factors, const ZZ& N) {factorT<ZZ>(factors, N);}

template<class zz> static void phiNT(zz &phin, vector<zz> &facts, const zz &N)
{
  if (facts.size()==0) factorize(facts,N);

  zz n = N;
  conv(phin,1); // initialize phiN=1
  for (unsigned i=0; i<facts.size(); i++) {
    zz p = facts[i];
    phin *= (p-1); // first factor of p
    for (n /= p; (n%p)==0; n /= p) phin *= p; // multiple factors of p
  } 
}
// Specific template instantiations for long and ZZ
void phiN(long &pN, vector<long> &fs, long N)  { phiNT<long>(pN,fs,N); }
void phiN(ZZ &pN, vector<ZZ> &fs, const ZZ &N) { phiNT<ZZ>(pN,fs,N);   }

/* Compute Phi(N) */
int phi_N(int N)
{
  int phiN=1,p,e;
  PrimeSeq s;
  while (N!=1)
    { p=s.next();
      e=0;
      while ((N%p)==0) { N=N/p; e++; }
      if (e!=0)
        { phiN=phiN*(p-1)*power_long(p,e-1); }
    }
  return phiN;
}

// finding e-th root of unity modulo the current modulus
template<class zp,class zz> void FindPrimRootT(zp &root, unsigned e)
{
  zz phiQ;            // phi(q)
  { vector<zz> facts; // just so we can call the PhiN routine
    zz q = zp::modulus();
    phiN(phiQ, facts, q);
  }
  if ((phiQ%e)!=0) { // error, no e'th roots of unity exist
    conv(root,0); // set root = 0 and return
    return;
  }
  vector<long> facts;
  factorize(facts,e); // factorization of e
  zz exp = phiQ/e;

  // generate random s^exp, then check if they have order e
  for (int i=0; i<1000; i++) { // similar to while(true) { ... }
    zp s;
    random(s);  // s = random in Z_p
    power(root, s, exp); // try root = s^{phi(q)/e} mod p
    power(s, root, e);   // this should be 1, if we chose (s,Q)==1
    if (s != 1) continue;

    // check that s^{e/p} != 1 for any prime divisor p of e
    for (unsigned i=0; i<facts.size(); i++) {
      long e2 = e/facts[i];
      power(s, root, e2);   // s = root^{e/p}
      if (s == 1) break;    // not an e-th root
    }
    if (s != 1) return; // root is indeed an e-th root
  }
  // we should never get here
  Error("FindPrimitiveRoot(): gave up after 1000 trials");
}
// instantiations of the template
void FindPrimitiveRoot(zz_p &r, unsigned e) {FindPrimRootT<zz_p,long>(r,e);}
void FindPrimitiveRoot(ZZ_p &r, unsigned e) {FindPrimRootT<ZZ_p,ZZ>(r,e);}

/* Compute mobius function (naive method as n is small) */
int mobius(int n)
{
  int p,e,arity=0;
  PrimeSeq s;
  while (n!=1)
    { p=s.next();
      e=0;
      while ((n%p==0)) { n=n/p; e++; }
      if (e>1) { return 0; }
      if (e!=0) { arity^=1; }
    }     
  if (arity==0) { return 1; }
  return -1;
}



/* Compute cyclotomic polynomial */
ZZX Cyclotomic(int N)
{
  ZZX Num,Den,G,F;
  set(Num); set(Den);
  int m,d;
  for (d=1; d<=N; d++)
    { if ((N%d)==0)
         { clear(G);
           SetCoeff(G,N/d,1); SetCoeff(G,0,-1);
           m=mobius(d);
           if (m==1)       { Num*=G; }
           else if (m==-1) { Den*=G; }
         }
    } 
  F=Num/Den;
  return F;
}



/* Find a primitive root modulo N */
int primroot(int N,int phiN)
{
  int g=2,p;
  PrimeSeq s;
  bool flag=false;

  while (flag==false)
    { flag=true;
      s.reset(1);
      do
        { p=s.next();
          if ((phiN%p)==0)
            { if (PowerMod(g,phiN/p,N)==1)
                { flag=false; }
            }
        }
      while (p<phiN && flag);
      if (flag==false) { g++; }
    }
  return g;
}


int ord(int N,int p)
{
  int o=0;
  while ((N%p)==0)
    { o++;
      N/=p;
    }
  return o;
}




ZZX RandPoly(int n,const ZZ& p)
{ 
  ZZX F; F.SetMaxLength(n);
  ZZ p2;  p2=p>>1;
  for (int i=0; i<n; i++)
    { SetCoeff(F,i,RandomBnd(p)-p2); }
  return F;
}


/* When q=2 maintains the same sign as the input */
void PolyRed(ZZX& out, const ZZX& in, const ZZ& q, bool abs)
{
  // ensure that out has the same degree as in
  out.SetMaxLength(deg(in)+1);               // allocate space if needed
  if (deg(out)>deg(in)) trunc(out,out,deg(in)+1); // remove high degrees

  ZZ q2; q2=q>>1;
  for (int i=0; i<=deg(in); i++)
    { ZZ c=coeff(in,i);
      c %= q;
      if (abs) {
        if (c<0) c += q;
      } 
      else if (q!=2) {
        if (c>q2)  { c=c-q; }
          else if (c<-q2) { c=c+q; }
      }
      else // q=2
        { if (sign(coeff(in,i))!=sign(c))
	    { c=-c; }
        }
      SetCoeff(out,i,c);
    }
}


void PolyRed(ZZX& out, const ZZX& in, int q, bool abs)
{
  // ensure that out has the same degree as in
  out.SetMaxLength(deg(in)+1);               // allocate space if needed
  if (deg(out)>deg(in)) trunc(out,out,deg(in)+1); // remove high degrees

  int q2; q2=q>>1;
  for (int i=0; i<=deg(in); i++)
    { int c=coeff(in,i)%q;
      if (abs)
        { if (c<0) { c=c+q; } }
      else if (q==2)
        { if (coeff(in,i)<0) { c=-c; } }
      else
        { if (c>q2)  { c=c-q; }
          else if (c<-q2) { c=c+q; }
	}
      SetCoeff(out,i,c);
    }
}



GF2X to_GF2X(const ZZX& a)
{
  GF2X ans; 
  if (!IsZero(a)) 
    { ans.SetMaxLength(deg(a));
      for (int i=0; i<=deg(a); i++)
        { SetCoeff(ans,i,coeff(a,i)%2); }
    }
  return ans;
}



ZZX to_ZZX(const GF2X& a)
{
  ZZX ans; 
  if (!IsZero(a)) 
    { ans.SetMaxLength(deg(a));
      for (int i=0; i<=deg(a); i++)
        { SetCoeff(ans,i,rep(coeff(a,i))); }
    }
  return ans;
}



int is_in(int x,int* X,int sz)
{
  for (int i=0; i<sz; i++)
    { if (x==X[i]) { return i; } }
  return -1;
}

/* Incremental integer CRT for vectors. Expects co-primes p>0,q>0 with q odd,
 * and such that all the entries in vp are in [-p/2,p/2) and all entries in
 * vq are in [0,q-1). Returns in vp the CRT of vp mod p and vq mod q, as
 * integers in [-pq/2, pq/2). Uses the formula:
 *
 *   CRT(vp,p,vq,q) = vp + p*[ (vq-vp)*p^{-1} ]_q
 *
 * where [...]_q means reduction to the interval [-q/2,q/2). As q is odd then
 * this is the same as reducing to [-(q-1)/2,(q-1)/2], hence [...]_q * p is
 * in [-p(q-1)/2, p(q-1)/2], and since vp is in [-p/2,p/2) then the sum is
 * indeed in [-pq/2,pq/2).
 *
 * Returns true if both vectors are of the same length, false otherwise
 */
template <class zzvec>
bool intVecCRT(vec_ZZ& vp, const ZZ& p, const zzvec& vq, long q)
{
  long pInv = InvMod(rem(p,q), q); // p^{-1} mod q
  long n = min(vp.length(),vq.length());
  long q_over_2 = q/2;
  ZZ tmp;
  long vqi;
  for (long i=0; i<n; i++) {
    conv(vqi, vq[i]); // convert to single precision
    long vq_minus_vp_mod_q = SubMod(vqi, rem(vp[i],q), q);

    long delta_times_pInv = MulMod(vq_minus_vp_mod_q, pInv, q);
    if (delta_times_pInv > q_over_2) delta_times_pInv -= q;

    mul(tmp, delta_times_pInv, p); // tmp = [(vq_i-vp_i)*p^{-1}]_q * p
    vp[i] += tmp;
  }
  // other entries (if any) are 0 mod q
  for (long i=vq.length(); i<vp.length(); i++) {
    long minus_vp_mod_q = NegateMod(rem(vp[i],q), q);

    long delta_times_pInv = MulMod(minus_vp_mod_q, pInv, q);
    if (delta_times_pInv > q_over_2) delta_times_pInv -= q;

    mul(tmp, delta_times_pInv, p); // tmp = [(vq_i-vp_i)*p^{-1}]_q * p
    vp[i] += tmp;
  }
  return (vp.length()==vq.length());
}
// specific instantiations: vq can be vec_long or vec_ZZ
template bool intVecCRT(vec_ZZ&, const ZZ&, const vec_ZZ&, long);
template bool intVecCRT(vec_ZZ&, const ZZ&, const vec_long&, long);

void sampleHWt(ZZX &poly, long Hwt, long n)
{
  if (n<=0) n=deg(poly)+1; if (n<=0) return;
  clear(poly);          // initialize to zero
  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial

  long b,u,i=0;
  if (Hwt>n) Hwt=n;
  while (i<Hwt) {  // continue until exactly Hwt nonzero coefficients
    u=lrand48()%n; // The next coefficient to choose
    if (IsZero(coeff(poly,u))) { // if we didn't choose it already
      b = lrand48()&2; // b random in {0,2}
      b--;             //   random in {-1,1}
      SetCoeff(poly,u,b);

      i++; // count another nonzero coefficient
    }
  }
  poly.normalize(); // need to call this after we work on the coeffs
}

void sampleSmall(ZZX &poly, long n)
{
  if (n<=0) n=deg(poly)+1; if (n<=0) return;
  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial

  for (int i=0; i<n; i++) {    // Chosse coefficients, one by one
    long u = lrand48();
    if (u&1) {                 // with prob. 1/2 choose between -1 and +1
      u = (u & 2) -1;
      SetCoeff(poly, i, u);
    }
    else SetCoeff(poly, i, 0); // with ptob. 1/2 set to 0
  }
  poly.normalize(); // need to call this after we work on the coeffs
}

void sampleGaussian(ZZX &poly, long n, double stdev)
{
  static double const Pi=4.0*atan(1.0); // Pi=3.1415..
  static long const bignum = 0xfffffff;

  if (n<=0) n=deg(poly)+1; if (n<=0) return;
  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial
  for (int i=0; i<n; i++) SetCoeff(poly, i, ZZ::zero());

  // Uses the Box-Muller method to get two Normal(0,stdev^2) variables
  for (int i=0; i<n; i+=2) {
    double r1 = (1+RandomBnd(bignum))/((double)bignum+1);
    double r2 = (1+RandomBnd(bignum))/((double)bignum+1);
    double theta=2*Pi*r1;
    double rr= sqrt(-2.0*log(r2))*stdev;

    assert(rr < 8*stdev); // sanity-check, no more than 8 standard deviations

    // Generate two Gaussians RV's, rounded to integers
    long x = (long) floor(rr*cos(theta) +0.5);
    SetCoeff(poly, i, x);
    if (i+1 < n) {
      x = (long) floor(rr*sin(theta) +0.5);
      SetCoeff(poly, i+1, x);
    }
  }
  poly.normalize(); // need to call this after we work on the coeffs
}

// ModComp: a pretty lame implementation

void ModComp(ZZX& res, const ZZX& g, const ZZX& h, const ZZX& f)
{
  assert(LeadCoeff(f) == 1);

  ZZX hh = h % f;
  ZZX r = to_ZZX(0);

  for (long i = deg(g); i >= 0; i--) 
    r = (r*hh + coeff(g, i)) % f; 

  res = r;
}

ZZ largestCoeff(const ZZX& f)
{
  ZZ mx = ZZ::zero();
  for (long i=0; i<=deg(f); i++) {
    if (mx < abs(coeff(f,i)))
      mx = abs(coeff(f,i));
  }
  return mx;    
}
