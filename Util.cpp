#include "Util.h"

void Reduce(ZZ &val, unsigned logQ, bool positive) {
  /* This can all be pre-computed, given constant q */
  ZZ signMask, dataMask;
  signMask = 1;
  dataMask = 1;
 
  signMask <<= logQ - 1;
  dataMask <<= logQ;
  dataMask -= 1;
  /***********************************************/
  //val += q;
  
  if(sign(val) == -1) {
    ZZ factor = val>>logQ;
    //cout << "val/q: " << factor << endl;
    val += (1-(val>>logQ)) << logQ;
    //cout << "New val: " << val << endl;
  }
  val &= dataMask;
  if(!positive) {
    val ^= signMask;
    val -= signMask;
  }
}

void ReduceCoefficients(ZZX &poly, unsigned logQ, bool positive) {
  for (long i = 0; i <= deg(poly); i++) {
    Reduce(poly.rep[i], logQ, positive);
  }
  poly.normalize();
}

void ReduceCoefficientsSlow(ZZX &poly, const ZZ &modulus, bool positive) {
  for (long i = 0; i <= deg(poly); i++) {
    poly.rep[i] %= modulus;
    if (!positive && poly.rep[i] > modulus/2) {
      poly.rep[i] -= modulus;
    }
  }
  poly.normalize();
}

void ReduceCoefficientsSlow(ZZX &poly, unsigned modulus, bool positive) {
  ReduceCoefficientsSlow(poly, to_ZZ(modulus), positive);
}

void SampleRandom(ZZX &poly, const ZZ &modulus, unsigned deg) {
  ZZ offset = modulus / 2;
  poly.SetMaxLength(deg);
  for (unsigned i = 0; i < deg; i++) {
    SetCoeff(poly, i, RandomBnd(modulus) - offset);
  }
}

void GetConstantTerm(ZZ &res, ZZ_pX &poly) {
  if (deg(poly) == -1) {
    res = 0;
  } else {
    res = rep(poly.rep[0]);
  }
}

