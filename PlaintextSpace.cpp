#ifndef _PLAINTEXT_SPACE_H_
#define _PLAINTEXT_SPACE_H_

#include "PlaintextSpace.h"
#include "Util.h"
#include <algorithm>
#include "assert.h"

static inline unsigned ZZToLong(const ZZ &num) {
  unsigned data;
  BytesFromZZ((unsigned char *)&data, num, sizeof(unsigned));
  return data;
}

PlaintextSpace::~PlaintextSpace() {
  for (unsigned i = 0; i < factors.size(); i++) {
    delete factors[i];
    delete crtCoeffs[i];
  }
}

void PlaintextSpace::Init(const ZZX &PhiX, const ZZ &p) {
  ZZ_p::init(p);

  this->PhiX = to_ZZ_pX(PhiX);
  this->p = p;
  
  vec_ZZ_pX crtFactors;
  SFCanZass(crtFactors, this->PhiX);
  totalSlots = crtFactors.length();

  factors.resize(totalSlots);
  for (unsigned i = 0; i < totalSlots; i++) {
    factors[i] = new ZZ_pX(crtFactors[i]);
  }
  
  usableSlots = 1;
  unsigned tmp = totalSlots;
  while (tmp > 1) {
    usableSlots <<= 1;
    tmp >>= 1;
  }
  
  crtCoeffs.resize(totalSlots);
  for (unsigned i = 0; i < totalSlots; i++) {
    crtCoeffs[i] = new ZZ_pX();
    ZZ_pX te = this->PhiX / crtFactors[i];
    te %= crtFactors[i];
    InvMod(*crtCoeffs[i], te, crtFactors[i]);
    *crtCoeffs[i] *= (this->PhiX / crtFactors[i]);
  }
  
  FindSlots(totalSlots);
}

void PlaintextSpace::Init(const ZZX &PhiX, const ZZ &p, unsigned generator) {
  this->generator = generator;
  Init(PhiX, p);
}

unsigned PlaintextSpace::GetUsableSlots() const {
  return usableSlots;
}

unsigned PlaintextSpace::GetTotalSlots() const {
  return totalSlots;
}

void PlaintextSpace::FindSlots(unsigned nSlots) {
  vector<ZZ_pX> crt(nSlots);
  for (unsigned i = 0; i < nSlots; i++) {
    crt[i] = to_ZZ_pX(i+1);
  }
  ZZ_pX permX;
  
  EmbedInSlots(permX, crt, false);
  FrobeniusMap(permX, generator);
  DecodeSlots(crt, permX, false);
  
  vector<unsigned> perm(nSlots);
  for (unsigned i = 0; i < nSlots; i++) {
    perm[i] = ZZToLong(rep(crt[i].rep[0])) - 1;
  }
  ReorderSlots(perm);
}

void PlaintextSpace::ReorderSlots(vector<unsigned> &perm) {
  // EDF doesn't always return the factors in a particular order, so
  // we place the 0 slot first to make it deterministic
  unsigned zeroInd = 0;
  for (; zeroInd < perm.size() && perm[zeroInd] != 0; zeroInd++);
  
  vector<ZZ_pX*> factorsP;
  vector<ZZ_pX*> crtCoeffsP;
  
  factorsP.push_back(factors[zeroInd]);
  crtCoeffsP.push_back(crtCoeffs[zeroInd]);
  for (unsigned i = perm[zeroInd]; i != zeroInd; i = perm[i]) {
    factorsP.push_back(factors[i]);
    crtCoeffsP.push_back(crtCoeffs[i]);
  }
  
  assert(factorsP.size() == perm.size());
  
  swap(crtCoeffs, crtCoeffsP);
  swap(factors, factorsP);
  
  crtCoeffsP.clear();
  factorsP.clear();
}

void PlaintextSpace::EmbedInSlots(ZZ_pX &embedded, const vector<ZZ_pX> &msgs,
                                  bool onlyUsable) const {
  embedded = ZZ_pX::zero();
  unsigned msgInd = 0;
  for (unsigned i = 0; i < totalSlots && msgInd < msgs.size(); i++) {
    if (onlyUsable && i >= usableSlots) break;
    embedded += *crtCoeffs[i] * msgs[msgInd++];
  }
  embedded %= PhiX;
}

void PlaintextSpace::DecodeSlots(vector<ZZ_pX> &msgBatch, const ZZ_pX &msg,
                                 bool onlyUsable) const {
  msgBatch.resize(totalSlots);
  for (unsigned i = 0; i < totalSlots; i++) {
    if (onlyUsable && i >= usableSlots) break;
    DecodeSlot(msgBatch[i], msg, i);
  }
}

void PlaintextSpace::DecodeSlot(ZZ_pX &val, const ZZ_pX &msg, unsigned ind) const {
  rem(val, msg, *factors[ind]);
}

void PlaintextSpace::FrobeniusMap(ZZ_pX &poly, long k) {
  for (long i = deg(poly); i >= 0; i--) {
    SetCoeff(poly, k*i, poly.rep[i]);
    if (i > 0) {
      poly.rep[i] = 0;
    }
  }
  rem(poly, poly, PhiX);
}

#endif
