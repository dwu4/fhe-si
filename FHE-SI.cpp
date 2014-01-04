#include "DoubleCRT.h"
#include "FHE-SI.h"
#include <vector>
#include <NTL/vec_ZZ.h>
#include "Plaintext.h"
#include "Util.h"
#include "Ciphertext.h"
#include "Serialization.h"

void FHESIPubKey::Encrypt(Ciphertext &ctxt_, const Plaintext &ptxt) const {
  Ciphertext &ctxt = (Ciphertext &)ctxt_;
  ctxt.Initialize(2, context);
  
  ZZX small;
  for (unsigned i = 0; i < context.zMstar.phiM(); i++) {
    SetCoeff(small, i, RandomBnd(2));
  }

  DoubleCRT r(small, context);
  DoubleCRT e(context);
  
  vector<DoubleCRT> ciphertext = publicKey;
  for (size_t i = 0; i < ciphertext.size(); i++) {
    e.sampleGaussian();
    e *= context.ModulusP();
  
    ciphertext[i] *= r;
    ciphertext[i] += e;
    ciphertext[i].toPoly(ctxt[i].poly);
  }
  ctxt[0] += (context.modulusQ / context.ModulusP()) * to_ZZX(ptxt.message);

  for (size_t i = 0; i < ciphertext.size(); i++) {
    ReduceCoefficients(ctxt[i].poly, context.logQ);
  }
}

const FHEcontext &FHESIPubKey::GetContext() const {
  return context;
}

void FHESIPubKey::Init(const FHESISecKey &secKey_) {
  const FHESISecKey &secKey = (const FHESISecKey &)secKey_;
  ZZX c0, c1;
  sampleGaussian(c0, context.zMstar.phiM(), context.stdev);
  SampleRandom(c1, context.modulusQ, context.zMstar.phiM());

  ZZX tmp;
  secKey.GetRepresentation()[1].toPoly(tmp);
  tmp *= c1;

  c0 += tmp;
  rem(c0, c0, context.zMstar.PhimX());
  c1 *= -1;
  
  ReduceCoefficients(c0, context.logQ);
  ReduceCoefficients(c1, context.logQ);

  publicKey.clear();
  publicKey.push_back(DoubleCRT(c0));
  publicKey.push_back(DoubleCRT(c1));
}

const vector<DoubleCRT> &FHESIPubKey::GetRepresentation() const {
  return publicKey;
}

void FHESIPubKey::UpdateRepresentation(const vector<DoubleCRT> &rep) {
  publicKey = rep;
}

void FHESIPubKey::Export(ofstream &out) const {
  ::Export(out, publicKey);
}

void FHESIPubKey::Import(ifstream &in) {
  ::Import(in, publicKey);
}

ostream &operator<<(ostream &os, const FHESIPubKey &pubKey_) {
  const FHESIPubKey &pubKey = (const FHESIPubKey &)pubKey_;
  return (os << pubKey.publicKey[0] << ", " <<
	  pubKey.publicKey[1]);
}

void FHESISecKey::Init(const FHEcontext &context) {
  sKeys.assign(2, DoubleCRT(context));
  sKeys[0] = 1;
  
  sKeys[1].sampleHWt(64);
}

void FHESISecKey::Decrypt(Plaintext &ptxt, const Ciphertext &ctxt_) const {
  const Ciphertext &ctxt = (const Ciphertext &)ctxt_;

  vector<DoubleCRT> ctxtPolys, sKeyPolys;
  for (size_t i = 0; i < sKeys.size(); i++) {
    ctxtPolys.push_back(DoubleCRT(ctxt.GetPart(i).poly));
    sKeyPolys.push_back(DoubleCRT(sKeys[i]));
  }
  
  DoubleCRT tmp;
  DotProduct(tmp, ctxtPolys, sKeyPolys);
  
  ZZX zPtxt;
  tmp.toPoly(zPtxt);

  ZZ modulusP;
  modulusP = context.ModulusP();

  ZZ q = context.modulusQ;
  ZZ q2 = 2*q;
  for(long i = 0; i <= deg(zPtxt); i++) {
    zPtxt.rep[i] *= 2*modulusP;
    zPtxt.rep[i] += q;
    zPtxt.rep[i] /= q2;
    SetCoeff(ptxt.message, i, to_ZZ_p(zPtxt.rep[i]));
  }
}

const vector<DoubleCRT> &FHESISecKey::GetRepresentation() const {
  return sKeys;
}

const FHEcontext &FHESISecKey::GetContext() const {
  return context;
}

size_t FHESISecKey::GetSize() const {
  return sKeys.size();
}

void FHESISecKey::UpdateRepresentation(vector<DoubleCRT> &rep) {
  sKeys = rep;
}

void FHESISecKey::Export(ofstream &out) const {
  ::Export(out, sKeys);
}

void FHESISecKey::Import(ifstream &in) {
  ::Import(in, sKeys);
}

ostream &operator<<(ostream &os, const FHESISecKey &secKey_) {
  const FHESISecKey &secKey = (const FHESISecKey &)secKey_;
  for (auto it = secKey.sKeys.begin(); it != secKey.sKeys.end(); it++) {
    os << *it;
  }
  return os;
}

void KeySwitchSI::Init(const FHESISecKey &src_, const FHESISecKey &dst_) {
  const FHESISecKey &src = (const FHESISecKey &)src_;
  const FHESISecKey &dst = (const FHESISecKey &)dst_;

  #ifdef DEBUG
  assert(&src.GetContext() == &dst.GetContext());
  assert(dst.GetSize() == 2); // Target keys should be in canonical form: (1,t)
  #endif
  
  vector<DoubleCRT> s = src.GetRepresentation();
  vector<ZZX> sCoeff(s.size());
  for (size_t i = 0; i < sCoeff.size(); i++) {
    s[i].toPoly(sCoeff[i]);
  }

  DoubleCRT t = dst.GetRepresentation()[1];
  size_t n = src.GetSize();

  vector<DoubleCRT> A(context.ndigits*n);
  vector<DoubleCRT> b(context.ndigits*n);
  vector<ZZX> errors(context.ndigits*n);

  unsigned ind = 0;
  for (unsigned i = 0; i < n; i++) {
    for (unsigned j = 0; j < context.ndigits; j++, ind++) {
      ZZX poly;
      SampleRandom(poly, context.modulusQ, context.zMstar.phiM());

      A[ind] = DoubleCRT(poly, context);
      b[ind] = A[ind];
      A[ind] *= -1;

      b[ind] *= t;
      
      ZZX bCoeff;
      b[ind].toPoly(bCoeff);

      ZZX err;
      sampleGaussian(err, context.zMstar.phiM(), context.stdev);

      bCoeff += err;
      errors[ind] = err;
      bCoeff += sCoeff[i];
      
      for (long k = 0; k <= deg(sCoeff[i]); k++) {
        sCoeff[i].rep[k] <<= (8*context.decompSize);
      }

      ReduceCoefficients(bCoeff, context.logQ);
      b[ind] = DoubleCRT(bCoeff);
    }
  }
  
  keySwitchMatrix.resize(2);
  keySwitchMatrix[0] = b;
  keySwitchMatrix[1] = A;
}

void KeySwitchSI::InitS2(const FHESISecKey &s_) {
  const FHESISecKey &s = (const FHESISecKey &)s_;
  vector<DoubleCRT> sKeys = s.GetRepresentation();

  vector<DoubleCRT> tKeys;
  tKeys.assign(sKeys.size()*2-1, sKeys[1]);
  tKeys[0] = sKeys[0];
  
  for (unsigned i = 2; i < tKeys.size(); i++) {
    tKeys[i] *= tKeys[i-1];
  }
  
  FHESISecKey tensoredKey = FHESISecKey(s.GetContext());
  tensoredKey.UpdateRepresentation(tKeys);
  
  Init(tensoredKey, s);
}

void KeySwitchSI::InitAutomorph(const FHESISecKey &s_, unsigned k) {
  const FHESISecKey &s = (const FHESISecKey &)s_;
  vector<DoubleCRT> sKeys = s.GetRepresentation();
  
  FHESISecKey automorphedKey(s.GetContext());
  for (unsigned i = 0; i < sKeys.size(); i++) {
    sKeys[i].automorph(k);
  }
  automorphedKey.UpdateRepresentation(sKeys);
  Init(automorphedKey, s);
}

void KeySwitchSI::ApplyKeySwitch(Ciphertext &ctxt_) const {
  Ciphertext &ctxt = (Ciphertext &)ctxt_;
  ctxt.ScaleDown();
  ctxt.ByteDecomp();

  vector<DoubleCRT> byteDecomp(ctxt.parts.size());
  for (size_t i = 0; i < byteDecomp.size(); i++) {
    byteDecomp[i] = DoubleCRT(ctxt.parts[i].poly);
  }

  vector<CiphertextPart> newCtxt(keySwitchMatrix.size());
  for (size_t i = 0; i < keySwitchMatrix.size(); i++) {
    DoubleCRT dotProd;
    DotProduct(dotProd, keySwitchMatrix[i], byteDecomp);
    dotProd.toPoly(newCtxt[i].poly);
    ReduceCoefficients(newCtxt[i].poly, context.logQ);
  }
  
  ctxt.parts = newCtxt;
}

const vector<vector<DoubleCRT>> &KeySwitchSI::GetRepresentation() const {
  return keySwitchMatrix;
}

void KeySwitchSI::UpdateRepresentation(const vector<vector<DoubleCRT>> &rep) {
  keySwitchMatrix = rep;
}

void KeySwitchSI::Export(ofstream &out) const {
  ::Export(out, keySwitchMatrix);
}

void KeySwitchSI::Import(ifstream &in) {
  ::Import(in, keySwitchMatrix);
}

KeySwitchSI &KeySwitchSI::operator=(const KeySwitchSI &other) {
  if (&context != &other.context) {
    Error("KeySwitchSI assignment: context mismatch");
  }
  keySwitchMatrix = other.keySwitchMatrix;
  return *this;
}

ostream &operator<<(ostream &os, const KeySwitchSI &keySwitch) {
  PrintVector(keySwitch.keySwitchMatrix, os);
  return os;
}
