#include <NTL/ZZ.h>

#include "FHEContext.h"
#include "Ciphertext.h"
#include "Util.h"

bool CiphertextPart::operator==(const CiphertextPart &other) const {
  return poly == other.poly;
}

CiphertextPart &CiphertextPart::operator+=(const ZZX &other) {
  poly += other;
  return *this;
}

CiphertextPart &CiphertextPart::operator+=(const CiphertextPart &other) {
  poly += other.poly;
  return *this;
}

CiphertextPart &CiphertextPart::operator*=(long l) {
  for (long i = 0; i <= deg(poly); i++) {
    poly.rep[i] *= l;
    Reduce(poly.rep[i], context.logQ);
  }
  return *this;
}

CiphertextPart &CiphertextPart::operator*=(const ZZX &other) {
  poly *= other;
  rem(poly, poly, context.zMstar.PhimX());
  for (long i = 0; i <= deg(poly); i++) {
    Reduce(poly.rep[i], context.logQ);
  }
  return *this;
}

CiphertextPart &CiphertextPart::operator*=(const ZZ_pX &other) {
  return operator*=(to_ZZX(other));
}

CiphertextPart &CiphertextPart::operator*=(const CiphertextPart &other) {
  poly *= other.poly;
  return *this;
}

CiphertextPart &CiphertextPart::operator%=(const ZZ &mod) {
  for (long i = 0; i <= deg(poly); i++) {
    poly.rep[i] %= mod;
  }
  return *this;
}

CiphertextPart &CiphertextPart::operator>>=(long k) {
  DoubleCRT tmp(poly);
  tmp >>= k;
  tmp.toPoly(poly);
  return *this;
}

CiphertextPart &CiphertextPart::operator=(const CiphertextPart &other) {
  if (&context != &other.context) {
    Error("Incompatible contexts.");
  }
  this->poly = other.poly;
  return *this;
}

ostream &operator<<(ostream &os, const CiphertextPart &ctxt) {
  return (os << ctxt.poly);
}

void Ciphertext::Initialize(unsigned n, const FHEcontext &context) {
  this->context = &context;
  parts.assign(n, CiphertextPart(context));
}

unsigned Ciphertext::size() const {
  return scaledUp ? tProd.size() : parts.size();
}

void Ciphertext::ByteDecompPart(vector<CiphertextPart>::iterator &decompPosition,
				const CiphertextPart &part) {
  long nbytes = context->ndigits * context->decompSize;
  unsigned char bytes[nbytes];
  
  unsigned char zeroes[context->decompSize];
  memset(zeroes, 0, sizeof(zeroes));
  
  long degree = deg(part.poly);
  
  for (long i = 0; i <= degree; i++) {
    ZZ coef = coeff(part.poly, i);
    Reduce(coef, context->logQ, true);
    BytesFromZZ(bytes, coef, nbytes);
    
    for(unsigned digit = 0; digit < context->ndigits; digit++) {
      unsigned char *num = &bytes[digit*context->decompSize];
      
      if(memcmp(num, zeroes, context->decompSize)) {
        SetCoeff(decompPosition[digit].poly, i, ZZFromBytes(num, context->decompSize));
      }
    }
  }
}

Ciphertext &Ciphertext::ByteDecomp() {
  vector<CiphertextPart> origParts = parts;

  parts.clear();
  parts.resize(origParts.size()*context->ndigits);

  vector<CiphertextPart>::iterator decompPosition = parts.end();
  
  for(vector<CiphertextPart>::reverse_iterator rit = origParts.rbegin();
      rit != origParts.rend(); ++rit) {
    decompPosition -= context->ndigits;
    ByteDecompPart(decompPosition, *rit);
  }
  return (Ciphertext&)*this;
}

Ciphertext &Ciphertext::operator+=(const Ciphertext &other_) {
  const Ciphertext &other = (const Ciphertext &)other_;
  
  if (!scaledUp) {
    for (unsigned i = 0; i < parts.size(); i++) {
      parts[i] += other.parts[i];
      ReduceCoefficients(parts[i].poly, context->logQ);
    }
  } else{
    for (unsigned i = 0; i < tProd.size(); i++) {
      tProd[i] += other.tProd[i];
    }
  }
  return *(Ciphertext* const)this;
}

Ciphertext &Ciphertext::operator+=(const ZZX &other) {
  ZZX scaledConstant(other);
  for (int i = 0; i <= deg(scaledConstant); i++) {
    scaledConstant.rep[i] <<= context->logQ;
    scaledConstant.rep[i] /= to_ZZ(context->ModulusP());
  }

  if (!scaledUp) {
    parts[0] += scaledConstant;
    ReduceCoefficients(parts[0].poly, context->logQ);
  } else {
    tProd[0] += scaledConstant;
  }
  return *(Ciphertext* const)this;
}

Ciphertext &Ciphertext::operator+=(const ZZ_pX &other) {
  return operator+=(to_ZZX(other));
}

Ciphertext &Ciphertext::operator*=(const Ciphertext &other_) {
  const Ciphertext &other = (const Ciphertext &)other_;
  // Convert from coefficient to DoubleCRT
  vector<DoubleCRT> c1(parts.size()), c2(other.parts.size());
  for (unsigned i = 0; i < parts.size(); i++) {
    c1[i] = DoubleCRT(parts[i].poly * context->ModulusP());
  }

  for (unsigned i = 0; i < other.parts.size(); i++) {
    c2[i] = DoubleCRT(other.parts[i].poly);
  }

  // Evaluate tensor product using DoubleCRT representation
  tProd.assign(c1.size()+c2.size()-1, DoubleCRT(*context));
  for (unsigned i = 0; i < c1.size(); i++) {
    for (unsigned j = 0; j < c2.size(); j++) {
      DoubleCRT tmp = c1[i];
      tmp *= c2[j];
      tProd[i+j] += tmp;
    }
  }
  
  parts.clear();
  scaledUp = true;

  return *(Ciphertext* const)this;
}

Ciphertext &Ciphertext::operator*=(Ciphertext &other) {
  return operator*=((const Ciphertext &)other);
}

void Ciphertext::ScaleDown() {
  if (!scaledUp) return;
  
  ZZ q = activeContext->modulusQ;
  ZZ q2 = 2*q;
  
  // Convert back to coefficient to round
  parts.clear();
  for (unsigned i = 0; i < tProd.size(); i++) {
    ZZX part;
    tProd[i].toPoly(part);

    for (long j = 0; j <= deg(part); j++) {
      part.rep[j] *= 2;
      part.rep[j] += q;
      part.rep[j] /= q2;
    }

    part.normalize();
    ReduceCoefficients(part, context->logQ);
    parts.push_back(CiphertextPart(part));
  }
  scaledUp = false;
  tProd.clear();
}

void Ciphertext::Clear() {
  tProd.clear();
  scaledUp = false;
  
  parts.clear();
}

Ciphertext &Ciphertext::operator*=(long l) {
  if (!scaledUp) {
    for (unsigned i = 0; i < parts.size(); i++) {
      parts[i] *= l;
    }
  } else {
    for (unsigned i = 0; i < tProd.size(); i++) {
      tProd[i] *= l;
    }
  }
  return *(Ciphertext* const)this;
}

Ciphertext &Ciphertext::operator*=(const ZZX &other) {
  if (!scaledUp) { // May want to consider doing this in dCRT
    for (unsigned i = 0; i < parts.size(); i++) {
      parts[i] *= other;
    }
  } else {
    DoubleCRT otherCRT(other, *context);
    for (unsigned i = 0; i < tProd.size(); i++) {
      tProd[i] *= otherCRT;
    }
  }
  return *(Ciphertext * const)this;
}

Ciphertext &Ciphertext::operator*=(const ZZ_pX &other) {
  return operator*=(to_ZZX(other));
}

Ciphertext &Ciphertext::operator>>=(long k) {
  if (!scaledUp) {
    for (unsigned i = 0; i < parts.size(); i++) {
      parts[i] >>= k;
    }
  } else {
    for (unsigned i = 0; i < tProd.size(); i++) {
      tProd[i] >>= k;
    }
  }
  return *(Ciphertext* const)this;
}

Ciphertext &Ciphertext::operator=(const Ciphertext &other_) {
  const Ciphertext &other = (const Ciphertext &)other_;
  this->context = other.context;
  this->parts = other.parts;
  this->tProd = other.tProd;
  this->scaledUp = other.scaledUp;
  return *(Ciphertext* const)this;
}

CiphertextPart Ciphertext::GetPart(unsigned ind) const {
  return parts[ind];
}

CiphertextPart &Ciphertext::operator[](unsigned ind) {
  return parts[ind];
}

ostream &operator<<(ostream &os, const Ciphertext &ctxt) {
  if (!ctxt.scaledUp) {
    for (unsigned i = 0; i < ctxt.size(); i++) {
      os << ctxt.GetPart(i) << ", ";
    }
  } else {
    PrintVector(ctxt.tProd);
  }
  return os;
}
