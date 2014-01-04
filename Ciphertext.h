#ifndef _CIPHERTEXT_H_
#define _CIPHERTEXT_H_

#include <vector>
#include <NTL/xdouble.h>
#include "DoubleCRT.h"
#include "Util.h"
#include "FHE-SI.h"

class FHESIPubKey;

class CiphertextPart {
  const FHEcontext &context;

 public:
  ZZX poly;

  CiphertextPart() : context(*activeContext) {}
  CiphertextPart(const FHEcontext &context) : context(context) {}
  CiphertextPart(const long val) : context(*activeContext) {
    poly = to_ZZX(val);
  }
  explicit CiphertextPart(const ZZX &poly) : context(*activeContext) {
    this->poly = poly;
  }

  bool operator==(const CiphertextPart &other) const;

  CiphertextPart &operator+=(const ZZX &other);
  CiphertextPart &operator+=(const CiphertextPart &other);

  CiphertextPart &operator*=(const ZZX &other);
  CiphertextPart &operator*=(const ZZ_pX &other);
  CiphertextPart &operator*=(const CiphertextPart &other);
  CiphertextPart &operator*=(long l);

  CiphertextPart &operator%=(const ZZ &mod);
  
  CiphertextPart &operator>>=(long k);

  CiphertextPart &operator=(const CiphertextPart &other);
  
  friend ostream &operator<<(ostream &os, const CiphertextPart &ctxt);
};

class Ciphertext {
  const FHEcontext *context;
  void ByteDecompPart(vector<CiphertextPart>::iterator &smallPolys,
                      const CiphertextPart &part);
                      
  vector<DoubleCRT> tProd;
  bool scaledUp;
 public:
  // WARNING: You should pass in a public key to construct ciphertexts! These constructors
  //          are provided for serialization purposes!
  Ciphertext() : context(activeContext) {
    scaledUp = false;
  }
 
  Ciphertext(const FHEcontext &context) : context(&context) {
    scaledUp = false;
  }
 
  Ciphertext(const FHESIPubKey &pk) : context(&pk.GetContext()) {
    scaledUp = false;
  }
  
  vector<CiphertextPart> parts;

  void Initialize(unsigned n, const FHEcontext &context);
  unsigned size() const;

  Ciphertext &ByteDecomp();
  
  Ciphertext &operator+=(const Ciphertext &other);
  Ciphertext &operator+=(const ZZX &other);
  Ciphertext &operator+=(const ZZ_pX &other);
  
  Ciphertext &operator*=(const Ciphertext &other);
  Ciphertext &operator*=(long l);

  Ciphertext &operator*=(const ZZX &other);
  Ciphertext &operator*=(const ZZ_pX &other);
  
  Ciphertext &operator>>=(long k);
  
  void Clear();
  void ScaleDown();
  void SetTensorRepresentation(vector<DoubleCRT> &repr);

  CiphertextPart GetPart(unsigned ind) const;
  CiphertextPart &operator[](unsigned ind);
  
  Ciphertext &operator=(const Ciphertext &other);
  
  friend ostream &operator<<(ostream &os, const Ciphertext &ctxt);
};

#endif
