#ifndef _FHE_SI_H_
#define _FHE_SI_H_

#include "DoubleCRT.h"
#include "Util.h"
#include <vector>
#include <fstream>
#include <NTL/ZZX.h>
#include "Plaintext.h"
#include <NTL/ZZ.h>

class Ciphertext;

class FHESISecKey {
 public:
  FHESISecKey() : context(*activeContext) {
     Init(*activeContext);
  }

  FHESISecKey(const FHEcontext &context) : context(context) {
    Init(context);
  }
    
  void Init(const FHEcontext &context);

  void Decrypt(Plaintext &plaintext, const Ciphertext &ciphertext) const;

  const FHEcontext &GetContext() const;
  size_t GetSize() const;

  const vector<DoubleCRT> &GetRepresentation() const;
  void UpdateRepresentation(vector<DoubleCRT> &rep);
  
  void Export(ofstream &out) const;
  void Import(ifstream &in);
  
  friend ostream &operator<<(ostream &os, const FHESISecKey &secKey);

 private:
  vector<DoubleCRT> sKeys;
  const FHEcontext &context;
};

class FHESIPubKey {
  friend class Ciphertext;
  const FHEcontext &context;
  vector<DoubleCRT> publicKey;
 public:
  FHESIPubKey(const FHEcontext &context) : context(context) {}
 
  FHESIPubKey(const FHESISecKey &secKey) : context(secKey.GetContext()) {
    Init(secKey);
  }

  FHESIPubKey(const FHESISecKey &secKey, const FHEcontext &context) : context(context) {
    Init(secKey);
  }
  
  void Encrypt(Ciphertext &ctxt, const Plaintext &ptxt) const;
  void Init(const FHESISecKey &secKey);

  const FHEcontext &GetContext() const;
  
  const vector<DoubleCRT> &GetRepresentation() const;
  void UpdateRepresentation(const vector<DoubleCRT> &rep);
  
  void Export(ofstream &out) const;
  void Import(ifstream &in);

  friend ostream &operator<<(ostream &os, const FHESIPubKey &secKey);
};

class KeySwitchSI {
 public:
  KeySwitchSI() : context(*activeContext) {}
  
  KeySwitchSI(FHEcontext &context) : context(context) {}
  
  KeySwitchSI(const FHESISecKey &src, const FHESISecKey &dst) : context(*activeContext) {
     Init(src, dst);
  }
  KeySwitchSI(const FHESISecKey &src, const FHESISecKey &dst, const FHEcontext &context) : context(context) {
     Init(src, dst);
  }

  KeySwitchSI(const FHESISecKey &s) : context(*activeContext) {
    InitS2(s);
  }
  KeySwitchSI(const FHESISecKey &s, const FHEcontext &context) : context(context) {
    InitS2(s);
  }
  
  KeySwitchSI(const FHESISecKey &s, unsigned k) : context(*activeContext) {
    InitAutomorph(s, k);
  }
  KeySwitchSI(const FHESISecKey &s, const FHEcontext &context, unsigned k) : context(context) {
    InitAutomorph(s, k);
  }
  
  void Init(const FHESISecKey &src, const FHESISecKey &dst);
  void InitS2(const FHESISecKey &s);
  void InitAutomorph(const FHESISecKey &s, unsigned k);
  void ApplyKeySwitch(Ciphertext &ctxt) const;
  
  const vector<vector<DoubleCRT>> &GetRepresentation() const;
  void UpdateRepresentation(const vector<vector<DoubleCRT>> &rep);
  
  void Export(ofstream &out) const;
  void Import(ifstream &in);
  
  KeySwitchSI &operator=(const KeySwitchSI &other);
  
  friend ostream &operator<<(ostream &os, const KeySwitchSI &keySwitch);
 private:
  const FHEcontext &context;
  vector<vector<DoubleCRT>> keySwitchMatrix;
};

#endif
