#ifndef _PLAINTEXT_H
#define _PLAINTEXT_H

#include "NTL/ZZX.h"
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <vector>
#include "FHEContext.h"

class Plaintext {
  public:
    Plaintext() : context(*activeContext) {}
    Plaintext(const FHEcontext &context) : context(context) {}
    
    Plaintext(const ZZ_pX &msg) : context(*activeContext) { Init(msg); }
    Plaintext(const FHEcontext &context, const ZZ_pX &msg) : context(context) { Init(msg); }
    
    template <typename T>
    Plaintext(const T &msg) : context(*activeContext) { Init(to_ZZ_pX(msg)); }
    
    template <typename T>
    Plaintext(const FHEcontext &context, const T &msg) : context(context) { Init(to_ZZ_pX(msg)); }
    
    Plaintext(const vector<ZZ_pX> &msgs) : context(*activeContext) { Init(msgs); }
    Plaintext(const FHEcontext &context, const vector<ZZ_pX> &msgs) : context(context) { Init(msgs); }
    
    template <typename T>
    Plaintext(const vector<T> &msgs) : context(*activeContext) { Init(msgs); }
    
    template <typename T>
    Plaintext(const FHEcontext &context, const vector<T> &msgs) : context(context) { Init(msgs); }
    
    void Init();
    void Init(const ZZ_pX &msg);
    
    void Init(const vector<ZZ_pX> &msgs);
    
    template <typename T>
    void Init(const vector<T> &msgs) {
      vector<ZZ_pX> newMsgs(msgs.size());
      for (unsigned i = 0; i < msgs.size(); i++) {
        newMsgs[i] = to_ZZ_pX(msgs[i]);
      }
      
      EmbedInSlots(newMsgs);
    }
    
    void EmbedInSlots(const vector<ZZ_pX> &msgs, bool onlyUsable = true);
    void DecodeSlots(vector<ZZ_pX> &msgBatch, bool onlyUsable = true);
    void DecodeSlot(ZZ_pX &val, unsigned slot);
    
    Plaintext &operator=(const Plaintext &other);
    
    bool operator==(const Plaintext &other) const;
    friend ostream &operator<<(ostream &os, const Plaintext &ptxt);
    
    ZZ_pX message;
    
    /*DEBUG functions*/
    
    void Randomize() {
      random(message, deg(context.zMstar.PhimX())-1);
    }
    static Plaintext Random(const FHEcontext &context) {
      Plaintext ret(context);
      ret.Randomize();
      return ret;
    }
    
    Plaintext &operator+=(const Plaintext &other) {
      assert(&context == &other.context);
      message += other.message;
      return *this;
    }

    Plaintext &operator-=(const Plaintext &other) {
      assert(&context == &other.context);
      message -= other.message;
      return *this;
    }

    Plaintext &operator*=(const Plaintext &other) {
      assert(&context == &other.context);
      MulMod(message, message, other.message, to_ZZ_pX(context.zMstar.PhimX()));
      return *this;
    }

    Plaintext &operator>>=(long k) {
      vector<ZZ_pX> plaintextArray;
      DecodeSlots(plaintextArray, false);
      vector<ZZ_pX> rotatedArray = plaintextArray;
      for (unsigned i = 0; i < plaintextArray.size(); i++) {
        rotatedArray[(i + plaintextArray.size() - k) % plaintextArray.size()] = plaintextArray[i];
      }
      EmbedInSlots(rotatedArray, false);
      return *this;
    }

    Plaintext operator+(const Plaintext &other) {
      assert(&context == &other.context);
      return Plaintext(context, message+other.message);
    }

    Plaintext operator*(const Plaintext &other) {
      assert(&context == &other.context);
      return Plaintext(context, MulMod(message, other.message, to_ZZ_pX(context.zMstar.PhimX())));
    }
    
  private:
    const FHEcontext &context;
};

#endif
