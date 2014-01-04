#ifndef _PLAINTEXT_SPACE_H
#define _PLAINTEXT_SPACE_H

#include "PAlgebra.h"
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <vector>

class PlaintextSpace {
  public:
    PlaintextSpace() {};
    ~PlaintextSpace();
    
    void Init(const ZZX &PhiX, const ZZ &p);
    void Init(const ZZX &PhiX, const ZZ &p, unsigned generator);
    
    unsigned GetUsableSlots() const;
    unsigned GetTotalSlots() const;
    
    void EmbedInSlots(ZZ_pX &embedded, const vector<ZZ_pX> &msgs, bool onlyUsable = true) const;
    void DecodeSlots(vector<ZZ_pX> &msgBatch, const ZZ_pX &msg, bool onlyUsable = true) const;
    void DecodeSlot(ZZ_pX &val, const ZZ_pX &msg, unsigned ind) const;
  private:
    ZZ p;
    unsigned generator;
  
    unsigned totalSlots;
    unsigned usableSlots;
    
    ZZ_pX PhiX;
    
    vector<ZZ_pX *> factors;
    vector<ZZ_pX *> crtCoeffs;
    
    void FindSlots(unsigned nSlots);
    void ReorderSlots(vector<unsigned> &perm);
    
    void FrobeniusMap(ZZ_pX &poly, long k);
  
  friend class FHEcontext;
};

#endif