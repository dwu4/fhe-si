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
 
#ifndef _FHEcontext_H_
#define _FHEcontext_H_

#include <fstream>
#include <NTL/xdouble.h>
#include "PAlgebra.h"
#include "CModulus.h"
#include "IndexSet.h"
#include "Util.h"
#include "PlaintextSpace.h"

class FHEcontext;

// Convenience routines

// Adds to the chain primes whose product is at least e^totalSize, 
// returns natural log of the product of all added primes
double AddPrimesBySize(FHEcontext& context, double totalSize,
		       bool special=false);

// Adds nPrimes primes to the chain
// returns natural log of the product of all added primes
double AddPrimesByNumber(FHEcontext& context, long nPrimes, 
			 long startAt=1,
			 bool special=false);

extern FHEcontext* activeContext; // Points to the "current" context

class FHEcontext {
  vector<Cmodulus> moduli;    // Cmodulus objects for the different primes
  // This is private since the implementation assumes that the list of
  // primes only grows and no prime is ever modified or removed.

  PlaintextSpace ptxtSpace;
public:
  // FHEContext is meant for convenience, not encapsulation: Most data
  // members are public and can be initialized by the application program.

  PAlgebra zMstar;       // The structure of Zm*
  //PAlgebraModTwo modTwo; // The structure of Z[X]/(Phi_m(X),2)
  //PAlgebraMod2r mod2r;   // The structure of Z[X]/(Phi_m(X),2^r)
  
  // The public encryption key and "fresh" ciphertexts are encrypted relative
  // to only a subset of the prime, to allow for mod-UP during key-switching.
  // In ctxtPrimes we keep the indexes of this subset. Namely, for a ciphertext
  // part p in a fresh ciphertext we have p.getMap().getIndexSet()==ctxtPrimes.
  // For convenience, we also keep in specialPrimes the complemeting subset,
  // i.e., specialPrimes = [0,numPrimes()-1] \setminus ctxtPrimes
  IndexSet ctxtPrimes;
  IndexSet specialPrimes;

  // The different columns in any key-switching matrix contain encryptions
  // of multiplies of the secret key, sk, B1*sk, B2*B1*sk, B3*B2*B1*sk,...
  // with each Bi a product of a few "non-special" primes in the chain. The
  // digits data member indicate which primes correspond to each of the Bi's.
  // These are all IndexSet objects, whose union is the subset ctxtPrimes.
  // The number of Bi's is one less than the number of columns in the key
  // switching matrices (since the 1st column encrypts sk, without any Bi's),
  // but we keep in the digits vector also an entry for the primes that do
  // not participate in any Bi (so digits.size() is the same as the number
  // of columns in the key switching matrices).

  vector<IndexSet> digits; // digits of ctxt/columns of key-switching matrix

  double stdev; // stdev is sqrt(variance) of the LWE error (default=3.2)

  ZZ modulusQ;
  unsigned logQ;

  unsigned decompSize; // number of bytes in each digit
  unsigned ndigits;

  unsigned primesNeeded;
  
 FHEcontext(unsigned m, unsigned logQ, unsigned p, 
            unsigned generator, unsigned decompSize = 3) {
    Init(m, logQ, to_ZZ(p), generator, decompSize);
  }

 FHEcontext(unsigned m, unsigned logQ, const ZZ &p, 
            unsigned generator, unsigned decompSize = 3) {
    Init(m, logQ, p, generator, decompSize);
  }
  
  FHEcontext(ifstream &in) {
    ImportSIContext(in);
  }
  
  void Init(unsigned m, unsigned logQ, const ZZ &p, unsigned generator, unsigned decompSize = 3) {
    stdev = 3.2;
    
    zMstar.init(m, generator);

    this->logQ = logQ;
    this->modulusQ = 1;
    this->modulusQ <<= logQ;

    this->decompSize = decompSize;
    this->ndigits = (logQ + 8*decompSize - 1) / (8*decompSize);
    
    ptxtSpace.Init(zMstar.PhimX(), p, generator);
  }
  
  void ExportSIContext(ofstream &out);
  void ImportSIContext(ifstream &in);
  
  void SetUpSIContext(long xi = 1);
  void SetUpBGVContext(unsigned L, unsigned c, long w, long pSize, long p0Size);
  
  unsigned Generator() const {
    return ptxtSpace.generator;
  }
  
  const ZZ &ModulusP() const {
    return ptxtSpace.p;
  }
  
  const PlaintextSpace &GetPlaintextSpace() const {
    return ptxtSpace;
  }
  
  void SetP(const ZZ &p) {
    ZZ_p::init(p);
    ptxtSpace.Init(zMstar.PhimX(), p);
  }
  
  long ithPrime(unsigned i) const {
    return (i<moduli.size()) ? moduli[i].getQ() : 0;
  }

  const Cmodulus& ithModulus(unsigned i) const { return moduli[i]; }

  long numPrimes() const { return moduli.size(); }

  // is num divisible by primes in the chain?
  bool isZeroDivisor(const ZZ& num) const {
    for (unsigned i=0; i<moduli.size(); i++) 
      if (divide(num,moduli[i].getQ())) return true;
    return false;
  }

  bool inChain(long p) const {
    for (unsigned i=0; i<moduli.size(); i++) 
      if (p==moduli[i].getQ()) return true;
    return false;
  }

  // The product of all the primes in the given set
  void productOfPrimes(ZZ& p, const IndexSet& s) const;

  ZZ productOfPrimes(const IndexSet& s) const {
    ZZ p;
    productOfPrimes(p,s);
    return p;
  }

  ZZ productOfPrimes() const {
    return productOfPrimes(ctxtPrimes);
  }

  // FIXME: run-time error when ithPrime(i) returns 0
  double logOfPrime(unsigned i) const { return log(ithPrime(i)); }

  // returns the natural log of productOfPrimes(s)
  double logOfProduct(const IndexSet& s) const {
    if (s.last() >= numPrimes())
      Error("FHEContext::logOfProduct: IndexSet has too many rows");

    double ans = 0.0;
    for (long i = s.first(); i <= s.last(); i = s.next(i))
      ans += logOfPrime(i);
    return ans;
  }

  void AddPrime(long p, bool special, long root = 0); 
  
  friend ostream &operator<<(ostream &os, const FHEcontext &context) {
    os << "logQ: " << context.logQ << endl
       << "p: " << context.ModulusP() << endl
       << "g: " << context.Generator() << endl
       << "primes: [";
    for (unsigned i = 0; i < context.moduli.size(); i++) {
      os << context.moduli[i].getQ() << ", ";
    }
    os << "]" << endl;
    
    return os;
  }
  
};

#endif
