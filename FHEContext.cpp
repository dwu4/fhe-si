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
 
#include "FHEContext.h"
#include "Serialization.h"

// A global variable, pointing to the "current" context
FHEcontext* activeContext = NULL;

void FHEcontext::productOfPrimes(ZZ& p, const IndexSet& s) const {
  p = 1;
  for (long i = s.first(); i <= s.last(); i = s.next(i))
    p *= ithPrime(i);

}

void FHEcontext::AddPrime(long p, bool special, long root) {
  long twoM = 2 * zMstar.M();  // ensure p-1 is divisible by 2m
  
  //if (special) twoM *= context.ModulusP();
  assert( ProbPrime(p) && p % twoM == 1 && !inChain(p) );

  long i = moduli.size();
  moduli.push_back( Cmodulus(zMstar, p, root) );

  if (special)
    specialPrimes.insert(i);
  else
    ctxtPrimes.insert(i);
}

void FHEcontext::ExportSIContext(ofstream &out) {
  Export(out, zMstar.M());
  Export(out, logQ);
  Export(out, ModulusP());
  Export(out, Generator());
  Export(out, decompSize);
  
  uint32_t size = moduli.size();
  Export(out, size);
  for (unsigned i = 0; i < moduli.size(); i++) {
    long q = moduli[i].getQ();
    long root = moduli[i].getRoot();
    Export(out, q);
    Export(out, root);
  }
}

void FHEcontext::ImportSIContext(ifstream &in) {
  unsigned m, logQ, generator, decompSize;
  ZZ p;
  
  Import(in, m);
  Import(in, logQ);
  Import(in, p);
  Import(in, generator);
  Import(in, decompSize);
  
  Init(m, logQ, p, generator, decompSize);
  uint32_t size;
  Import(in, size);
  long q, root;
  for (unsigned i = 0; i < size; i++) {
    Import(in, q);
    Import(in, root);
    AddPrime(q, false, root);
  }
}

void FHEcontext::SetUpSIContext(long xi) {
  AddPrimesBySize(*this, log(modulusQ)*2+log(ModulusP())+log(zMstar.phiM())*2+log(2)+log(xi), false);
}

// Adds to the chain primes whose product is at least totalSize bits
double AddPrimesBySize(FHEcontext& context, double totalSize, bool special) {
  if (!context.zMstar.M() || context.zMstar.M() > (1<<20)) // sanity checks
    Error("AddModuli1: m undefined or larger than 2^20");

  long p = (1UL << NTL_SP_NBITS)-1;   // Start from as large prime as possible
  long twoM = 2 * context.zMstar.M(); // make p-1 divisible by 2m
  //if(special) twoM *= context.ModulusP();  // ensure congruent to 1 mod ptxtSpace
  p -= (p%twoM); // 0 mod 2m
  p += twoM +1;  // 1 mod 2m, a 2m quantity is subtracted below

  bool lastPrime = false;
  double sizeLeft = totalSize; // how much is left to do
  while (sizeLeft > 0.0) {
    // A bit of convoluted logic attempting not to overshoot totalSize by much
    if (sizeLeft < log((double)p) && !lastPrime) { // decrease p
      lastPrime = true;
      p = ceil(exp(sizeLeft));
      p -= (p%twoM)-1; // make p-1 divisible by 2m
      twoM = -twoM;    // increase p below, rather than decreasing it
    }
    do { p -= twoM; } while (!ProbPrime(p)); // next prime
    if (!context.inChain(p)) {
      context.AddPrime(p, special);
      sizeLeft -= log((double)p);
    }
  }
  return totalSize-sizeLeft;
}


// Adds nPrimes primes to the chain, returns the bitsize of the product of
// all primes in the chain.
double AddPrimesByNumber(FHEcontext& context, long nPrimes, 
			 long p, bool special) {
  if (!context.zMstar.M() || context.zMstar.M() > (1<<20))  // sanity checks
    Error("FHEcontext::AddModuli2: m undefined or larger than 2^20");

  long twoM = 2 * context.zMstar.M();

  // make sure that p>0 and that p-1 is divisible by m
  if (p<1) p = 1;
  p -= (p % twoM) -1;

  double sizeSoFar = 0.0;
  while (nPrimes>0) {
    do { p += twoM; } while (!ProbPrime(p)); // next prime
    if (!context.inChain(p)) {
      context.AddPrime(p, special);
      nPrimes -= 1;
      sizeSoFar += log((double)p);
    }
  }
  return sizeSoFar;
}
