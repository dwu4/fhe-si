#ifndef _UTIL_H_
#define _UTIL_H_

#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>
#include <vector>
#include <cstring>
#include <iostream>

using namespace NTL;
using namespace std;

void Reduce(ZZ &val, unsigned logQ, bool positive = false);
void ReduceCoefficients(ZZX &poly, unsigned logQ,
          bool positive = false);
void ReduceCoefficientsSlow(ZZX &poly, const ZZ &modulus,
			    bool positive = false);
void ReduceCoefficientsSlow(ZZX &poly, unsigned modulus,
			    bool positive = false);


/*
  Reduces val to integer congruent (mod q)  in the range (-q/2,q/2]
  Assumes q is a power of two.
*/
void ReduceMod(ZZ &val, unsigned logQ);

/*
  Base decomposition for bases of the form 2^(8*size)
*/
void ByteDecomp(vector<ZZX> &smallPolys, const ZZX &poly,
		unsigned logQ, unsigned size = 1);

void ByteDecomp(vector<vector<ZZX>> &smallPolys, const vector<ZZX> &polys,
		unsigned logQ, unsigned size = 1);

/*
  First n powers of 2^(8*i) times poly
*/
void Powers(vector<ZZX> &powers, const ZZX &poly, long n, unsigned size = 1);

/*
  byteDecomps is the vector of Byte Decomps to match components size with.
*/
void Powers(vector<vector<ZZX>> &powers, const vector<ZZX> &polys,
	    vector<vector<ZZX>> byteDecomps, unsigned size = 1);

void SampleRandom(ZZX &poly, const ZZ &modulus, unsigned deg);

void GetConstantTerm(ZZ &res, ZZ_pX &poly);

template<typename T>
void PrintVector(const vector<T> &vec, ostream &out = std::cout) {
  for (unsigned i = 0; i < vec.size(); i++) {
    out << vec[i] << " ";
  }
}

template<typename T>
void PrintVector(const vector<vector<T>> &vec, ostream &out = std::cout) {
  for (unsigned i = 0; i < vec.size(); i++) {
    PrintVector(vec[i]);
    out << endl;
  }
}

template <typename T>
unsigned ComputeLog(T val) {
  unsigned log = 0;
  while (val != 0) {
    val >>= 1;
    log++;
  }

  return log-1;
}

template<typename T>
static void DotProduct(T &res, const vector<T> &v1,
		       const vector<T> &v2) {
  #ifdef DEBUG
  assert(v1.size() == v2.size());
  #endif

  if (v1.size() == 0) {
    return;
  }

  res = v1[0];
  res *= v2[0];

  for (unsigned i = 1; i < v1.size(); i++) {
    T val = v1[i];
    val *= v2[i];
    res += val;
  }
}

template<typename T>
static void TensorProduct(vector<T> &res, const vector<T> &v1, const vector<T> &v2) {
  res.resize(v1.size()*v2.size());
  
  unsigned ind = 0;
  for (unsigned i = 0; i < v1.size(); i++) {
    for (unsigned j = 0; j < v2.size(); j++) {
      res[ind] = v1[i];
      res[ind++] *= v2[j];
    }
  }
}

#endif
