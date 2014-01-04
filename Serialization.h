#ifndef _SERIALIZATION_H_
#define _SERIALIZATION_H_

#include <NTL/ZZ.h>
#include <fstream>

#include "DoubleCRT.h"
#include "Ciphertext.h"
#include "Matrix.h"

void Export(ofstream &out, const ZZ &val);
void Export(ofstream &out, const ZZX &val);
void Export(ofstream &out, const DoubleCRT &val);
void Export(ofstream &out, const vec_long &vec);

void Import(ifstream &in, ZZ &val);
void Import(ifstream &in, ZZX &val);
void Import(ifstream &in, DoubleCRT &val);
void Import(ifstream &in, vec_long &vec);

void Export(ofstream &out, const CiphertextPart &part);
void Export(ofstream &out, const Ciphertext &ctxt);

void Import(ifstream &in, CiphertextPart &part);
void Import(ifstream &in, Ciphertext &ctxt);

template<typename T>
void Export(ofstream &out, const T &val) {
  out.write((char *) &val, sizeof(T));
}

template<typename T>
void Import(ifstream &in, T &val) {
  in.read((char *) &val, sizeof(T));
}

template<typename T>
void Export(ofstream &out, const vector<T> &vec) {
  uint32_t size = vec.size();
  Export(out, size);
  for (unsigned i = 0; i < vec.size(); i++) {
    Export(out, vec[i]);
  }
}

template<typename T>
void Import(ifstream &in, vector<T> &vec) {
  uint32_t size;
  Import(in, size);
  
  vec.resize(size);
  for (unsigned i = 0; i < size; i++) {
    Import(in, vec[i]);
  }
}

template<typename T>
void Export(ofstream &out, const Matrix<T> &mat) {
  Export(out, mat.NumRows());
  Export(out, mat.NumCols());
  
  for (unsigned i = 0; i < mat.NumRows(); i++) {
    for (unsigned j = 0; j < mat.NumCols(); j++) {
      Export(out, mat(i,j));
    }
  }
}

template<typename T>
void Import(ifstream &in, Matrix<T> &mat) {
  uint32_t nRows, nCols;
  Import(in, nRows);
  Import(in, nCols);
  
  mat.Resize(nRows, nCols);
  for (unsigned i = 0; i < nRows; i++) {
    for (unsigned j = 0; j < nCols; j++) {
      Import(in, mat(i,j));
    }
  }
}

#endif