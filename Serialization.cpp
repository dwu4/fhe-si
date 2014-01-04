#include "Serialization.h"

void Export(ofstream &out, const ZZ &val) {
  uint32_t nBytes = NumBytes(val);
  out.write((char *) &nBytes, sizeof(uint32_t));
  
  bool neg = (val < 0);
  out.write((char *) &neg, sizeof(bool));
  
  unsigned char data[nBytes];
  BytesFromZZ(data, val, nBytes);
  out.write((char *) data, nBytes);
}

void Import(ifstream &in, ZZ &val) {
  uint32_t nBytes;
  in.read((char *) &nBytes, sizeof(uint32_t));
  
  bool neg;
  in.read((char *) &neg, sizeof(bool));
  
  unsigned char data[nBytes];
  in.read((char *) data, nBytes);
  ZZFromBytes(val, data, nBytes);
  
  if (neg) val *= -1;
}

void Export(ofstream &out, const ZZX &poly) {
  int32_t degree = deg(poly);
  
  out.write((char *) &degree, sizeof(int32_t));
  for (int i = 0; i <= degree; i++) {
    Export(out, poly.rep[i]);
  }
}

void Import(ifstream &in, ZZX &poly) {
  poly = ZZX::zero();
  
  int32_t degree;
  
  in.read((char *) &degree, sizeof(int32_t));
  if (degree == -1) {
    return;
  }
  
  poly.SetMaxLength(degree + 1);
  for (int i = 0; i <= degree; i++) {
    ZZ coeff;
    Import(in, coeff);
    SetCoeff(poly, i, coeff);
  }
}

void Export(ofstream &out, const DoubleCRT &poly) {
  IndexMap<vec_long> map = poly.getMap();
  
  uint32_t size = map.getIndexSet().card();
  Export(out, size);
  for (long i = map.first(); i <= map.last(); i = map.next(i)) {
    Export(out, i);
    Export(out, map[i]);
  }
}

void Import(ifstream &in, DoubleCRT &poly) {
  IndexMap<vec_long> map;
  
  uint32_t size;
  Import(in, size);
  
  for (unsigned i = 0; i < size; i++) {
    long key;
    Import(in, key);
    map.insert(key);
    Import(in, map[key]);
  }
  
  poly.setMap(map);
}

void Export(ofstream &out, const vec_long &vec) {
  uint32_t len = vec.length();
  Export(out, len);
  for (int i = 0; i < vec.length(); i++) {
    Export(out, vec[i]);
  }
}

void Import(ifstream &in, vec_long &vec) {
  uint32_t size;
  Import(in, size);
  vec.SetLength(size);
  
  for (uint32_t i = 0; i < size; i++) {
    Import(in, vec[i]);
  }
}

void Export(ofstream &out, const CiphertextPart &part) {
  Export(out, part.poly);
}

void Import(ifstream &in, CiphertextPart &part) {
  Import(in, part.poly);
}

void Export(ofstream &out, const Ciphertext &ctxt) {
  Ciphertext copy = ctxt;
  
  copy.ScaleDown();
  Export(out, copy.parts);
}

void Import(ifstream &in, Ciphertext &ctxt) {
  ctxt.Clear();
  Import(in, ctxt.parts);
}
