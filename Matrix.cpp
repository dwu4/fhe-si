#ifndef _MATRIX_CPP
#define _MATRIX_CPP

#include "Matrix.h"

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &other) const {
  Matrix<T> newMatrix = *this;
  newMatrix += other;
  return newMatrix;
}

template<class T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &other) {
  for (unsigned i = 0; i < NumRows(); i++) {
    for (unsigned j = 0; j < NumCols(); j++) {
      ElemAt(i,j) += other(i,j);
    }
  }
  return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &other) const {
  Matrix<T> newMatrix = *this;
  newMatrix -= other;
  return newMatrix;
}

template<class T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &other) {
  for (unsigned i = 0; i < NumRows(); i++) {
    for (unsigned j = 0; j < NumCols(); j++) {
      T tmp = other(i,j);
      tmp *= -1;
      ElemAt(i,j) += tmp;
    }
  }
  return *this;
}

template<class T>
Matrix<T> Matrix<T>::operator*(Matrix<T> &other) const {
  Matrix<T> newMatrix = *this;
  newMatrix *= other;
  return newMatrix;
}

template<class T>
Matrix<T> Matrix<T>::operator*(vector<T> &other) const {
  Matrix<T> newMatrix = *this;
  newMatrix *= other;
  return newMatrix;
}

template<class T>
Matrix<T> &Matrix<T>::operator*=(Matrix &other) {
  if (mat.empty()) return *this;
  Matrix newMatrix(NumRows(), other.NumCols(), dummy);
  
  for (unsigned i = 0; i < NumRows(); i++) {
    for (unsigned j = 0; j < other.NumCols(); j++) {
      newMatrix(i,j) = ElemAt(i,0);
      newMatrix(i,j) *= other(0,j);
      
      for (unsigned k = 1; k < NumCols(); k++) {
        T tmp = ElemAt(i, k);
        tmp *= other(k, j);
        
        newMatrix(i,j) += tmp;
      }
    }
  }
  
  swap(newMatrix.mat, mat);
  transpose = false;
  return *this;
}

template<class T>
Matrix<T> &Matrix<T>::operator*=(vector<T> &other) {
  if (mat.empty()) return *this;
  Matrix newMatrix(NumRows(), 1, dummy);
  
  for (unsigned i = 0; i < NumRows(); i++) {
    ElemAt(i,0) *= other[0];
    newMatrix(i,0) = ElemAt(i,0);
    for (unsigned j = 1; j < NumCols(); j++) {
      ElemAt(i,j) *= other[j];
      newMatrix(i,0) += ElemAt(i,j);
    }
  }
  
  swap(mat, newMatrix.mat);
  transpose = false;
  return *this;
}

template<class T>
Matrix<T> &Matrix<T>::operator*=(T &other) {
  if (mat.empty()) return *this;
  for (unsigned i = 0; i < NumRows(); i++) {
    for (unsigned j = 0; j < NumCols(); j++) {
      ElemAt(i,j) *= other;
    }
  }
  return *this;
}

template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &other) {
  mat = other.mat;
  transpose = other.transpose;
  dummy = other.dummy;
  
  return *this;
}

template<class T>
T &Matrix<T>::operator()(unsigned row, unsigned col) {
  return ElemAt(row, col);
}

template<class T>
T const &Matrix<T>::operator()(unsigned row, unsigned col) const {
  return ElemAt(row, col);
}

template<class T>
vector<T> &Matrix<T>::operator[](unsigned row) {
  return mat[row];
}

template<class T>
vector<T> const &Matrix<T>::operator[](unsigned row) const {
  return mat[row];
}

template<class T>
T &Matrix<T>::ElemAt(unsigned row, unsigned col) {
  return transpose ? mat[col][row] : mat[row][col];
}

template<class T>
T const &Matrix<T>::ElemAt(unsigned row, unsigned col) const {
  return transpose ? mat[col][row] : mat[row][col];
}

template<class T>
void Matrix<T>::MultByTranspose() {
  if (mat.empty()) return;
  Matrix newMatrix(NumRows(), NumRows(), dummy);
  
  for (unsigned i = 0; i < NumRows(); i++) {
    for (unsigned j = i; j < NumRows(); j++) {
      newMatrix(i,j) = ElemAt(i,0);
      newMatrix(i,j) *= ElemAt(j,0);
      
      for (unsigned k = 1; k < NumCols(); k++) {
        T tmp = ElemAt(i, k);
        tmp *= ElemAt(j, k);
        
        newMatrix(i,j) += tmp;
      }
      
      if (i != j) {
        newMatrix(j,i) = newMatrix(i,j);
      }
    }
  }
  
  swap(newMatrix.mat, mat);
  transpose = false;
}

template<class T>
void Matrix<T>::Transpose() {
  transpose = !transpose;
}

template<class T>
void Matrix<T>::Invert(T& det, std::function<void(T&)> reduce) {
  unsigned dim = NumRows();
  
  Matrix<T> adj(dim, dim, dummy);
  vector<bool> usedRows(dim), usedCols(dim);
  for (unsigned i = 0; i < dim; i++) {
    for (unsigned j = 0; j < dim; j++) {
      usedRows[i] = usedCols[j] = true;
      Determinant(adj(j,i), usedRows, usedCols, dim-1, reduce);
      usedRows[i] = usedCols[j] = false;
      if ((i+j) % 2 == 1) {
        adj(j,i) *= -1;
      }
    }
  }
  
  // Since we already computed the adjugate matrix, we can
  // reuse the results to compute the determinant  
  det = ElemAt(0,0);
  det *= adj(0,0);
  for (unsigned i = 1; i < dim; i++) {
    T tmp = ElemAt(0,i);

    tmp *= adj(i,0);
    det += tmp;
  }
  if (reduce) {
    reduce(det);
  }
  
  swap(adj.mat, mat);
  transpose = false;
}

template<class T>
void Matrix<T>::Determinant(T &det, std::function<void(T&)> reduce) const {
  unsigned dim = NumRows();
  vector<bool> usedRows(dim), usedCols(dim);
  Determinant(det, usedRows, usedCols, dim, reduce);
}

template<class T>
void Matrix<T>::Determinant(T &det, vector<bool> &usedRows, 
                            vector<bool> &usedCols, unsigned dim,
                            std::function<void(T&)> reduce) const {
  unsigned matDim = NumRows();
  unsigned row = 0;
  while (usedRows[row]) row++;
  
  bool negative = false;
  bool first = true;
  for (unsigned col = 0; col < matDim; col++) {
    if (usedCols[col]) continue;
 
    if (dim == 1) {
      det = ElemAt(row, col);
      return;
    }
 
    T tmp = ElemAt(row, col);
    if (negative) tmp *= -1;
    negative = !negative;
    
    usedRows[row] = usedCols[col] = true;
    T tmp2(dummy);
    Determinant(tmp2, usedRows, usedCols, dim-1, reduce);
    usedRows[row] = usedCols[col] = false;
    
    tmp *= tmp2;
    
    if (first) {
      det = tmp;
      first = false;
    } else {
      det += tmp;
    }
  }
  
  if (reduce) {
    reduce(det);
  }
}

template<class T>
unsigned Matrix<T>::NumRows() const {
  if (mat.empty()) return 0;
  return transpose ? mat[0].size() : mat.size();
}

template<class T>
unsigned Matrix<T>::NumCols() const {
  if (mat.empty()) return 0;
  return transpose ? mat.size() : mat[0].size();
}

template<class T>
void Matrix<T>::Resize(unsigned nRows, unsigned nCols) {
  Clear();
  transpose = false;
  
  mat.resize(nRows);
  for (unsigned i = 0; i < nRows; i++) {
    mat[i].resize(nCols, dummy);
  }
}

template<class T>
void Matrix<T>::AddRow(vector<T> &row) {
  if (transpose) return; // No support for adding to transposed matrix
  mat.push_back(row);
}

template<class T>
void Matrix<T>::Concatenate(Matrix<T> &other) {
  if (transpose) return; // No support for concatenating to transposed matrix
  mat.insert(mat.end(), other.mat.begin(), other.mat.end());
}

template<class T>
void Matrix<T>::Clear() {
  mat.clear();
}

template<class T>
void Matrix<T>::MapAll(std::function<void(T&)> func) {
  for (unsigned i = 0; i < mat.size(); i++) {
    for (unsigned j = 0; j < mat[i].size(); j++) {
      func(mat[i][j]);
    }
  }
}

template<class T>
ostream &operator<<(ostream &os, const Matrix<T> &m) {
  for (unsigned i = 0; i < m.NumRows(); i++) {
    for (unsigned j = 0; j < m.NumCols(); j++) {
      os << m(i,j) << " ";
    }
    os << endl;
  }
  return os;
}

#endif
