#ifndef _MATRIX_H
#define _MATRIX_H

#include <vector>
#include <iostream>
#include <functional>

using namespace std;

template<class T>
class Matrix;

template <typename T>
ostream& operator<<(ostream &, const Matrix<T> &);

template<class T>
class Matrix {
  T dummy;
  public:
    Matrix() : dummy(T()) {
      transpose = false;
    }
    Matrix(const T &dummy) : dummy(dummy) {
      transpose = false;
    }
    
    Matrix(unsigned nRows, unsigned nCols, const T &dummy) : dummy(dummy) {
      Resize(nRows, nCols);
    }
    
    Matrix(unsigned nRows, unsigned nCols) : dummy(T()) {
      Resize(nRows, nCols);
    }
    
    Matrix operator+(const Matrix &other) const;
    Matrix &operator+=(const Matrix &other);
    
    Matrix operator-(const Matrix &other) const;
    Matrix &operator-=(const Matrix &other);
    
    Matrix operator*(Matrix &other) const;
    Matrix operator*(vector<T> &other) const;
    Matrix &operator*=(Matrix &other);
    Matrix &operator*=(vector<T> &other);
    Matrix &operator*=(T &other);
    
    Matrix &operator=(const Matrix &other);
    
    T &operator()(unsigned row, unsigned col);
    T const &operator()(unsigned row, unsigned col) const;
    
    vector<T> &operator[](unsigned row);
    vector<T> const &operator[](unsigned row) const;
    
    void MultByTranspose();
    void Transpose();
    void Invert(T& det, std::function<void(T&)> reduce = NULL);
    
    void Determinant(T& det, std::function<void(T&)> reduce = NULL) const;
    
    unsigned NumRows() const;
    unsigned NumCols() const;
    
    void Resize(unsigned nRows, unsigned nCols);
    
    void AddRow(vector<T> &row);
    void Concatenate(Matrix<T> &other);
    void Clear();
    
    void MapAll(std::function<void(T&)> func);
    
    friend ostream &operator<< <>(ostream &os, const Matrix &m);
  private:
    vector<vector<T>> mat;
    bool transpose;
    
    T &ElemAt(unsigned row, unsigned col);
    T const &ElemAt(unsigned row, unsigned col) const;
    
    void Determinant(T& det, vector<bool> &usedRows,
                     vector<bool> &usedCols, unsigned dim,
                     std::function<void(T&)> reduce = NULL) const;
};

#endif
