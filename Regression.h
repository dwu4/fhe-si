#ifndef _REGRESSION_H_
#define _REGRESSION_H_

#include "FHEContext.h"
#include "Plaintext.h"
#include "FHE-SI.h"
#include "Ciphertext.h"
#include <vector>
#include <fstream>

#include "Matrix.h"
#include "Matrix.cpp"

bool LoadData(Matrix<ZZ> &rawData, vector<ZZ> &labels,
              unsigned &dim, const string &filename) {
  ifstream fin;
  fin.open(filename);
  if (!fin) {
    cout << "Unable to read data file." << endl;
    return false;
  }
  
  rawData.Clear();
  labels.clear();

  int label, n;

  fin >> dim >> n;
  vector<ZZ> data(dim);
  for (int i = 0; i < n; i++) {
    for (unsigned j = 0; j < dim; j++) {
      fin >> data[j];
    }
    fin >> label;

    rawData.AddRow(data);
    labels.push_back(to_ZZ(label));
  }
  
  return true;
}

double BatchData(vector<vector<Plaintext>> &ptxtData, vector<Plaintext> &ptxtLabels,
             const Matrix<ZZ> &rawData, const vector<ZZ> &labels, const FHEcontext &context) {
  double start(clock());
  
  ZZ p = to_ZZ(context.ModulusP());
  ptxtData.clear();
  ptxtLabels.clear();
  unsigned batchSize = context.GetPlaintextSpace().GetUsableSlots();
  for (unsigned i = 0; i < rawData.NumRows(); i += batchSize) {
    vector<Plaintext> row;
    vector<ZZ> curLabel;
    for (unsigned j = 0; j < rawData.NumCols(); j++) {
      vector<ZZ> columnBatch;
      for (unsigned k = i; ((k < i + batchSize) && (k < rawData.NumRows())); k++) {
        columnBatch.push_back(to_ZZ(rawData[k][j]) % p);
        if(j==0) curLabel.push_back(to_ZZ(labels[k]) % p);
      }
      row.push_back(Plaintext(context, columnBatch));
    }
    ptxtData.push_back(row);
    ptxtLabels.push_back(Plaintext(context, curLabel));
  }
  return (clock()-start)/CLOCKS_PER_SEC;
}

class Regression {
 public:
  Regression(const FHEcontext &context) : context(context), secretKey(context), publicKey(secretKey), keySwitch(secretKey), data(Ciphertext(publicKey)) {
    int k = context.Generator();
    
    int nSlots = context.GetPlaintextSpace().GetUsableSlots();
    while (nSlots > 1) {
      autoKeySwitch.push_back(KeySwitchSI(secretKey, k));
    
      nSlots >>= 1;
      k *= k;
      k %= context.zMstar.M();
    }
  }
  
  void AddData(const vector<vector<Plaintext>> &ptxtData, const vector<Plaintext> &ptxtLabels) {
    for (unsigned i = 0; i < ptxtData.size(); i++) {
      vector<Ciphertext> encExample(ptxtData[i].size(), Ciphertext(publicKey));
      for (unsigned j = 0; j < ptxtData[i].size(); j++) {
        publicKey.Encrypt(encExample[j], ptxtData[i][j]);
      }
      Ciphertext encLabel(publicKey);
      publicKey.Encrypt(encLabel, ptxtLabels[i]);
      
      data.AddRow(encExample);
      labels.push_back(encLabel); 
    }
  }

  void Clear() {
    data.Clear();
    labels.clear();
  }

  void Regress(vector<Ciphertext> &theta, Ciphertext &det) const {
    Matrix<Ciphertext> dataCopy = data;
    vector<Ciphertext> labels = this->labels;
    dataCopy.Transpose();

    Matrix<Ciphertext> last = dataCopy*labels;
    dataCopy.MultByTranspose();
    
    auto processFunc = [this](Ciphertext &ctxt) {
      this->keySwitch.ApplyKeySwitch(ctxt);
      this->SumBatchedData(ctxt);
    };
    last.MapAll(processFunc); 
    dataCopy.MapAll(processFunc);

    // Hack for when dimension is 1 (so we don't need to provide
    // an encryption of 1)
    if (data.NumCols() == 1) {
      det = dataCopy(0,0);
      theta.assign(1, last(0,0));
      
      return;
    }

    dataCopy.Invert(det, [this](Ciphertext &ctxt) {
      this->keySwitch.ApplyKeySwitch(ctxt);
    });
    
    dataCopy *= last;
    dataCopy.MapAll([this](Ciphertext &ctxt) {
      this->keySwitch.ApplyKeySwitch(ctxt);
    });
    
    theta.resize(dataCopy.NumRows(), Ciphertext(publicKey));
    for (unsigned i = 0; i < dataCopy.NumRows(); i++) {
      theta[i] = dataCopy(i,0);
    }
    
    // Add uniformly random values to all but the first slot
    Ciphertext noise(publicKey);
    for (unsigned i = 0; i < theta.size(); i++) {
      GenerateNoise(noise);
      theta[i] += noise;
    }
    
    GenerateNoise(noise);
    det += noise;
  }
  
  FHESIPubKey &GetPublicKey() { return publicKey; }
  FHESISecKey &GetSecretKey() { return secretKey; }
  
  vector<Ciphertext> labels;
 private:
  const FHEcontext &context;
 
  FHESISecKey secretKey;
  FHESIPubKey publicKey;
  
  KeySwitchSI keySwitch;
  vector<KeySwitchSI> autoKeySwitch;

  Matrix<Ciphertext> data;
  
  void SumBatchedData(Ciphertext &batchedData) const {
    int k = context.Generator();
    
    for (unsigned i = 0; i < autoKeySwitch.size(); i++) {
      Ciphertext tmp = batchedData;
      tmp >>= k;
      autoKeySwitch[i].ApplyKeySwitch(tmp);
      batchedData += tmp;
      
      k *= k;
      k %= context.zMstar.M();
    }
  }
  
  void GenerateNoise(Ciphertext &noiseCtxt) const {
    vector<ZZ_pX> randVec(context.GetPlaintextSpace().GetTotalSlots());
    randVec[0] = to_ZZ_pX(0);
    for (unsigned i = 1; i < randVec.size(); i++) {
      randVec[i] = to_ZZ_pX(random_ZZ_p());
    }
    
    Plaintext ptxt(context);
    ptxt.EmbedInSlots(randVec, false);
    publicKey.Encrypt(noiseCtxt, ptxt);
  }
};

void RegressPT(vector<ZZ> &theta, ZZ &det, Matrix<ZZ> &data,
               vector<ZZ> &labels) {
  Matrix<ZZ> A = data;
  A.Transpose();
  
  Matrix<ZZ> tmp = A*labels;

  A.MultByTranspose();
  if (data.NumCols() == 1) {
    det = A(0,0);
    theta.assign(1, tmp(0,0));
    return;
  }
  
  A.Invert(det);
  A *= tmp;
  
  theta.resize(A.NumRows());
  for (unsigned i = 0; i < A.NumRows(); i++) {
    theta[i] = A(i,0);
  }
}

#endif
