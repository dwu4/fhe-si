#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#include "FHEContext.h"
#include <vector>

#include "FHE-SI.h"
#include "Ciphertext.h"
#include "Matrix.h"
#include "Matrix.cpp"

class Statistics {
 public:
  Statistics(const FHEcontext &context) : context(context), secretKey(context), 
                                          publicKey(secretKey), keySwitch(secretKey),
                                          data(Ciphertext(publicKey)) {
    int k = context.Generator();
    
    int nSlots = context.GetPlaintextSpace().GetUsableSlots();
    while (nSlots > 1) {
      autoKeySwitch.push_back(KeySwitchSI(secretKey, k));
    
      nSlots >>= 1;
      k *= k;
      k %= context.zMstar.M();
    }
  }
  
  void AddData(const Matrix<Plaintext> &blocks, const vector<Plaintext> &blockSizes) {
    for (unsigned i = 0; i < blocks.NumRows(); i++) {
      vector<Ciphertext> encExample(blocks[i].size(), Ciphertext(publicKey));
      for (unsigned j = 0; j < blocks[i].size(); j++) {
        publicKey.Encrypt(encExample[j], blocks[i][j]);
      }
      Ciphertext encN(publicKey);
      publicKey.Encrypt(encN, blockSizes[i]);
      
      data.AddRow(encExample);
      nElems.push_back(encN); 
    }
  }

  void Clear() {
    data.Clear();
    nElems.clear();
  }

  void ComputeNthMoment(vector<Ciphertext> &moment, Ciphertext &denom, unsigned n) {
    if (n < 1 || n > 2) { // Not supported currently
      return;
    }
    
    moment.resize(data.NumCols(), Ciphertext(publicKey));
    denom = nElems[0];
    for (unsigned j = 0; j < data.NumCols(); j++) {
      moment[j] = data(0, j);
      if (n == 2) {
        moment[j] *= moment[j];
      }
      
      for (unsigned i = 1; i < data.NumRows(); i++) {
        if (j == 0) {
          denom += nElems[i];
        }
        
        Ciphertext tmp = data(i, j);
        if (n == 2) {
          tmp *= tmp;
        }
        moment[j] += tmp;
      }
      
      if (n == 2) {
        keySwitch.ApplyKeySwitch(moment[j]);
      }
      SumBatchedData(moment[j]);
    }
    
    Ciphertext noise(publicKey);
    for (unsigned i = 0; i < moment.size(); i++) {
      GenerateNoise(noise);
      moment[i] += noise;
    }
  }
  
  void ComputeCovariance(Matrix<Ciphertext> &cov, vector<Ciphertext> &mu, Ciphertext &n, Ciphertext &n2) {
    ComputeNthMoment(mu, n, 1);
    
    Ciphertext dummy(publicKey);
    Matrix<Ciphertext> muMat(dummy);
    muMat.AddRow(mu);
    muMat.Transpose();
    muMat.MultByTranspose();
    
    // Symmetric matrix so just process half the entries
    for (unsigned i = 0; i < muMat.NumRows(); i++) {
      for (unsigned j = i; j < muMat.NumCols(); j++) {
        keySwitch.ApplyKeySwitch(muMat(i,j));
        muMat(i,j) *= -1;
      }
    }
    
    cov = data;
    cov.Transpose();
    cov.MultByTranspose();
    
    Ciphertext noise(publicKey);
    // Again, symmetric matrix so we just process half the entries
    for (unsigned i = 0; i < cov.NumRows(); i++) {
      for (unsigned j = i; j < cov.NumCols(); j++) {
        keySwitch.ApplyKeySwitch(cov(i,j));
        SumBatchedData(cov(i,j));
        cov(i,j) *= n;
        keySwitch.ApplyKeySwitch(cov(i,j));
        cov(i,j) += muMat(i,j);
        
        GenerateNoise(noise);
        cov(i,j) += noise;
        
        cov(j,i) = cov(i,j);
      }
    }

    n2 = n;
    n2 *= n2;
    
    keySwitch.ApplyKeySwitch(n2);
  }
  
  // For debugging purposes
  FHESISecKey &GetSecretKey() { return secretKey; }
  FHESIPubKey &GetPublicKey() { return publicKey; }
  
 private:
  const FHEcontext &context;
 
  FHESISecKey secretKey;
  FHESIPubKey publicKey;
  
  KeySwitchSI keySwitch;
  vector<KeySwitchSI> autoKeySwitch;

  Matrix<Ciphertext> data;
  vector<Ciphertext> nElems;
  
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

void ComputeNthMomentPT(vector<ZZ> &moment, unsigned n, const Matrix<ZZ> &data) {
  moment.resize(data.NumCols());
  
  for (unsigned j = 0; j < data.NumCols(); j++) {
    moment[j] = 0;
    ZZ val;
    for (unsigned i = 0; i < data.NumRows(); i++) {
      power(val, data(i,j), n);
      moment[j] += val;
    }
  }
}

void ComputeMomentsPT(vector<ZZ> &sum, vector<ZZ> &sqSum, const Matrix<ZZ> &data) {
  ComputeNthMomentPT(sum, 1, data);
  ComputeNthMomentPT(sqSum, 2, data);
}

void ComputeCovariancePT(Matrix<ZZ> &cov, const Matrix<ZZ> &data) {
  cov = data;
  cov.Transpose();
  cov.MultByTranspose();
  
  ZZ n = to_ZZ(data.NumRows());
  cov *= n;
  
  vector<ZZ> muVec;
  ComputeNthMomentPT(muVec, 1, data);
  
  Matrix<ZZ> mu;
  mu.AddRow(muVec);
  mu.Transpose();
  mu.MultByTranspose();
  
  cov -= mu;
}

#endif
