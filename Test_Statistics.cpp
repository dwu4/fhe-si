#include "FHE-SI.h"

#include "Statistics.h"
#include <time.h>
#include <ctime>
#include <fstream>
#include "Util.h"

static bool LoadData(Matrix<ZZ> &data, unsigned &dim, const string &filename) {
  ifstream fin;
  fin.open(filename);
  if (!fin) {
    cout << "Unable to read data file." << endl;
    return false;
  }
  
  data.Clear();

  unsigned n, junk;
  fin >> dim >> n;
  vector<ZZ> vec(dim);
  
  for (unsigned i = 0; i < n; i++) {
    for (unsigned j = 0; j < dim; j++) {
      fin >> vec[j];
    }
    fin >> junk;
    
    data.AddRow(vec);
  }
  
  return true;
}

double BatchData(Matrix<Plaintext> &ptxtBlocks, vector<Plaintext> &ptxtBlockSizes,
                 const Matrix<ZZ> &rawData, const FHEcontext &context) {
                 
  double start(clock());
  
  ZZ p = to_ZZ(context.ModulusP());
  unsigned batchSize = context.GetPlaintextSpace().GetUsableSlots();
  
  ptxtBlocks.Clear();
  ptxtBlockSizes.clear();
  
  for (unsigned i = 0; i < rawData.NumRows(); i += batchSize) {
    vector<Plaintext> row;
    unsigned n = 0;
    
    for (unsigned j = 0; j < rawData.NumCols(); j++) {
      vector<ZZ> columnBatch;
      for (unsigned k = i; ((k < i + batchSize) && (k < rawData.NumRows())); k++) {
        columnBatch.push_back(to_ZZ(rawData[k][j]) % p);
      }
      n = columnBatch.size();
      row.push_back(Plaintext(context, columnBatch));
    }
    
    ptxtBlocks.AddRow(row);
    ptxtBlockSizes.push_back(Plaintext(context, n));
  }
  
  return (clock()-start)/CLOCKS_PER_SEC;
}

static int RunStatisticsTest(Matrix<ZZ> &data, unsigned p, const FHEcontext &context) {
  vector<ZZ> mean;
  Matrix<ZZ> cov;
  
  ComputeNthMomentPT(mean, 1, data);
  ComputeCovariancePT(cov, data);
  
  cout << "True values: " << endl;
  cout << "  Mean: ";
  for (unsigned i = 0; i < mean.size(); i++) {
    cout << mean[i] << ", ";
  }
  cout << endl;
  
  cout << "  Covariance: " << endl;
  for (unsigned i = 0; i < cov.NumRows(); i++) {
    cout << "  ";
    for (unsigned j = 0; j < cov.NumCols(); j++) {
      cout << cov(i,j) << " ";
    }
    cout << endl;
  }
  cout << endl;
  
  cout << "Expected values: " << endl;
  cout << "  Mean: ";
  for (unsigned i = 0; i < mean.size(); i++) {
    cout << mean[i] % to_ZZ(p) << ", ";
  }
  cout << endl;
  
  cout << "  N: " << data.NumRows() % p << endl << endl;
  cov.MapAll([p](ZZ &elem) {
    elem %= to_ZZ(p);
  });
  cout << "  Covariance: " << endl;
  for (unsigned i = 0; i < cov.NumRows(); i++) {
    cout << "  ";
    for (unsigned j = 0; j < cov.NumCols(); j++) {
      cout << cov(i,j) % to_ZZ(p) << " ";
    }
    cout << endl;
  }
  cout << endl;
  cout << "  N^2: " << ((data.NumRows() % p) * (data.NumRows() % p)) % p << endl << endl;
  
  double start(clock());
  
  Statistics stats(context);
  cout << "Setup time: " << (clock()-start)/CLOCKS_PER_SEC << endl;
  
  Matrix<Plaintext> ptxtBlocks;
  vector<Plaintext> ptxtBlockSizes;
  
  double batchTime = BatchData(ptxtBlocks, ptxtBlockSizes, data, context);
  
  cout << "Batch time: " << batchTime << endl;

  double encryptionStart(clock());
  stats.AddData(ptxtBlocks, ptxtBlockSizes);
  cout << "Encryption time: " << (clock()-encryptionStart)/CLOCKS_PER_SEC << endl;
  
  vector<Ciphertext> encMean;
  Matrix<Ciphertext> encCov(Ciphertext(stats.GetPublicKey()));
  Ciphertext encN(stats.GetPublicKey());
  Ciphertext encN2(stats.GetPublicKey());
  
  double computationStart(clock());
  stats.ComputeCovariance(encCov, encMean, encN, encN2);
  
  cout << "Computation time: " << (clock()-computationStart)/CLOCKS_PER_SEC << endl;
  
  FHESISecKey secretKey = stats.GetSecretKey();
  
  Plaintext tmp(context);
  vector<ZZ_pX> msgs;
  ZZ_pX msg;

  double decStart(clock());
  cout << endl << "Computed values: " << endl;
  cout << "  Mean: ";
  for (unsigned i = 0; i < encMean.size(); i++) {
    secretKey.Decrypt(tmp, encMean[i]);
    tmp.DecodeSlots(msgs);
    cout << msgs[0] << ", ";
  }
  cout << endl;

  secretKey.Decrypt(tmp, encN);
  cout << "  N: " << tmp << endl << endl;
  
  cout << "  Covariance: " << endl;
  Matrix<ZZ_pX> decCov(encCov.NumRows(), encCov.NumCols());
  for (unsigned i = 0; i < encCov.NumRows(); i++) {
    for (unsigned j = i; j < encCov.NumCols(); j++) {
      secretKey.Decrypt(tmp, encCov(i,j));
      tmp.DecodeSlots(msgs);
      decCov(i,j) = decCov(j,i) = msgs[0];
    }
  }
  
  cout << decCov << endl;
  
  secretKey.Decrypt(tmp, encN2);
  cout << "  N^2: " << tmp << endl << endl;

  cout << "Decryption time: " << (clock()-decStart)/CLOCKS_PER_SEC << endl;
  cout << "Total time: " << (clock()-start)/CLOCKS_PER_SEC << endl;
  
  return 0;
}

int main(int argc, char *argv[]) {
  srand48(time(NULL));
  SetSeed(to_ZZ(time(NULL)));

  unsigned p, g, logQ = 0;
  char *datafile;
  
  if (argc >= 4) {
    datafile = argv[1];
    p = atoi(argv[2]);
    g = atoi(argv[3]);
  } else {
    cout << "usage: Test_Statistics_x datafile p generator" << endl;
    return 1;
  }
  
  unsigned blockSize = 1;
  unsigned val = (p-1)/2;
  while (val > 1) {
    blockSize <<= 1;
    val >>= 1;
  }
  
  Matrix<ZZ> data;
  unsigned dim;
  
  if (!LoadData(data, dim, datafile)) {
    cout << "Failed to load data file." << endl;
    return 1;
  }
  
  unsigned n = (p-1)/2-1;
  unsigned nBlocks = data.NumRows() / blockSize;
  if (data.NumRows() % blockSize != 0) {
    nBlocks++;
  }
  unsigned xi = max(nBlocks, dim);
  
  double lgQ = 6.5*log(n)+log(xi);
  logQ = (unsigned) ceil(lgQ / log(2) + 36.1);
  
  cout << "================================================" << endl
     << "Running statistics test using Brakerski system." << endl
     << "================================================" << endl;
  
  cout << "Parameters: " << endl
     << "  data file: " << datafile << endl
     << "  logQ: " << logQ << endl
     << "  p: " << p << endl
     << "  generator: " << g << endl
     << "  block size: " << blockSize << endl
     << "  num blocks: " << nBlocks << endl;
     
  
  cout << "Performing statistical analysis on " << data.NumRows()
         << " values (dimension " << dim << ") in " << (data.NumRows()+blockSize-1)/blockSize
         << " blocks, modulo prime " << p << endl << endl;
  
  vector<ZZ> mean;
  Matrix<ZZ> cov;
  
  FHEcontext context(p-1, logQ, p, g, 3);
  activeContext = &context;
  
  context.SetUpSIContext(xi);
  return RunStatisticsTest(data, p, context);
}
