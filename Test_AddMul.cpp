#define N_TESTS 5000

#include "Plaintext.h"
#include "DoubleCRT.h"
#include "Ciphertext.h"
#include "NTL/ZZ_pX.h"

#include <time.h>
#include "FHE-SI.h"

bool runTest(bool disp, long long seed, unsigned p, FHEcontext &context) {
  ZZ seedZZ;
  seedZZ = seed;

  srand48(seed);
  SetSeed(seedZZ);

  FHESISecKey secretKey(context);
  const FHESIPubKey &publicKey(secretKey);
  
  long phim = context.zMstar.phiM();
  
  ZZ_pX ptxt1Poly, ptxt2Poly, sum, sumMult, prod, prod2, sumQuad;
  Plaintext resSum, resSumMult, resProd, resProdSwitch, resProd2, resSumQuad;

  ptxt1Poly.rep.SetLength(phim);
  ptxt2Poly.rep.SetLength(phim);
  for (long i=0; i < phim; i++) {
    ptxt1Poly.rep[i] = RandomBnd(p);
    ptxt2Poly.rep[i] = RandomBnd(p);
  }

  ptxt1Poly.normalize();
  ptxt2Poly.normalize();
  
  sum = ptxt1Poly + ptxt2Poly;
  sumMult = ptxt2Poly * 7;
  prod = ptxt1Poly * ptxt2Poly;
  prod2 = prod * prod;
  sumQuad = prod2 * prod2 * 9;

  rem(prod, prod, to_ZZ_pX(context.zMstar.PhimX()));
  rem(prod2, prod2, to_ZZ_pX(context.zMstar.PhimX()));
  rem(sumQuad, sumQuad, to_ZZ_pX(context.zMstar.PhimX()));
  
  Plaintext ptxt1(context, ptxt1Poly), ptxt2(context, ptxt2Poly);
  Ciphertext ctxt1(publicKey), ctxt2(publicKey);
  publicKey.Encrypt(ctxt1, ptxt1);
  publicKey.Encrypt(ctxt2, ptxt2);

  Ciphertext cSum = ctxt1;
  cSum += ctxt2;

  Ciphertext cSumMult = ctxt2;
  for (int i = 1; i < 7; i++) {
    cSumMult += ctxt2;
  }

  Ciphertext cProd = ctxt1;
  cProd *= ctxt2;

  secretKey.Decrypt(resSum, cSum);
  secretKey.Decrypt(resSumMult, cSumMult);
  
  KeySwitchSI keySwitch(secretKey);
  keySwitch.ApplyKeySwitch(cProd);
  secretKey.Decrypt(resProd, cProd);

  cProd *= cProd;
  Ciphertext tmp = cProd;
  Ciphertext cSumQuad = cProd;
  
  keySwitch.ApplyKeySwitch(cProd);
  secretKey.Decrypt(resProd2, cProd);

  for (int i = 0; i < 8; i++) {
    cSumQuad += tmp;
  }
  keySwitch.ApplyKeySwitch(cSumQuad);
  cSumQuad *= cProd;
  keySwitch.ApplyKeySwitch(cSumQuad);
  secretKey.Decrypt(resSumQuad, cSumQuad);
  
  bool success = ((resSum.message == sum) && (resSumMult.message == sumMult) &&
                  (resProd.message == prod) && (resProd2.message == prod2) &&
                  (resSumQuad == sumQuad));
  
  if (disp || !success) {
    cout << "Seed: " << seed << endl << endl;
    
    if (resSum.message != sum) {
      cout << "Add failed." << endl;
    }
    if (resSumMult.message != sumMult) {
      cout << "Adding multiple times failed." << endl;
    }
    if (resProd.message != prod) {
      cout << "Multiply failed." << endl;
    }
    if (resProd2.message != prod2) {
      cout << "Squaring failed." << endl;
    }
    if (resSumQuad.message != sumQuad) {
      cout << "Sum and quad failed." << endl;
    }
  }
  
  if (disp || !success) {
    cout << "Test " << (success ? "SUCCEEDED" : "FAILED") << endl;
  }

  return success;
}

int main(int argc, char *argv[]) {
  unsigned p, g, logQ;

  if (argc >= 4) {
    logQ = atoi(argv[1]);
    p = atoi(argv[2]);
    g = atoi(argv[3]);
  } else {
    cout << "usage: Test_AddMul_x logQ p generator [seed]" << endl;
    return 1;
  }
  
  cout << "==================================================" << endl
       << "Running add/multiply tests using Brakerski system." << endl
       << "==================================================" << endl;

  FHEcontext context(p-1, logQ, p, g, 3);
  activeContext = &context;
  
  context.SetUpSIContext();
  
  cout << "Finished setting up context." << endl;

  if (argc >= 5) {
    long long seed = atoi(argv[4]);
    bool res = runTest(true, seed, p, context);
    
    if (!res) {
      cout << "Failed test with seed " << seed << endl;
      return 1;
    }
    return 0;
  }

  int testsFailed = 0;

  srand48(time(NULL));
  long long start = lrand48();
  for (int iter = 0; iter < N_TESTS; iter++) {
    if (!runTest(false, start+iter, p, context)) {
      testsFailed++;
    }

    if (iter % 100 == 0) {
      cout << "." << flush;
    }
  }
  cout << endl;

  if (testsFailed == 0) {
    cout << "All tests SUCCEEDED!" << endl;
  } else{
    cout << testsFailed << " of " << N_TESTS << " failed." << endl;
  }

  return testsFailed;
}
