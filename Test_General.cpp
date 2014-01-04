#include "FHE-SI.h"
#include "Ciphertext.h"
#include "Plaintext.h"

void randomizePlaintext(Plaintext &ptxt, unsigned deg, unsigned p) {
	ZZ_pX poly;
	poly.rep.SetLength(deg);
	for(unsigned i = 0; i < deg; i++) {
		poly.rep[i] = RandomBnd(p);
	}

	poly.normalize();
	ptxt.Init(poly);
}

int main(int argc, char *argv[]) {
	long seed = time(NULL);

	srand48(seed);
	SetSeed(to_ZZ(seed));

	unsigned p = 2027;
	unsigned g = 3;
	unsigned logQ = 120;

	FHEcontext context(p-1, logQ, p, g);
	activeContext = &context;
	context.SetUpSIContext();

	FHESISecKey secretKey(context);
	const FHESIPubKey &publicKey(secretKey);
	KeySwitchSI keySwitch(secretKey);

	long phim = context.zMstar.phiM();
	long numSlots = context.GetPlaintextSpace().GetTotalSlots();

	long rotAmt = rand() % numSlots;
	long rotDeg = 1;
	for (int i = 0; i < rotAmt; i++) {
		rotDeg *= context.Generator();
		rotDeg %= context.zMstar.M();
	}
	KeySwitchSI automorphKeySwitch(secretKey, rotDeg);

	Plaintext p0, p1, p2, p3;
	randomizePlaintext(p0, phim, p);
	randomizePlaintext(p1, phim, p);
	randomizePlaintext(p2, phim, p);
	randomizePlaintext(p3, phim, p);

	Plaintext const1, const2;
	randomizePlaintext(const1, phim, p);
	randomizePlaintext(const2, phim, p);

	Ciphertext c0(publicKey);
	Ciphertext c1(publicKey);
	Ciphertext c2(publicKey);
	Ciphertext c3(publicKey);

	publicKey.Encrypt(c0, p0);
	publicKey.Encrypt(c1, p1);
	publicKey.Encrypt(c2, p2);
	publicKey.Encrypt(c3, p3);

	p1 *= p2;
	p0 += const1;
	p2 *= const2;
	p3 >>= rotAmt;
	p1 *= -1;
	p3 *= p2;
	p0 -= p3;

	c1 *= c2;
	keySwitch.ApplyKeySwitch(c1);
	c0 += const1.message;

	c2 *= const2.message;

	c3 >>= rotDeg;
	automorphKeySwitch.ApplyKeySwitch(c3);

	c1 *= -1;
	c3 *= c2;
	keySwitch.ApplyKeySwitch(c3);

	Ciphertext tmp(c3);
	tmp *= -1;
	c0 += tmp;

	Plaintext pp0, pp1, pp2, pp3;
	secretKey.Decrypt(pp0, c0);
	secretKey.Decrypt(pp1, c1);
	secretKey.Decrypt(pp2, c2);
	secretKey.Decrypt(pp3, c3);

    if (!(pp0 == p0)) cerr << "oops 0" << endl;
    if (!(pp1 == p1)) cerr << "oops 1" << endl;
    if (!(pp2 == p2)) cerr << "oops 2" << endl;
    if (!(pp3 == p3)) cerr << "oops 3" << endl;
    cout << "All tests finished." << endl;
}