#pragma once

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>

#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)
#define NULL 0

//defining prime polynomials in GF^x for your use
#define PRIM(x) PRIM_ ## x
#define PRIM_0 0
#define PRIM_1 0
#define PRIM_2 0
#define PRIM_3 0xb
#define PRIM_4 0x13
#define PRIM_5 0x25
#define PRIM_6 0x43
#define PRIM_7 0x83
#define PRIM_8 0x11d
#define PRIM_9 0x211
#define PRIM_10 0x409
#define PRIM_11 0x805
#define PRIM_12 0x1053
#define PRIM_13 0x201b
#define PRIM_14 0x402b
#define PRIM_15 0x8003
#define PRIM_16 0x1002d
#define PRIM_17 0x20009
#define PRIM_18 0x40027
#define PRIM_19 0x80027
#define PRIM_20 0x100009
#define PRIM_21 0x200005
#define PRIM_22 0x400003
#define PRIM_23 0x800021
#define PRIM_24 0x100001b
#define PRIM_25 0x2000009
#define PRIM_26 0x4000047
#define PRIM_27 0x8000027
#define PRIM_28 0x10000009
#define PRIM_29 0x20000005
#define PRIM_30 0x40000053
#define PRIM_31 0x80000009

// The underlying data type used can be changed depending on your space requirements
// just change this typedef
typedef unsigned long RS_WORD;

using namespace std;

void FindPrimePolys(ostream* out, int fieldPower, int limit);

class GaloisField
{
public:
	RS_WORD* powTable, *logTable;
	int fieldPower;
	RS_WORD characteristic, primitivePoly;
	
	GaloisField(int fieldPower);
	~GaloisField();
	RS_WORD multNoLUT(RS_WORD a, RS_WORD b);
	inline RS_WORD mult(RS_WORD a, RS_WORD b);
	inline RS_WORD div(RS_WORD a, RS_WORD b);
	inline RS_WORD pow(RS_WORD x, RS_WORD power);
	inline RS_WORD inv(RS_WORD x);
	inline RS_WORD sqrt(RS_WORD x);
};

void Init();


class Poly
{
public:
	int n;
	RS_WORD* coef;
	Poly();
	Poly(int n, RS_WORD* coef);
	~Poly();

	void init();
	void setCopy(int n, RS_WORD* coef);
	void setRef(int n, RS_WORD* coef);
	void print();
};

Poly* Poly_Create(int n, RS_WORD* coef);
void Poly_Free(Poly* poly);
void Poly_Add(Poly* out, Poly* a, Poly* b);
void Poly_Scale(Poly* out, Poly* in, RS_WORD scale, GaloisField* gf);
void Poly_Mult(Poly* out, Poly* a, Poly* b, GaloisField* gf);
void Poly_Div(Poly* result, Poly* quotient, Poly* remainder, Poly* a, Poly* b, GaloisField* gf);
RS_WORD Poly_Eval(Poly* poly, RS_WORD x, GaloisField* gf);
void Poly_ChienSearch(vector<unsigned int>* out, Poly* poly, int max, GaloisField* gf);
void Poly_Pad(Poly* poly, int left, int right);
void Poly_Trim(Poly* poly, int left, int right);
void Poly_Append(Poly* out, Poly* a, Poly* b);
void Poly_Reverse(Poly* out, Poly* in);

class ReedSolomon
{
public:
	GaloisField gf;
	ReedSolomon(int fieldPower);

	void createGenerator(Poly* out, int nsym);
	void encode(RS_WORD* out, RS_WORD* data, int k, int nsym);
	void calcSyndromes(Poly* out, Poly* msg, int nsym);
	bool checkSyndromes(Poly* synd);
	void findErrataLocator(Poly* out, vector<unsigned int>* errorPos);
	void findErrorEvaluator(Poly* out, Poly* synd, Poly* errLoc, int nsym);
	bool correctErrata(Poly* msg, Poly* synd, vector<unsigned int>* errPos);
	bool findErrorLocator(Poly* out, Poly* synd, int nsym, Poly* eraseLoc, int eraseCount);
	bool findErrors(vector<unsigned int>* out, Poly* errLoc, int n);
	void forneySyndromes(Poly* out, Poly* synd, vector<unsigned int>* pos, int n);
	bool decode(RS_WORD* wholeOut, RS_WORD* out, RS_WORD* data, int k, int nsym, vector<unsigned int>* erasePos, bool debug);
};

