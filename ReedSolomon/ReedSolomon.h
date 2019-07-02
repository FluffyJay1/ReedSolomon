#pragma once

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>

#define PRIMITIVE_POLY 0x11D
#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)
#define NULL 0

using namespace std;

void Init();
inline unsigned char GF_Mult(unsigned char a, unsigned char b);
inline unsigned char GF_Div(unsigned char a, unsigned char b);
inline unsigned char GF_Pow(unsigned char x, unsigned char power);
inline unsigned char GF_Inv(unsigned char x);
inline unsigned char GF_Sqrt(unsigned char x);

class Poly
{
public:
	int n;
	unsigned char* coef;
	Poly();
	Poly(int n, unsigned char* coef);
	~Poly();

	void Init();
	void SetCopy(int n, unsigned char* coef);
	void SetRef(int n, unsigned char* coef);
};

Poly* Poly_Create(int n, unsigned char* coef);
void Poly_Free(Poly* poly);
void Poly_Add(Poly* out, Poly* a, Poly* b);
void Poly_Scale(Poly* out, Poly* in, unsigned char scale);
void Poly_Div(Poly* result, Poly* quotient, Poly* remainder, Poly* a, Poly* b);
unsigned char Poly_Eval(Poly* poly, unsigned char x);
void Poly_ChienSearch(vector<unsigned char>* out, Poly* poly, int max);
void Poly_Pad(Poly* poly, int left, int right);
void Poly_Trim(Poly* poly, int left, int right);
void Poly_Append(Poly* out, Poly* a, Poly* b);
void Poly_Reverse(Poly* out, Poly* in);

void RS_CreateGenerator(Poly* out, int nsym);
void RS_Encode(unsigned char* out, unsigned char* data, int k, int nsym);
void RS_CalcSyndromes(Poly* out, Poly* msg, int nsym);
void RS_FindErrataLocator(Poly* out, vector<unsigned char>* errorPos);
void RS_FindErrorEvaluator(Poly* out, Poly* synd, Poly* errLoc, int nsym);
bool RS_CorrectErrata(Poly* msg, Poly* synd, vector<unsigned char>* errPos);
bool RS_FindErrorLocator(Poly* out, Poly* synd, int nsym, Poly* eraseLoc, int eraseCount);
void RS_FindErrors(vector<unsigned char>* out, Poly* errLoc, int n);
void RS_ForneySyndromes(Poly* out, Poly* synd, vector<unsigned char>* pos, int n);
bool RS_Decode(unsigned char* wholeOut, unsigned char* out, unsigned char* data, int k, int nsym, vector<unsigned char>* erasePos, bool debug);