#include "ReedSolomon.h"
using namespace std;

unsigned char logTable[256], powTable[512];
bool hasInit = false;


void Init()
{
	//init the tables 
	int val = 1;
	powTable[0] = val;
	for (int i = 1; i < 256; i++)
	{
		val <<= 1;
		if (val & 0x100)
		{
			val ^= PRIMITIVE_POLY;
		}
		powTable[i] = (unsigned char)val;
		logTable[(unsigned char)val] = i;
	}
	for (int i = 256; i < 512; i++)
	{
		powTable[i] = powTable[i - 255];
	}
	hasInit = true;
}

inline unsigned char GF_Mult(unsigned char a, unsigned char b)
{
	return (a == 0 || b == 0) ? 0 : powTable[logTable[a] + logTable[b]];
}
inline unsigned char GF_Div(unsigned char a, unsigned char b)
{
	return(b == 0) ? -1 : powTable[logTable[a] - logTable[b] + 255];
}

inline unsigned char GF_Pow(unsigned char x, unsigned char power)
{
	return powTable[(logTable[x] * power) % 255];
}

inline unsigned char GF_Inv(unsigned char x)
{
	return powTable[255 - logTable[x]];
}

inline unsigned char GF_Sqrt(unsigned char x)
{
	return logTable[x] % 2 ? powTable[(logTable[x] + 255) / 2] : powTable[logTable[x] / 2];
}

Poly::Poly()
{
	this->Init();
}

Poly::Poly(int n, unsigned char* data)
{
	this->Init();
	this->SetCopy(n, data);
}

Poly::~Poly()
{
	if (this->coef)
	{
		free(this->coef);
	}
}

void Poly::Init()
{
	this->n = 0;
	this->coef = nullptr;
}

void Poly::SetCopy(int n, unsigned char* coef)
{
	if (n > this->n)
	{
		if (this->coef)
		{
			free(this->coef);
		}
		this->coef = (unsigned char*)malloc(sizeof(char) * n);
	}
	this->n = n;
	if (coef)
	{
		memcpy(this->coef, coef, sizeof(char) * n);
	} else
	{
		memset(this->coef, 0, sizeof(char) * n);
	}
}

void Poly::SetRef(int n, unsigned char* coef)
{
	if (this->coef)
	{
		free(this->coef);
	}
	this->n = n;
	this->coef = coef;
}

Poly* Poly_Create(int n, unsigned char* coef)
{
	Poly* poly = (Poly*)malloc(sizeof(Poly));
	poly->Init();
	poly->SetCopy(n, coef);
	return poly;
}

void Poly_Free(Poly* poly)
{
	free(poly->coef);
	free(poly);
}

void Poly_Add(Poly* out, Poly* a, Poly* b)
{
	int n = MAX(a->n, b->n);
	unsigned char* temp = (unsigned char*)malloc(sizeof(char) * n);
	memset(temp, 0, sizeof(char) * n);
	for (int i = 0; i < a->n; i++)
	{
		temp[i + n - a->n] = a->coef[i];
	}
	for (int i = 0; i < b->n; i++)
	{
		temp[i + n - b->n] ^= b->coef[i];
	}
	out->SetRef(n, temp);
}

void Poly_Scale(Poly* out, Poly* in, unsigned char scale)
{
	if (out == in)
	{
		for (int i = 0; i < in->n; i++)
		{
			in->coef[i] = GF_Mult(in->coef[i], scale);
		}
	} else
	{
		unsigned char* temp = (unsigned char*)malloc(sizeof(char) * in->n);
		for (int i = 0; i < in->n; i++)
		{
			temp[i] = GF_Mult(in->coef[i], scale);
		}
		out->SetRef(in->n, temp);
	}
}

void Poly_Mult(Poly* out, Poly* a, Poly* b)
{
	int n = a->n + b->n - 1;
	unsigned char* temp = (unsigned char*)malloc(sizeof(char) * n);
	memset(temp, 0, sizeof(char) * n);
	for (int i = 0; i < a->n; i++)
	{
		for (int j = 0; j < b->n; j++)
		{
			temp[i + j] ^= GF_Mult(a->coef[i], b->coef[j]);
		}
	}
	out->SetRef(n, temp);
}

void Poly_Div(Poly* result, Poly* quotient, Poly* remainder, Poly* a, Poly* b)
{
	unsigned char* temp = (unsigned char*)malloc(sizeof(char) * a->n);
	unsigned char normalizer = b->coef[0];
	memcpy(temp, a->coef, sizeof(char) * a->n);
	for (int i = 0; i < a->n - b->n + 1; i++)
	{
		temp[i] = GF_Div(temp[i], normalizer);
		unsigned char coef = temp[i];
		if (coef != 0)
		{
			for (int j = 1; j < b->n; j++)
			{
				if (b->coef[j] != 0)
				{
					temp[i + j] ^= GF_Mult(b->coef[j], coef);
				}
			}
		}
	}
	if (result)
	{
		result->SetCopy(a->n, temp);
	}
	int separator = a->n - b->n + 1;
	if (quotient)
	{
		quotient->SetCopy(separator, temp);
	}
	if (remainder)
	{
		remainder->SetCopy(b->n - 1, temp + separator);
	}
	free(temp);
}

unsigned char Poly_Eval(Poly* poly, unsigned char x)
{
	unsigned char y = poly->coef[0];
	for (int i = 1; i < poly->n; i++)
	{
		y = GF_Mult(y, x) ^ poly->coef[i];
	}
	return y;
}

void Poly_ChienSearch(vector<unsigned char>* out, Poly* poly, int max)
{
	//this seems unnecessary because all multiplications are performed via lookup table anyway
	unsigned char* temp = (unsigned char*)malloc(sizeof(char) * poly->n);
	memcpy(temp, poly->coef, poly->n);
	for (int i = 0; i < max; i++)
	{
		unsigned char sum = 0;
		for (int j = 0; j < poly->n; j++)
		{
			sum ^= temp[j];
			temp[j] = GF_Mult(temp[j], powTable[poly->n - j - 1]);
		}
		if (!sum)
		{
			out->push_back(i);
		}
	}
	free(temp);
}

void Poly_Pad(Poly* poly, int left, int right)
{
	int n = poly->n + left + right;
	unsigned char* temp = (unsigned char*)malloc(sizeof(char) * n);
	memset(temp, 0, sizeof(char) * left);
	memcpy(temp + left, poly->coef, poly->n);
	memset(temp + left + poly->n, 0, sizeof(char) * right);
	poly->SetRef(n, temp);
}

void Poly_Trim(Poly* poly, int left, int right)
{
	int n = poly->n - left - right;
	unsigned char* temp = (unsigned char*)malloc(sizeof(char) * n);
	memcpy(temp, poly->coef + left, n);
	poly->SetRef(n, temp);
}

void Poly_Append(Poly* out, Poly* a, Poly* b)
{
	int n = a->n + b->n;
	unsigned char* temp = (unsigned char*)malloc(sizeof(char) * n);
	memcpy(temp, a->coef, sizeof(char) * a->n);
	memcpy(temp + a->n, b->coef, sizeof(char) * b->n);
	out->SetRef(n, temp);
}

void Poly_Reverse(Poly* out, Poly* in)
{
	unsigned char* temp = (unsigned char*)malloc(sizeof(char) * in->n);
	for (int i = 0; i < in->n; i++)
	{
		temp[i] = in->coef[in->n - i - 1];
	}
	out->SetRef(in->n, temp);
}

void RS_CreateGenerator(Poly* out, int nsym)
{
	out->SetCopy(1, nullptr);
	out->coef[0] = 1;
	Poly factor(2, nullptr);
	factor.coef[0] = 1;
	for (int i = 0; i < nsym; i++)
	{
		factor.coef[1] = powTable[i];
		Poly_Mult(out, out, &factor);
	}
}

void RS_Encode(unsigned char* out, unsigned char* data, int k, int nsym)
{
	if (!hasInit)
	{
		Init();
	}
	Poly msg(k, data);
	Poly generator, remainder;
	RS_CreateGenerator(&generator, nsym);
	Poly_Pad(&msg, 0, nsym);
	Poly_Div(nullptr, nullptr, &remainder, &msg, &generator);
	memcpy(out, data, sizeof(char) * k);
	memcpy(out + k, remainder.coef, sizeof(char) * remainder.n);
}

void RS_CalcSyndromes(Poly* out, Poly* msg, int nsym)
{
	unsigned char* temp = (unsigned char*)malloc(sizeof(char) * (nsym + 1));
	for (int i = 0; i < nsym; i++)
	{
		temp[nsym - i - 1] = Poly_Eval(msg, powTable[i]);
	}
	temp[nsym] = 0; //pad
	out->SetRef(nsym + 1, temp);
}

bool RS_CheckSyndromes(Poly* synd)
{
	for (int i = 0; i < synd->n; i++)
	{
		if (synd->coef[i])
		{
			return false;
		}
	}
	return true;
}

void RS_FindErrataLocator(Poly* out, vector<unsigned char>* errPos)
{
	out->SetCopy(1, nullptr);
	out->coef[0] = 1;
	Poly factor(2, nullptr);
	factor.coef[1] = 1;
	for (int i : *errPos)
	{
		factor.coef[0] = powTable[i];
		Poly_Mult(out, out, &factor);
	}
}

void RS_FindErrorEvaluator(Poly* out, Poly* synd, Poly* errLoc, int nsym)
{
	Poly_Mult(out, synd, errLoc); //synd lul
	Poly_Trim(out, out->n - nsym, 0);
}

bool RS_CorrectErrata(Poly* msg, Poly* synd, vector<unsigned char>* errPos)
{
	vector<unsigned char> coefPos(0);
	for (int i : *errPos)
	{
		coefPos.push_back(msg->n - 1 - i);
	}
	Poly errLoc, errEval;
	RS_FindErrataLocator(&errLoc, &coefPos);
	RS_FindErrorEvaluator(&errEval, synd, &errLoc, errLoc.n); //TODO determine if correct
	//Poly_Reverse(errEval, errEval); //reverse it for later use
	vector<unsigned char> x(coefPos.size());
	for (int i = 0; i < x.size(); i++)
	{
		x[i] = powTable[coefPos[i]]; //TODO determine if correct
	}
	Poly e(msg->n, nullptr);
	for (int i = 0; i < x.size(); i++)
	{
		unsigned char xi = powTable[255 - coefPos[i]]; //TODO determine if equivalent to GF_Inv(x[i])
		unsigned char errLocPrime = 1;
		for (int j = 0; j < x.size(); j++)
		{
			if (j != i)
			{
				errLocPrime = GF_Mult(errLocPrime, 1 ^ GF_Mult(xi, x[j]));
			}
		}
		if (errLocPrime == 0)
		{
			return false;
		}
		unsigned char y = GF_Mult(x[i], Poly_Eval(&errEval, xi)); //errEval is still reversed from earlier
		//TODO determine if equivalent to GF_Mult(GF_Pow(xi, 1), y)
		e.coef[errPos->at(i)] = GF_Div(y, errLocPrime); //magnitude
	}
	Poly_Add(msg, msg, &e);
	return true;
}

bool RS_FindErrorLocator(Poly* out, Poly* synd, int nsym, Poly* eraseLoc, int eraseCount)
{
	//this spits out a polynomial in reverse order but i dont know why
	unsigned char init = 1;
	Poly errLoc(1, &init);
	Poly oldLoc(1, &init);
	Poly temp;
	if (eraseLoc)
	{
		errLoc.SetCopy(eraseLoc->n, eraseLoc->coef);
		oldLoc.SetCopy(eraseLoc->n, eraseLoc->coef);
	}
	int syndShift = 0; //TODO optimize
	for (int i = nsym - eraseCount - 1; i >= 0; i--)
	{
		int K = i + syndShift + eraseCount;
		unsigned char delta = synd->coef[K];
		for (int j = 1; j < errLoc.n; j++)
		{
			delta ^= GF_Mult(errLoc.coef[errLoc.n - j - 1], synd->coef[K + j]);
		}
		Poly_Pad(&oldLoc, 0, 1); //TODO optimize
		if (delta != 0)
		{
			if (oldLoc.n > errLoc.n)
			{ //TODO isn't this always the case?
				Poly_Scale(&temp, &oldLoc, delta);
				Poly_Scale(&oldLoc, &errLoc, GF_Inv(delta));
				errLoc.SetCopy(temp.n, temp.coef);
			}
			Poly_Scale(&temp, &oldLoc, delta);
			Poly_Add(&errLoc, &errLoc, &temp);
		}
	}
	int leading = 0;
	for (; errLoc.coef[leading] == 0; leading++);
	Poly_Trim(&errLoc, leading, 0);
	int errs = errLoc.n - 1;
	out->SetCopy(errLoc.n, errLoc.coef);
	if (errs * 2 - eraseCount > nsym)
	{
		return false;
	}
	return true;
}

void RS_FindErrors(vector<unsigned char>* out, Poly* errLoc, int n)
{
	int errs = errLoc->n - 1;
	Poly revErrLoc;
	Poly_Reverse(&revErrLoc, errLoc);
	if (errLoc->n == 2)
	{ //linear equation
		out->push_back(logTable[GF_Div(errLoc->coef[0], errLoc->coef[1])]);
	} else
	{
		Poly_ChienSearch(out, &revErrLoc, n);
	}
	//map to string pos
	for (int i = 0; i < out->size(); i++)
	{
		(*out)[i] = n - out->at(i) - 1;
	}
}

void RS_ForneySyndromes(Poly* out, Poly* synd, vector<unsigned char>* pos, int n)
{
	Poly fsynd(synd->n - 1, synd->coef);
	if (pos)
	{
		for (unsigned char i : *pos)
		{
			unsigned char rev = n - i - 1;
			unsigned char x = powTable[rev];
			for (int j = fsynd.n - 2; j >= 0; j--)
			{
				fsynd.coef[j + 1] = GF_Mult(fsynd.coef[j + 1], x) ^ fsynd.coef[j];
			}
		}
	}
	//fsynd.coef[fsynd.n - 1] = 0;
	out->SetCopy(fsynd.n, fsynd.coef);
}

bool RS_Decode(unsigned char* wholeOut, unsigned char* out, unsigned char* data, int k, int nsym, vector<unsigned char>* erasePos, bool debug)
{
	if (!hasInit)
	{
		Init();
	}
	Poly synd;
	Poly msg(k + nsym, data);
	if (erasePos)
	{
		if (erasePos->size() > nsym)
		{
			if (debug) cout << "too many erasures to be corrected" << endl;
			return false;
		} else
		{
			for (unsigned char i : *erasePos)
			{
				msg.coef[i] = 0;
			}
		}
	}
	RS_CalcSyndromes(&synd, &msg, nsym);
	if (!RS_CheckSyndromes(&synd))
	{
		if(debug) cout << "errors detected, locating" << endl;
		Poly fsynd, errLoc;
		RS_ForneySyndromes(&fsynd, &synd, erasePos, k + nsym);
		RS_FindErrorLocator(&errLoc, &fsynd, nsym, nullptr, erasePos ? erasePos->size() : 0);
		vector<unsigned char> pos(0);
		RS_FindErrors(&pos, &errLoc, k + nsym);
		if (!pos.size())
		{
			if(debug) cout << "errors unable to be located" << endl;
			return false;
		}
		if (debug)
		{
			cout << "errors detected at ";
			for_each(pos.begin(), pos.end(), [](unsigned char e) {cout << (int)e << ", "; });
			cout << "correcting" << endl;
		}
		if (erasePos)
		{
			pos.insert(pos.begin(), erasePos->begin(), erasePos->end());
		}
		bool success = RS_CorrectErrata(&msg, &synd, &pos);
		if (!success)
		{
			if(debug) cout << "decode failure" << endl;
			return false;
		}
		if(debug) cout << "errors corrected!" << endl;
	}
	if (wholeOut)
	{
		memcpy(wholeOut, msg.coef, sizeof(char) * (k + nsym));
	}
	if (out)
	{
		memcpy(out, msg.coef, sizeof(char) * k);
	}
	return true;
}