#include "ReedSolomon.h"
#include <string.h>
using namespace std;

unsigned int primes[] = {
	PRIM(0), PRIM(1), PRIM(2), PRIM(3), 
	PRIM(4), PRIM(5), PRIM(6), PRIM(7), 
	PRIM(8), PRIM(9), PRIM(10), PRIM(11), 
	PRIM(12), PRIM(13), PRIM(14), PRIM(15),
	PRIM(16), PRIM(17), PRIM(18), PRIM(19), 
	PRIM(20), PRIM(21), PRIM(22), PRIM(23),
	PRIM(24), PRIM(25), PRIM(26), PRIM(27),
	PRIM(28), PRIM(29), PRIM(30), PRIM(31)
};

void FindPrimePolys(ostream* out, int fieldPower, int limit)
{
	GaloisField gf(fieldPower);
	int primesFound = 0;
	RS_WORD fieldCharacteristic = ((RS_WORD)1 << fieldPower) - 1, fieldCharacteristicNext = ((RS_WORD)1 << (fieldPower + 1)) - 1;
	for (RS_WORD i = fieldCharacteristic + 2; (fieldCharacteristicNext == 0 ? i > 0 : i < fieldCharacteristicNext); i += 2)
	{
		unsigned int x = 2; //skip first iteration
		bool conflict = false;
		for (unsigned int j = 1; j < fieldCharacteristic; j++)
		{
			x <<= 1;
			if (x > fieldCharacteristic)
			{
				x ^= i;
			}
			if (x == 2) //cyclical, 2 will always repeat first
			{
				conflict = true;
				break;
			}
		}
		if (!conflict)
		{
			*out << hex << i << endl;
			primesFound++;
			if (primesFound >= limit)
			{
				return;
			}
		}
	}
}

GaloisField::GaloisField(int fieldPower)
{
	this->characteristic = ((RS_WORD)1 << fieldPower) - 1;
	this->fieldPower = fieldPower;
	this->primitivePoly = primes[fieldPower];
	//init the tables 
	unsigned int val = 1;
	this->powTable = (RS_WORD*)malloc(sizeof(RS_WORD) * this->characteristic * 2);
	this->logTable = (RS_WORD*)malloc(sizeof(RS_WORD) * (this->characteristic + 1));
	powTable[0] = val;
	logTable[0] = 0;
	logTable[1] = 0;
	for (RS_WORD i = 1; i < this->characteristic; i++)
	{
		val <<= 1;
		if (val > this->characteristic)
		{
			val ^= this->primitivePoly;
		}
		powTable[i] = (RS_WORD)val;
		logTable[(RS_WORD)val] = i;
	}
	for (RS_WORD i = this->characteristic; i < this->characteristic * 2; i++)
	{
		powTable[i] = powTable[i - this->characteristic];
	}
	/*
	for (unsigned int i = 0; i < this->characteristic * 2; i++)
	{
		cout << "2^" << i << "=" << powTable[i] << endl;
	}
	for (unsigned int i = 0; i < this->characteristic + 1; i++)
	{
		cout << "log" << i << "=" << logTable[i] << endl;
	}
	*/
}

GaloisField::~GaloisField()
{
	free(this->powTable);
	free(this->logTable);
}

RS_WORD GaloisField::multNoLUT(RS_WORD a, RS_WORD b)
{
	RS_WORD ret = 0;
	while (b > 0)
	{
		if (b & 1) //if odd
		{
			ret ^= a;
		}
		b >>= 1;
		a <<= 1;
		if (a > this->characteristic)
		{
			a ^= this->primitivePoly;
		}
	}
	return ret;
}

inline RS_WORD GaloisField::mult(RS_WORD a, RS_WORD b)
{
	return (a == 0 || b == 0) ? 0 : this->powTable[this->logTable[a] + this->logTable[b]];
}
inline RS_WORD GaloisField::div(RS_WORD a, RS_WORD b)
{
	return a == 0 ? 0 : (b == 0) ? -1 : this->powTable[this->logTable[a] - this->logTable[b] + this->characteristic];
}

inline RS_WORD GaloisField::pow(RS_WORD x, RS_WORD power)
{
	return this->powTable[(this->logTable[x] * power) % this->characteristic];
}

inline RS_WORD GaloisField::inv(RS_WORD x)
{
	return this->powTable[this->characteristic - this->logTable[x]];
}

inline RS_WORD GaloisField::sqrt(RS_WORD x)
{
	return logTable[x] % 2 ? this->powTable[(logTable[x] + this->characteristic) / 2] : this->powTable[logTable[x] / 2];
}

Poly::Poly()
{
	this->init();
}

Poly::Poly(int n, RS_WORD* data)
{
	this->init();
	this->setCopy(n, data);
}

Poly::~Poly()
{
	if (this->coef)
	{
		free(this->coef);
	}
}

void Poly::init()
{
	this->n = 0;
	this->coef = nullptr;
}

void Poly::setCopy(int n, RS_WORD* coef)
{
	if (n > this->n)
	{
		if (this->coef)
		{
			free(this->coef);
		}
		this->coef = (RS_WORD*)malloc(sizeof(RS_WORD) * n);
	}
	this->n = n;
	if (coef)
	{
		memcpy(this->coef, coef, sizeof(RS_WORD) * n);
	} else
	{
		memset(this->coef, 0, sizeof(RS_WORD) * n);
	}
}

void Poly::setRef(int n, RS_WORD* coef)
{
	if (this->coef)
	{
		free(this->coef);
	}
	this->n = n;
	this->coef = coef;
}

void Poly::print()
{
	cout << "Poly(n=" << this->n << ")";
	if (this->n > 0)
	{
		cout << "[" << setw(3) << hex << (RS_WORD)this->coef[0];
	}
	for (int i = 1; i < this->n; i++)
	{
		cout << ", " << setw(3) << hex << (RS_WORD)this->coef[i];
	}
	if (this->n > 0)
	{
		cout << "]" << dec << endl;
	}
}

Poly* Poly_Create(int n, RS_WORD* coef)
{
	Poly* poly = (Poly*)malloc(sizeof(Poly));
	poly->init();
	poly->setCopy(n, coef);
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
	RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD) * n);
	memset(temp, 0, sizeof(RS_WORD) * n);
	for (int i = 0; i < a->n; i++)
	{
		temp[i + n - a->n] = a->coef[i];
	}
	for (int i = 0; i < b->n; i++)
	{
		temp[i + n - b->n] ^= b->coef[i];
	}
	out->setRef(n, temp);
}

void Poly_Scale(Poly* out, Poly* in, RS_WORD scale, GaloisField* gf)
{
	if (out == in)
	{
		for (int i = 0; i < in->n; i++)
		{
			in->coef[i] = gf->mult(in->coef[i], scale);
		}
	} else
	{
		RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD) * in->n);
		for (int i = 0; i < in->n; i++)
		{
			temp[i] = gf->mult(in->coef[i], scale);
		}
		out->setRef(in->n, temp);
	}
}

void Poly_Mult(Poly* out, Poly* a, Poly* b, GaloisField* gf)
{
	int n = a->n + b->n - 1;
	RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD) * n);
	memset(temp, 0, sizeof(RS_WORD) * n);
	for (int i = 0; i < a->n; i++)
	{
		for (int j = 0; j < b->n; j++)
		{
			temp[i + j] ^= gf->mult(a->coef[i], b->coef[j]);
		}
	}
	out->setRef(n, temp);
}

void Poly_Div(Poly* result, Poly* quotient, Poly* remainder, Poly* a, Poly* b, GaloisField* gf)
{
	RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD)* a->n);
	RS_WORD normalizer = b->coef[0];
	memcpy(temp, a->coef, sizeof(RS_WORD) * a->n);
	for (int i = 0; i < a->n - b->n + 1; i++)
	{
		temp[i] = gf->div(temp[i], normalizer);
		RS_WORD coef = temp[i];
		if (coef != 0)
		{
			for (int j = 1; j < b->n; j++)
			{
				if (b->coef[j] != 0)
				{
					temp[i + j] ^= gf->mult(b->coef[j], coef);
				}
			}
		}
	}
	if (result)
	{
		result->setCopy(a->n, temp);
	}
	int separator = a->n - b->n + 1;
	if (quotient)
	{
		quotient->setCopy(separator, temp);
	}
	if (remainder)
	{
		remainder->setCopy(b->n - 1, temp + separator);
	}
	free(temp);
}

RS_WORD Poly_Eval(Poly* poly, RS_WORD x, GaloisField* gf)
{
	RS_WORD y = poly->coef[0];
	for (int i = 1; i < poly->n; i++)
	{
		y = gf->mult(y, x) ^ poly->coef[i];
	}
	return y;
}

void Poly_ChienSearch(vector<unsigned int>* out, Poly* poly, int max, GaloisField* gf)
{
	//this seems unnecessary because all multiplications are performed via lookup table anyway
	RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD)* poly->n);
	memcpy(temp, poly->coef, sizeof(RS_WORD) * poly->n);
	for (int i = 0; i < max; i++)
	{
		RS_WORD sum = 0;
		for (int j = 0; j < poly->n; j++)
		{
			sum ^= temp[j];
			temp[j] = gf->mult(temp[j], gf->powTable[poly->n - j - 1]);
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
	RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD)* n);
	memset(temp, 0, sizeof(RS_WORD) * left);
	memcpy(temp + left, poly->coef, sizeof(RS_WORD) * poly->n);
	memset(temp + (left + poly->n), 0, sizeof(RS_WORD) * right);
	poly->setRef(n, temp);
}

void Poly_Trim(Poly* poly, int left, int right)
{
	int n = poly->n - left - right;
	RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD)* n);
	memcpy(temp, poly->coef + left, sizeof(RS_WORD) * n);
	poly->setRef(n, temp);
}

void Poly_Append(Poly* out, Poly* a, Poly* b)
{
	int n = a->n + b->n;
	RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD)* n);
	memcpy(temp, a->coef, sizeof(RS_WORD)* a->n);
	memcpy(temp + a->n, b->coef, sizeof(RS_WORD)* b->n);
	out->setRef(n, temp);
}

void Poly_Reverse(Poly* out, Poly* in)
{
	RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD)* in->n);
	for (int i = 0; i < in->n; i++)
	{
		temp[i] = in->coef[in->n - i - 1];
	}
	out->setRef(in->n, temp);
}

ReedSolomon::ReedSolomon(int fieldPower) : gf(fieldPower)
{
	
}

void ReedSolomon::createGenerator(Poly* out, int nsym)
{
	out->setCopy(1, nullptr);
	out->coef[0] = 1;
	Poly factor(2, nullptr);
	factor.coef[0] = 1;
	for (int i = 0; i < nsym; i++)
	{
		factor.coef[1] = this->gf.powTable[i];
		Poly_Mult(out, out, &factor, &this->gf);
	}
}

void ReedSolomon::encode(RS_WORD* out, RS_WORD* data, int k, int nsym)
{
	Poly msg(k, data);
	Poly generator, remainder;
	this->createGenerator(&generator, nsym);
	Poly_Pad(&msg, 0, nsym);
	Poly_Div(nullptr, nullptr, &remainder, &msg, &generator, &this->gf);
	memcpy(out, data, sizeof(RS_WORD) * k);
	memcpy(out + k, remainder.coef, sizeof(RS_WORD) * remainder.n);
}

void ReedSolomon::calcSyndromes(Poly* out, Poly* msg, int nsym)
{
	RS_WORD* temp = (RS_WORD*)malloc(sizeof(RS_WORD) * (nsym + 1));
	for (int i = 0; i < nsym; i++)
	{
		temp[nsym - i - 1] = Poly_Eval(msg, this->gf.powTable[i], &this->gf);
	}
	temp[nsym] = 0; //pad
	out->setRef(nsym + 1, temp);
}

bool ReedSolomon::checkSyndromes(Poly* synd)
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

void ReedSolomon::findErrataLocator(Poly* out, vector<unsigned int>* errPos)
{
	out->setCopy(1, nullptr);
	out->coef[0] = 1;
	Poly factor(2, nullptr);
	factor.coef[1] = 1;
	for (unsigned int i : *errPos)
	{
		factor.coef[0] = this->gf.powTable[i];
		Poly_Mult(out, out, &factor, &this->gf);
	}
}

void ReedSolomon::findErrorEvaluator(Poly* out, Poly* synd, Poly* errLoc, int nsym)
{
	Poly_Mult(out, synd, errLoc, &this->gf); //synd lul
	Poly_Trim(out, out->n - nsym, 0);
}

bool ReedSolomon::correctErrata(Poly* msg, Poly* synd, vector<unsigned int>* errPos)
{
	vector<unsigned int> coefPos(0);
	for (unsigned int i : *errPos)
	{
		coefPos.push_back(msg->n - 1 - i);
	}
	Poly errLoc, errEval;
	this->findErrataLocator(&errLoc, &coefPos);
	this->findErrorEvaluator(&errEval, synd, &errLoc, errLoc.n); //TODO determine if correct
	//Poly_Reverse(errEval, errEval); //reverse it for later use
	vector<RS_WORD> x(coefPos.size());
	for (int i = 0; i < x.size(); i++)
	{
		x[i] = this->gf.powTable[coefPos[i]]; //TODO determine if correct
	}
	Poly e(msg->n, nullptr);
	for (int i = 0; i < x.size(); i++)
	{
		RS_WORD xi = this->gf.powTable[this->gf.characteristic - coefPos[i]]; //TODO determine if equivalent to GaloisField::Inv(x[i])
		RS_WORD errLocPrime = 1;
		for (int j = 0; j < x.size(); j++)
		{
			if (j != i)
			{
				errLocPrime = this->gf.mult(errLocPrime, 1 ^ this->gf.mult(xi, x[j]));
			}
		}
		if (errLocPrime == 0)
		{
			return false;
		}
		RS_WORD y = this->gf.mult(x[i], Poly_Eval(&errEval, xi, &this->gf)); //errEval is still reversed from earlier
		//TODO determine if equivalent to GaloisField::Mult(GaloisField::Pow(xi, 1), y)
		e.coef[errPos->at(i)] = this->gf.div(y, errLocPrime); //magnitude
	}
	Poly_Add(msg, msg, &e);
	return true;
}

bool ReedSolomon::findErrorLocator(Poly* out, Poly* synd, int nsym, Poly* eraseLoc, int eraseCount)
{
	//this spits out a polynomial in reverse order but i dont know why
	RS_WORD init = 1;
	Poly errLoc(1, &init);
	Poly oldLoc(1, &init);
	Poly temp;
	if (eraseLoc)
	{
		errLoc.setCopy(eraseLoc->n, eraseLoc->coef);
		oldLoc.setCopy(eraseLoc->n, eraseLoc->coef);
	}
	int syndShift = 0; //TODO optimize
	for (int i = nsym - eraseCount - 1; i >= 0; i--)
	{
		int K = i + syndShift + eraseCount;
		RS_WORD delta = synd->coef[K];
		for (int j = 1; j < errLoc.n; j++)
		{
			delta ^= this->gf.mult(errLoc.coef[errLoc.n - j - 1], synd->coef[K + j]);
		}
		Poly_Pad(&oldLoc, 0, 1); //TODO optimize
		if (delta != 0)
		{
			if (oldLoc.n > errLoc.n)
			{ //TODO isn't this always the case?
				Poly_Scale(&temp, &oldLoc, delta, &this->gf);
				Poly_Scale(&oldLoc, &errLoc, this->gf.inv(delta), &this->gf);
				errLoc.setCopy(temp.n, temp.coef);
			}
			Poly_Scale(&temp, &oldLoc, delta, &this->gf);
			Poly_Add(&errLoc, &errLoc, &temp);
		}
	}
	int leading = 0;
	for (; errLoc.coef[leading] == 0; leading++);
	Poly_Trim(&errLoc, leading, 0);
	int errs = errLoc.n - 1;
	out->setCopy(errLoc.n, errLoc.coef);
	if (errs * 2 - eraseCount > nsym)
	{
		return false;
	}
	return true;
}

bool ReedSolomon::findErrors(vector<unsigned int>* out, Poly* errLoc, int n)
{
	int errs = errLoc->n - 1;
	Poly revErrLoc;
	Poly_Reverse(&revErrLoc, errLoc);
	if (errLoc->n == 1)
	{
		//do something special here? idk
	}
	else if (errLoc->n == 2)
	{ //linear equation
		out->push_back(this->gf.logTable[this->gf.div(errLoc->coef[0], errLoc->coef[1])]);
	} 
	else
	{
		Poly_ChienSearch(out, &revErrLoc, n, &this->gf);
	}
	if (out->size() != errs)
	{
		// Too many (or few) errors found by Chien Search for the errata locator polynomial!
		return false;
	}
	//map to string pos
	for (RS_WORD i = 0; i < out->size(); i++)
	{
		if (out->at(i) >= n) //clearly something messed up
		{
			return false;
		}
		(*out)[i] = n - out->at(i) - 1;
	}
	return true;
}

void ReedSolomon::forneySyndromes(Poly* out, Poly* synd, vector<unsigned int>* pos, int n)
{
	Poly fsynd(synd->n - 1, synd->coef);
	if (pos)
	{
		for (unsigned int i : *pos)
		{
			RS_WORD rev = (RS_WORD)n - i - 1;
			RS_WORD x = this->gf.powTable[rev];
			for (int j = fsynd.n - 2; j >= 0; j--)
			{
				fsynd.coef[j + 1] = this->gf.mult(fsynd.coef[j + 1], x) ^ fsynd.coef[j];
			}
		}
	}
	//fsynd.coef[fsynd.n - 1] = 0;
	out->setCopy(fsynd.n, fsynd.coef);
}

bool ReedSolomon::decode(RS_WORD* wholeOut, RS_WORD* out, RS_WORD* data, int k, int nsym, vector<unsigned int>* erasePos, bool debug)
{
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
			for (unsigned int i : *erasePos)
			{
				msg.coef[i] = 0;
			}
		}
	}
	this->calcSyndromes(&synd, &msg, nsym);
	if (!this->checkSyndromes(&synd))
	{
		if(debug) cout << "errors detected, locating" << endl;
		Poly fsynd, errLoc;
		this->forneySyndromes(&fsynd, &synd, erasePos, k + nsym);
		bool canLocate = this->findErrorLocator(&errLoc, &fsynd, nsym, nullptr, erasePos ? erasePos->size() : 0);
		if (!canLocate)
		{
			if (debug) cout << "too many errors to locate!" << endl;
			return false;
		}
		vector<unsigned int> pos;
		canLocate = this->findErrors(&pos, &errLoc, k + nsym);
		if (!canLocate || !(pos.size() || (erasePos && erasePos->size())))
		{
			if(debug) cout << "errors unable to be located!" << endl;
			return false;
		}
		if (debug)
		{
			if (pos.size())
			{
				cout << "additional errors detected at ";
				for_each(pos.begin(), pos.end(), [](unsigned int e) {cout << (int)e << ", "; });
			}
			cout << "correcting" << endl;
		}
		if (erasePos)
		{
			pos.insert(pos.begin(), erasePos->begin(), erasePos->end());
		}
		bool success = this->correctErrata(&msg, &synd, &pos);
		if (!success)
		{
			if(debug) cout << "decode failure!" << endl;
			return false;
		}
		if(debug) cout << "errors corrected" << endl;
	}
	if (wholeOut)
	{
		memcpy(wholeOut, msg.coef, sizeof(RS_WORD) * (k + nsym));
	}
	if (out)
	{
		memcpy(out, msg.coef, sizeof(RS_WORD) * k);
	}
	return true;
}