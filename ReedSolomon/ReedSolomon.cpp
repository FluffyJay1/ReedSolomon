#include "ReedSolomon.h"
using namespace std;

int primes[] = {
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
	unsigned int fieldCharacteristic = (1 << fieldPower) - 1, fieldCharacteristicNext = (1 << (fieldPower + 1)) - 1;
	for (unsigned int i = fieldCharacteristic + 2; (fieldCharacteristicNext == 0 ? i > 0 : i < fieldCharacteristicNext); i += 2)
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
	this->characteristic = (1 << fieldPower) - 1;
	this->fieldPower = fieldPower;
	this->primitivePoly = primes[fieldPower];
	//init the tables 
	unsigned int val = 1;
	this->powTable = (unsigned int*)malloc(sizeof(int) * (this->characteristic + 1) * 2);
	this->logTable = (unsigned int*)malloc(sizeof(int) * (this->characteristic + 1));
	powTable[0] = val;
	logTable[0] = 0;
	logTable[1] = 0;
	for (unsigned int i = 1; i < this->characteristic; i++)
	{
		val <<= 1;
		if (val > this->characteristic)
		{
			val ^= this->primitivePoly;
		}
		powTable[i] = (unsigned int)val;
		logTable[(unsigned int)val] = i;
	}
	for (unsigned int i = this->characteristic; i < (this->characteristic + 1) * 2; i++)
	{
		powTable[i] = powTable[i - this->characteristic];
	}
	/*
	for (unsigned int i = 0; i < this->characteristic; i++)
	{
		cout << "2^" << i << "=" << powTable[i] << endl;
	}
	for (unsigned int i = 0; i < this->characteristic; i++)
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

unsigned int GaloisField::multNoLUT(int a, int b)
{
	unsigned int ret = 0;
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

inline unsigned int GaloisField::mult(unsigned int a, unsigned int b)
{
	return (a == 0 || b == 0) ? 0 : this->powTable[this->logTable[a] + this->logTable[b]];
}
inline unsigned int GaloisField::div(unsigned int a, unsigned int b)
{
	return a == 0 ? 0 : (b == 0) ? -1 : this->powTable[this->logTable[a] - this->logTable[b] + this->characteristic];
}

inline unsigned int GaloisField::pow(unsigned int x, unsigned int power)
{
	return this->powTable[(this->logTable[x] * power) % this->characteristic];
}

inline unsigned int GaloisField::inv(unsigned int x)
{
	return this->powTable[this->characteristic - this->logTable[x]];
}

inline unsigned int GaloisField::sqrt(unsigned int x)
{
	return logTable[x] % 2 ? this->powTable[(logTable[x] + this->characteristic) / 2] : this->powTable[logTable[x] / 2];
}

Poly::Poly()
{
	this->init();
}

Poly::Poly(int n, unsigned int* data)
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

void Poly::setCopy(int n, unsigned int* coef)
{
	if (n > this->n)
	{
		if (this->coef)
		{
			free(this->coef);
		}
		this->coef = (unsigned int*)malloc(sizeof(int) * n);
	}
	this->n = n;
	if (coef)
	{
		memcpy(this->coef, coef, sizeof(int) * n);
	} else
	{
		memset(this->coef, 0, sizeof(int) * n);
	}
}

void Poly::setRef(int n, unsigned int* coef)
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
		cout << "[" << setw(3) << hex << (int)this->coef[0];
	}
	for (int i = 1; i < this->n; i++)
	{
		cout << ", " << setw(3) << hex << (int)this->coef[i];
	}
	if (this->n > 0)
	{
		cout << "]" << dec << endl;
	}
}

Poly* Poly_Create(int n, unsigned int* coef)
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
	unsigned int* temp = (unsigned int*)malloc(sizeof(int) * n);
	memset(temp, 0, sizeof(int) * n);
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

void Poly_Scale(Poly* out, Poly* in, unsigned int scale, GaloisField* gf)
{
	if (out == in)
	{
		for (int i = 0; i < in->n; i++)
		{
			in->coef[i] = gf->mult(in->coef[i], scale);
		}
	} else
	{
		unsigned int* temp = (unsigned int*)malloc(sizeof(int) * in->n);
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
	unsigned int* temp = (unsigned int*)malloc(sizeof(int) * n);
	memset(temp, 0, sizeof(int) * n);
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
	unsigned int* temp = (unsigned int*)malloc(sizeof(int)* a->n);
	unsigned int normalizer = b->coef[0];
	memcpy(temp, a->coef, sizeof(int) * a->n);
	for (int i = 0; i < a->n - b->n + 1; i++)
	{
		temp[i] = gf->div(temp[i], normalizer);
		unsigned int coef = temp[i];
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

unsigned int Poly_Eval(Poly* poly, unsigned int x, GaloisField* gf)
{
	unsigned int y = poly->coef[0];
	for (int i = 1; i < poly->n; i++)
	{
		y = gf->mult(y, x) ^ poly->coef[i];
	}
	return y;
}

void Poly_ChienSearch(vector<unsigned int>* out, Poly* poly, int max, GaloisField* gf)
{
	//this seems unnecessary because all multiplications are performed via lookup table anyway
	unsigned int* temp = (unsigned int*)malloc(sizeof(int)* poly->n);
	memcpy(temp, poly->coef, sizeof(int) * poly->n);
	for (int i = 0; i < max; i++)
	{
		unsigned int sum = 0;
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
	unsigned int* temp = (unsigned int*)malloc(sizeof(int)* n);
	memset(temp, 0, sizeof(int) * left);
	memcpy(temp + left, poly->coef, sizeof(int) * poly->n);
	memset(temp + (left + poly->n), 0, sizeof(int) * right);
	poly->setRef(n, temp);
}

void Poly_Trim(Poly* poly, int left, int right)
{
	int n = poly->n - left - right;
	unsigned int* temp = (unsigned int*)malloc(sizeof(int)* n);
	memcpy(temp, poly->coef + left, sizeof(int) * n);
	poly->setRef(n, temp);
}

void Poly_Append(Poly* out, Poly* a, Poly* b)
{
	int n = a->n + b->n;
	unsigned int* temp = (unsigned int*)malloc(sizeof(int)* n);
	memcpy(temp, a->coef, sizeof(int)* a->n);
	memcpy(temp + a->n, b->coef, sizeof(int)* b->n);
	out->setRef(n, temp);
}

void Poly_Reverse(Poly* out, Poly* in)
{
	unsigned int* temp = (unsigned int*)malloc(sizeof(int)* in->n);
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

void ReedSolomon::encode(unsigned int* out, unsigned int* data, int k, int nsym)
{
	Poly msg(k, data);
	Poly generator, remainder;
	this->createGenerator(&generator, nsym);
	Poly_Pad(&msg, 0, nsym);
	Poly_Div(nullptr, nullptr, &remainder, &msg, &generator, &this->gf);
	memcpy(out, data, sizeof(int) * k);
	memcpy(out + k, remainder.coef, sizeof(int) * remainder.n);
}

void ReedSolomon::calcSyndromes(Poly* out, Poly* msg, int nsym)
{
	unsigned int* temp = (unsigned int*)malloc(sizeof(int)* (nsym + 1));
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
	for (int i : *errPos)
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
	for (int i : *errPos)
	{
		coefPos.push_back(msg->n - 1 - i);
	}
	Poly errLoc, errEval;
	this->findErrataLocator(&errLoc, &coefPos);
	this->findErrorEvaluator(&errEval, synd, &errLoc, errLoc.n); //TODO determine if correct
	//Poly_Reverse(errEval, errEval); //reverse it for later use
	vector<unsigned int> x(coefPos.size());
	for (unsigned int i = 0; i < x.size(); i++)
	{
		x[i] = this->gf.powTable[coefPos[i]]; //TODO determine if correct
	}
	Poly e(msg->n, nullptr);
	for (unsigned int i = 0; i < x.size(); i++)
	{
		unsigned int xi = this->gf.powTable[this->gf.characteristic - coefPos[i]]; //TODO determine if equivalent to GaloisField::Inv(x[i])
		unsigned int errLocPrime = 1;
		for (unsigned int j = 0; j < x.size(); j++)
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
		unsigned int y = this->gf.mult(x[i], Poly_Eval(&errEval, xi, &this->gf)); //errEval is still reversed from earlier
		//TODO determine if equivalent to GaloisField::Mult(GaloisField::Pow(xi, 1), y)
		e.coef[errPos->at(i)] = this->gf.div(y, errLocPrime); //magnitude
	}
	Poly_Add(msg, msg, &e);
	return true;
}

bool ReedSolomon::findErrorLocator(Poly* out, Poly* synd, int nsym, Poly* eraseLoc, int eraseCount)
{
	//this spits out a polynomial in reverse order but i dont know why
	unsigned int init = 1;
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
		unsigned int delta = synd->coef[K];
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
	for (unsigned int i = 0; i < out->size(); i++)
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
			unsigned int rev = n - i - 1;
			unsigned int x = this->gf.powTable[rev];
			for (int j = fsynd.n - 2; j >= 0; j--)
			{
				fsynd.coef[j + 1] = this->gf.mult(fsynd.coef[j + 1], x) ^ fsynd.coef[j];
			}
		}
	}
	//fsynd.coef[fsynd.n - 1] = 0;
	out->setCopy(fsynd.n, fsynd.coef);
}

bool ReedSolomon::decode(unsigned int* wholeOut, unsigned int* out, unsigned int* data, int k, int nsym, vector<unsigned int>* erasePos, bool debug)
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
		memcpy(wholeOut, msg.coef, sizeof(int) * (k + nsym));
	}
	if (out)
	{
		memcpy(out, msg.coef, sizeof(int) * k);
	}
	return true;
}