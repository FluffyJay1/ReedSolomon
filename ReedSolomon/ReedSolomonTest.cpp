#include "ReedSolomon.h"
#include <time.h>
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

void test()
{
	int bits = 8, k = 16, nsym = 10;
	RS_WORD data1[] = { 0x40, 0xd2, 0x75, 0x47, 0x76, 0x17, 0x32, 0x06, 0x27, 0x26, 0x96, 0xc6, 0xc6, 0x96, 0x70, 0xec };
	Poly msg(k, data1);
	Poly a(k + nsym, data1);
	ReedSolomon rs(bits);
	rs.encode(a.coef, data1, k, nsym);
	cout << "Original message: " << endl;
	msg.print();
	cout << endl << "Encoded message: " << endl;
	a.print();
	cout << endl;
	vector<unsigned int> corrPos, erasePos;
	erasePos.push_back(1);
	erasePos.push_back(2);
	erasePos.push_back(3);
	corrPos.push_back(4);
	corrPos.push_back(5);
	corrPos.push_back(6);
	srand(time(0));
	for (unsigned char i : erasePos)
	{
		a.coef[i] = 0;
	}
	for (unsigned char i : corrPos)
	{
		a.coef[i] = 1;
	}
	cout << endl << "Corrupted message: " << endl;
	a.print();
	cout << endl;
	bool success = rs.decode(a.coef, msg.coef, a.coef, k, nsym, &erasePos, true);
	if (!success)
	{
		cout << "Decoding failed!" << endl;
	} else
	{
		cout << "After decoding: " << endl;
		a.print();
		cout << endl << "Decoded message: " << endl;
		msg.print();
	}
}

bool test(int bits, int k, int nsym, int ncorr, bool print, bool* unmatch, int erasures)
{
	RS_WORD *originalData = (RS_WORD*)malloc(sizeof(RS_WORD) * (k + nsym));
	for (int i = 0; i < k; i++)
	{
		originalData[i] = rand() % (1 << bits);
	}
	Poly msg(k, originalData);
	Poly a(k + nsym, originalData);
	ReedSolomon rs(bits);
	rs.encode(a.coef, originalData, k, nsym);
	if (print)
	{
		cout << "Original message: " << endl;
		msg.print();
		cout << endl << "Encoded message: " << endl;
		a.print();
		cout << endl;
	}
	vector<unsigned int> possPos, corrPos, erasePos;
	for (int i = 0; i < k + nsym; i++)
	{
		possPos.push_back(i);
	}
	for (int i = 0; i < ncorr; i++)
	{
		int randInd = rand() % possPos.size();
		vector<unsigned int>::iterator it = next(possPos.begin(), randInd);
		corrPos.push_back(possPos[randInd]);
		possPos.erase(it);
	}
	if (print) cout << "Corrupting indices ";
	for (unsigned char i : corrPos)
	{
		a.coef[i] = rand() % (1 << bits);
		if (print) cout << (int)i << " ";
	}
	if (print)
	{
		cout << endl << "Corrupted message: " << endl;
		a.print();
		cout << endl;
	}
	for (int i = 0; i < erasures && corrPos.size() > 0; i++)
	{
		int randInd = rand() % corrPos.size();
		vector<unsigned int>::iterator it = next(corrPos.begin(), randInd);
		erasePos.push_back(corrPos[randInd]);
		corrPos.erase(it);
	}
	if (print)
	{
		cout << "Erasures at indices ";
		for (unsigned char i : erasePos)
		{
			cout << (int)i << " ";
		}
		cout << endl;
	}
	bool success = rs.decode(a.coef, msg.coef, a.coef, k, nsym, erasures ? &erasePos : nullptr, print);
	if (print)
	{
		if (!success)
		{
			cout << "Decoding failed!" << endl;
		} else
		{
			cout << "After decoding: " << endl;
			a.print();
			cout << endl << "Decoded message: " << endl;
			msg.print();
		}
	}
	for (int i = 0; i < k; i++)
	{
		if (msg.coef[i] != originalData[i])
		{
			if (unmatch)
			{
				*unmatch = true;
			}
			break;
		}
	}
	return success;
}

void testConfig(int k, int nsym, int bits)
{
	srand(time(0));
	for (int errors = 0; errors <= k + nsym; errors++)
	{
		int successes = 0, trials = 0, falsePositives = 0, truePositives = 0;
		for (int i = 0; i < 1000; i++)
		{
			bool fp = false;
			bool success = test(bits, k, nsym, errors, false, &fp, 0);
			if (success)
			{
				successes++;
			}
			if (success)
			{
				if (fp)
				{
					falsePositives++;
				} else
				{
					truePositives++;
				}
			}
			trials++;
		}
		cout << "E(" << errors << "):" << "decoded: " << ((double)successes / trials) << ", fp:" << ((double)falsePositives / trials)
			<< ", rs:" << ((double)truePositives / trials) << endl;
	}
}

void testFalsePositiveRate(int k, int maxNsym, int bits)
{
	for (int nsym = 2; nsym <= maxNsym; nsym++)
	{
		int successes = 0, trials = 0, falsePositives = 0, truePositives = 0;
		for (int i = 0; i < 1000; i++)
		{
			bool fp = false;
			bool success = test(bits, k, nsym, k + nsym, false, &fp, 0);
			if (success)
			{
				successes++;
			}
			if (success)
			{
				if (fp)
				{
					falsePositives++;
				} else
				{
					truePositives++;
				}
			}
			trials++;
		}
		cout << "E(" << nsym << "):" << "decoded: " << ((double)successes / trials) << ", fp:" << ((double)falsePositives / trials)
			<< ", rs:" << ((double)truePositives / trials) << endl;
	}
}

int main()
{
	srand(time(NULL));
	//test(3, 6, 2, 1, true, nullptr, 0);
	testConfig(2, 4, 3);
	//testFalsePositiveRate(12, 9, 6);
	//test();
	//FindPrimePolys(&cout, 5, 10);
}