#include "ReedSolomon.h"
#include <time.h>
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

void test()
{
	int bits = 10, k = 10, nsym = 7, ncorrupt = 10;
	unsigned int data1[] = { 0x40, 0xd2, 0x75, 0x47, 0x76, 0x17, 0x32, 0x06, 0x27, 0x26, 0x96, 0xc6, 0xc6, 0x96, 0x70, 0xec };
	Poly msg(k, data1);
	Poly a(k + nsym, data1);
	ReedSolomon rs(bits);
	rs.encode(a.coef, data1, k, nsym);
	cout << "Original message: " << endl;
	msg.print();
	cout << endl << "Encoded message: " << endl;
	a.print();
	cout << endl;
	vector<unsigned int> possPos, corrPos;
	for (int i = 0; i < k + nsym; i++)
	{
		possPos.push_back(i);
	}
	srand(time(0));
	for (int i = 0; i < ncorrupt; i++)
	{
		int randInd = rand() % possPos.size();
		vector<unsigned int>::iterator it = next(possPos.begin(), randInd);
		corrPos.push_back(possPos[randInd]);
		possPos.erase(it);
	}
	cout << "Corrupting indices ";
	for (unsigned char i : corrPos)
	{
		a.coef[i] = 0;
		cout << (int)i << " ";
	}
	cout << endl << "Corrupted message: " << endl;
	a.print();
	cout << endl;
	bool success = rs.decode(a.coef, msg.coef, a.coef, k, nsym, nullptr, true);
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

int main()
{
	test();
	//FindPrimePolys(&cout, 16, 10);
}