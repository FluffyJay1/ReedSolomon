#include "ReedSolomon.h"
#include <time.h>
#include <iomanip>

using namespace std;

void Poly_Print(Poly* poly)
{
	cout << "Poly(n=" << poly->n << ")";
	if (poly->n > 0)
	{
		cout << "[" << setw(3) << (int)poly->coef[0];
	}
	for (int i = 1; i < poly->n; i++)
	{
		cout << ", " << setw(3) << (int)poly->coef[i];
	}
	if (poly->n > 0)
	{
		cout << "]" << endl;
	}
}

int main()
{
	int k = 10, nsym = 7, ncorrupt = 3;
	unsigned char data1[] = { 0x40, 0xd2, 0x75, 0x47, 0x76, 0x17, 0x32, 0x06, 0x27, 0x26, 0x96, 0xc6, 0xc6, 0x96, 0x70, 0xec };
	Poly msg(k, data1);
	Poly a(k + nsym, data1);
	RS_Encode(a.coef, data1, k, nsym);
	cout << "Original message: " << endl;
	Poly_Print(&msg);
	cout << endl << "Encoded message: " << endl;
	Poly_Print(&a);
	cout << endl;
	vector<unsigned char> possPos, corrPos;
	for (int i = 0; i < k + nsym; i++)
	{
		possPos.push_back(i);
	}
	srand(time(0));
	for (int i = 0; i < ncorrupt; i++)
	{
		int randInd = rand() % possPos.size();
		vector<unsigned char>::iterator it = next(possPos.begin(), randInd);
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
	Poly_Print(&a);
	cout << endl;
	bool success = RS_Decode(a.coef, msg.coef, a.coef, k, nsym, nullptr, true);
	if (!success)
	{
		cout << "Decoding failed!" << endl;
	}
	cout << "After decoding: " << endl;
	Poly_Print(&a);
	cout << endl << "Decoded message: " << endl;
	Poly_Print(&msg);
}