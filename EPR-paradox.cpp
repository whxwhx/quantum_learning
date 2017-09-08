#include <cstdio>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#define rep(i,a,b) for(int i = a; i <= b; i++)
#define dep(i,a,b) for(int i = a; i >= b; i--)
#define Rep(i,a) for(int i = 0; i < a; i ++)
extern "C"
{
	#include <quantum.h>
}
using namespace std;

void classical() {
	double best_probability = 0;
	Rep(da0,2) Rep(da1,2) Rep(db0,2) Rep(db1,2) {
		int cnt = 0;
		rep(i,1,1000000) {
			int x = rand() & 1, y = rand() & 1;
			int a = (x == 0) ? da0 : da1;
			int b = (y == 0) ? db0 : db1;
			if ((a ^ b) == (x & y)) cnt++;
		}
		best_probability = max(best_probability, cnt / 1000000.0);
	}
	printf("Classical best probability : %.5lf\n", best_probability);
}

void init_bell(quantum_reg &share) {
	share = quantum_new_qureg(0, 2);
	quantum_hadamard(0, &share);
	quantum_cnot(0, 1, &share);
}

const double pi = acos(-1.0);
typedef unsigned long long ull;
void quantum() {
	int cnt = 0, t1 = 0, c1 = 0;
	rep(i,1,1000000) {
		quantum_reg share;
		init_bell(share);
		
		int x = rand() & 1, y = rand() & 1;
		if (x) quantum_r_y(0, pi / 4, &share);
		if (y) quantum_r_y(1,- pi / 4, &share);


		ull s = quantum_measure(share);	
		int a = s & 1;
		int b = (s >> 1) & 1;
		
		
		if ((a ^ b) == (x & y)) cnt++;
/*		if (x == 1 && y == 1) {
			t1++;
			if ((a ^ b) == (x & y)) c1++;
		}
*/
		quantum_delete_qureg(&share);
	}
//	cerr <<c1 * 1.0 / t1<<' '<<c1<<' '<<t1<<endl;
	printf("Quantum probability : %.5lf\n", cnt / 1000000.0);
}

int main ()
{
	srand(time(NULL));
	classical();
	quantum();
	return 0;
}


