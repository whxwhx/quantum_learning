//a sigle solution 3-CNF-SAT solver using Grover's searching algorithm 
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
extern "C"
{
	#include <quantum.h>
}
#define rep(i,a,b) for(int i = a; i <= b; i++)
#define dep(i,a,b) for(int i = a; i >= b; i--)
#define Rep(i,a) for(int i = 0; i < a; i++)
using namespace std;
const int N = 50 + 5;
const double pi = acos(-1.0);
int n, m, id[N][3], f[N][3];

typedef unsigned long long ull;

void quantum_or(int a, int b, int c, quantum_reg &reg) { //c ^= a | b
	quantum_sigma_x(a, &reg);
	quantum_sigma_x(b, &reg);
	quantum_toffoli(a, b, c, &reg);
	quantum_sigma_x(c, &reg);
	quantum_sigma_x(a, &reg);
	quantum_sigma_x(b, &reg);
}

#define pos(x) (x + n)
#define id(j) (pos(m + 3 + j))
void cal_f(quantum_reg &reg) {	
	rep(i,1,m) {
		int x = pos(m + 1), y = pos(m + 2);
		Rep(j,3) {
			quantum_cnot(id[i][j], id(j), &reg);
			if (f[i][j]) quantum_sigma_x(id(j), &reg);
		}
		quantum_or(id(0), id(1), x, reg);
		quantum_or(x, id(2), y, reg);
		if (i == 1) quantum_cnot(y, pos(i), &reg);
		else quantum_toffoli(pos(i - 1), y, pos(i), &reg);
		quantum_or(x, id(2), y, reg);
		quantum_or(id(0), id(1), x, reg);
		Rep(j,3) {
			if (f[i][j]) quantum_sigma_x(id(j), &reg);
			quantum_cnot(id[i][j], id(j), &reg);
		}
	}
	
	quantum_cnot(pos(m), n, &reg);
	
	//erase scratchpad
	dep(i,m,1){
		int x = pos(m + 1), y = pos(m + 2);
		Rep(j,3) {
			quantum_cnot(id[i][j], id(j), &reg);
			if (f[i][j]) quantum_sigma_x(id(j), &reg);
		}
		quantum_or(id(0), id(1), x, reg);
		quantum_or(x, id(2), y, reg);
		if (i == 1) quantum_cnot(y, pos(i), &reg);
		else quantum_toffoli(pos(i - 1), y, pos(i), &reg);
		quantum_or(x, id(2), y, reg);
		quantum_or(id(0), id(1), x, reg);
		Rep(j,3) {
			quantum_cnot(id[i][j], id(j), &reg);
			if (f[i][j]) quantum_sigma_x(id(j), &reg);
		}
	}
	
}

void cal_g(quantum_reg &reg) {
	rep(i,1,n) {
		quantum_sigma_x(i - 1, &reg);
		if (i > 1) quantum_toffoli(pos(i - 1), i - 1, pos(i), &reg);
		else quantum_cnot(i - 1, pos(i), &reg);
		quantum_sigma_x(i - 1, &reg);
	}
	quantum_cnot(pos(n), n, &reg);
	quantum_sigma_x(n, &reg);
	dep(i,n,1) {
		quantum_sigma_x(i - 1, &reg);
		if (i > 1) quantum_toffoli(pos(i - 1), i - 1, pos(i), &reg);
		else quantum_cnot(i - 1, pos(i), &reg);
		quantum_sigma_x(i - 1, &reg);
	}
}

int main() {
	freopen("SAT.in","r",stdin);
	srand(time(NULL));
	scanf("%d",&n); //n variables
	scanf("%d",&m); //m clauses
	rep(i,1,m) Rep(j,3) scanf("%d%d",&f[i][j], &id[i][j]); //f : 0->x_id 1->!x_id
	
	cout <<"Quantum sigle solution 3-CNF-SAT Solver"<<endl;
	rep(i,1,m) {
		cout <<"(";
		Rep(j,3) {
			if (f[i][j]) cout <<"!"; 
			cout <<"x_"<<id[i][j];
			if (j < 2) cout <<" or ";
		}
		cout <<")"; if (i < m) cout <<" and ";
	}
	cout <<endl;
	
	int cnt = n + max(n, m + 6); cout <<cnt<<" qubits required"<<endl;
	
	rep(t,0,100) {	
		quantum_reg reg;
		reg = quantum_new_qureg(0, cnt);
		quantum_walsh(n, &reg);
	
		int times = 1; rep(i,1,(n + 1) / 2) times <<= 1;
		times = times * pi / 2 + 10;
	
		rep(i,1,times) {
			//reflect around e
			cal_f(reg);
			quantum_sigma_z(n, &reg);
			cal_f(reg);
		
			quantum_walsh(n, &reg);
			cal_g(reg);
			quantum_sigma_z(n, &reg);
			cal_g(reg);
			quantum_walsh(n, &reg);	
		}
	
		ull S = quantum_measure(reg);
		static bool x[N];
		Rep(i,n) x[i] = (S >> i & 1);
		bool ok = true;
		rep(i,1,m) {
			bool flag = false;
			Rep(j,3) flag |= f[i][j] ? !x[id[i][j]] : x[id[i][j]];
			if (!flag) { cout <<"failed #"<<t<<endl; ok = false; break; }
		}
		if (ok) {
			Rep(i,n) cout <<"x_"<<i<<" : "<<x[i]<<endl;
			return 0;
		}
		
		quantum_delete_qureg(&reg);
	}
	cout <<"Can't find any solution"<<endl;
	return 0;
}
