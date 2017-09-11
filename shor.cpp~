#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
extern "C"{
	#include <quantum.h>
}
#define rep(i,a,b) for(int i = a; i <= b; i++)
#define dep(i,a,b) for(int i = a; i >= b; i--)
#define Rep(i,a) for(int i = 0; i < a; i ++)
using namespace std;
int n;

int pw(int a, int b, int mod) {
	int w = 1; 
	for(;b;b >>= 1, a = 1LL * a * a % mod) if (b & 1) w = 1LL * w * a % mod;
	return w; 
}

bool check(int a, int n) {
	int r = n - 1;
	while (!(r & 1)) r >>= 1;
	int t = pw(a, r, n); if (t == 1) return true;
	while (r < n - 1) {
		if (t == n - 1) return true;
		r <<= 1;
		t = 1LL * t * t % n;
	}
	return false;
}

bool miller_rabin(int n) {
	if (n == 2) return true;
	rep(i,1,40) {
		int base = rand() % (n - 2) + 2;
		if (!check(base, n)) return false;
	}
	return true;
}

vector<int> ans;
#define pb(a) push_back(a)

bool compare(int a, int b, int n) { //a ^ b < n
	int w = 1;
	while (b) {
		if (b & 1) { 
			if (n / w < a) return false;
			w = w * a;
		}
		b >>= 1; if (n / a < a && b) return false; a = a * a;
	}
	return w < n;
}

int gcd(int a, int b) { return b ? gcd(b, a % b) : a; }

int find(int n, int i) {
	int l = 0, r = n;
	while (l + 1 < r) {
		int mid = (l + r) >> 1;
		if (compare(mid, i, n)) l = mid; else r = mid;
	}
	return !pw(r, i, n) ? r : 0;
}

int Log(int n) {
	int j = 0;
	while ((1 << j) < n) j++;
	return j;
}

#define id(x) (x + l)
#define mp(a,b) make_pair(a,b)
typedef pair<int, int> pii;

typedef long long LL;


void Swap(int a, int b, quantum_reg &reg) {
	quantum_cnot(a, b, &reg);
	quantum_cnot(b, a, &reg);
	quantum_cnot(a, b, &reg);
}

const double pi = acos(-1.0);

void qft(int m, quantum_reg &reg) {
	Rep(i,m / 2) Swap(i, m - 1 - i, reg);
	Rep(i,m) {
		Rep(j,i) quantum_cond_phase_kick(i, j, 2 * pi / (1 << (i + 1)) * (1 << j), &reg);//use controled gate!
		quantum_hadamard(i, &reg);
	}
}

const long double eps = 1e-10;
inline int approximate(int a, int b, int n) {
	long double t = ((long double)(a)) / b;
	LL q0 = 0, q1 = 1;
	t -= int(t); 
	while (q1 < n && t) {
		t = 1 / t;
		LL q2 = 1LL * int(t) * q1 + q0; 
		q0 = q1, q1 = q2;
		t -= int(t); 
		if (t < eps) break;
	}
	return q1 < n ? q1 : q0;
}

int find_order(int a, int n) {
	int t = Log(n); 
	int m = 2 * t - 1, M = 1 << m;
	quantum_reg reg;
	reg = quantum_new_qureg(0, m);
	quantum_walsh(m, &reg);
	quantum_addscratch(3*t+2, &reg);
	quantum_exp_mod_n(n, a, m, t, &reg);
	Rep(i,3 * t + 2) quantum_bmeasure(0, &reg);
//	quantum_print_qureg(reg);
	qft(m, reg);
//	quantum_print_qureg(reg);
	int x = quantum_measure(reg);
	quantum_delete_qureg(&reg);
	return approximate(x, M, n);
}


void shor(int n) {
	if (n == 1 || n == 2) return;
	if (n % 2 == 0) { cout <<n<<" = "<<n / 2<<'*'<<2<<endl; ans.pb(2); shor(n / 2); return; }
	for(int i = 1; (1 << i) <= n; i++) {
		int root = find(n, i);
		if (root && miller_rabin(root)) {
			cout <<n<<" = "; 
			rep(j,1,i) cout <<root<<(j == i ? '\n' : '*'), ans.pb(root);
			return;
		}
	}
	bool flag = false;
	for(;!flag;) {
		int a = rand() % (n - 2) + 2;
		int g = gcd(a, n);
		if (g > 1) {
			cout <<n<<" = "<<g<<'*'<<n/g<<endl;
			flag = true;
			shor(g); shor(n / g);
		} else {
			int r = find_order(a, n);
			if (r > 0 && !(r & 1) && pw(a, r, n) == 1 && pw(a, r/2,n) != 1 && pw(a, r/2,n) != n-1) {
				r >>= 1; int t = pw(a, r, n);
				if (t != -1) {
					int g1 = gcd(t - 1, n);
					int g2 = gcd(t + 1, n);
					int g3 = gcd(g1, g2);
					g1 /= g3, g2 /= g3;
					flag = true;
					cout <<n<<" = "<<g1<<'*'<<g2<<'*'<<g3<<'*'<<n / g1 / g2 / g3<<endl;
					shor(g1), shor(g2), shor(g3), shor(n / g1 / g2 / g3);
				}
			}
		}
	}
}

int main() {
	srand(233);
	scanf("%d",&n);
	shor(n);
	sort(ans.begin(), ans.end());
	cout <<"Final result : "<<n<<" = ";
	int l = ans.size();
	Rep(j,l) cout <<ans[j]<<(j == l - 1 ? '\n' : '*');
	return 0;
}

