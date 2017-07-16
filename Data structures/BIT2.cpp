#include <bits/stdc++.h>
using namespace std;
const int N = (1 << 4);

typedef long long ll;

struct BIT {
	ll C[N];
	int n;

	void init(int n) {

		this->n = n = (1 << (int) ceil(log2(n)));

		memset(C, 0, n * sizeof(C[0]));
	}

	void add(int idx, ll valC) {
		++idx;
		while (idx <= n) {

			C[idx - 1] += valC;
			idx += (idx & (-idx));
		}
	}
	ll get(int idx) {
		ll valC = 0;
		int x = idx;
		++idx;
		while (idx) {

			valC += C[idx - 1];
			idx -= (idx & (-idx));
		}
		return valC;
	}

	int find(ll val) {
		int st = 0;
		for (int mid = (n >> 1); mid; mid >>= 1) {
			int idx = st + mid - 1;
			ll cur = C[idx];
			if (val > cur) {
				st += mid;
				val -= cur;
			}
		}
		return st;
	}

} bit;

int main() {

//	int tc;
//	scanf("%d", &tc);
//	while (tc--) {
//		int n, c;
//		scanf("%d%d", &n, &c);
//		bit.init(n);
//		int t, x, y;
//		ll v;
//		while (c--) {
//			scanf("%d", &t);
//			if (t) {
//				scanf("%d%d", &x, &y);
//				printf("%lld\n", bit.get(y) - bit.get(x - 1));
//			} else {
//				scanf("%d%d%lld", &x, &y, &v);
//				bit.addrange(x, y, v);
//			}
//		}
//
//	}
	bit.init(15);
	bit.add(5, 2);
	bit.add(6, 2);
	bit.add(7, 2);

	for (int i = 0; i < 9; i++)
		cout << i << " " << bit.find(i) << "\n";

}
