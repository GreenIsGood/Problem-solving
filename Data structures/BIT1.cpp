#include <bits/stdc++.h>
using namespace std;
const int N = (1 << 17);

typedef long long ll;

struct BIT {
	ll M[N], C[N];

	void init(int n) {
		n = (1 << (int) ceil(log2(n)));
		memset(M, 0, n * sizeof(M[0]));
		memset(C, 0, n * sizeof(C[0]));
	}

	void add(int idx, ll valM, ll valC) {
		++idx;
		while (idx <= N) {
			M[idx - 1] += valM;
			C[idx - 1] += valC;
			idx += (idx & (-idx));
		}
	}
	ll get(int idx) {
		ll valM = 0, valC = 0;
		int x = idx;
		++idx;
		while (idx) {
			valM += M[idx - 1];
			valC += C[idx - 1];
			idx -= (idx & (-idx));
		}
		return valM * x + valC;
	}

	void addrange(int s, int e, ll val) {
		add(s, val, (1 - s) * val);
		add(e + 1, -val, e * val);
	}

} bit;

int main() {

	int tc;
	scanf("%d", &tc);
	while (tc--) {
		int n, c;
		scanf("%d%d", &n, &c);
		bit.init(n);
		int t, x, y;
		ll v;
		while (c--) {
			scanf("%d", &t);
			if (t) {
				scanf("%d%d", &x, &y);
				printf("%lld\n", bit.get(y) - bit.get(x - 1));
			} else {
				scanf("%d%d%lld", &x, &y, &v);
				bit.addrange(x, y, v);
			}
		}

	}
}
