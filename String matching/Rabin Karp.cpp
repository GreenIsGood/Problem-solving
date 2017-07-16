#include <bits/stdc++.h>
using namespace std;

#define MAX (3 * (int) (1e4) + 4)
#define ll long long
#define MOD 1000000007
#define base 1009

int k;
char s[MAX];
ll hsh[MAX];
ll rhsh[MAX];

ll fast_pow(ll a) {
    ll ans = 1;
    ll base1 = base;
    while(a) {
        if(a & 1) {
            ans *= base1;
            ans %= MOD;
        }
        base1 *= base1;
        base1 %= MOD;
        a >>= 1;
    }
    return ans;
}
int hash_() {
    int n = strlen(s);
    ll h1 = 0, h2 = 0;
    for(int i = 0, j = n - 1; i < n; i++, j--) {
        h1 = (h1 * base) % MOD;
        h1 = (h1 + s[i]) % MOD;
        h2 = (h2 * base) % MOD;
        h2 = (h2 + s[j]) % MOD;
        hsh[i] = h1;
        rhsh[j] = h2;
    }
    ll ans = 0;
    ll v1, v2;
    for(int i = 0, j = k - 1; j < n; i++, j++) {
        if(i) {
            v1 = (hsh[j] - (hsh[i - 1] * fast_pow(k)) % MOD + MOD) % MOD;
        }else {
            v1 = hsh[j];
        }
        if(j + 1 == n) {
            v2 = rhsh[i];
        }else {
            v2 = (rhsh[i] - (rhsh[j + 1] * fast_pow(k)) % MOD + MOD) % MOD;
        }
        //cout << v1 << " " << v2 << endl;
        if(v1 == v2) {
            ans++;
        }
    }
    return ans;
}

int main() {
    scanf("%d", &k);
    scanf("%s", s);
    cout << hash_() << endl;
}
