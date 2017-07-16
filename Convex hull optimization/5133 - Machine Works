
C++ code for Machine Works, Author fegla

/*
 TASK:
 LANG: C++
 */
#include <algorithm>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <iterator>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <list>
#include <bitset>
#include <cstring>
#include <climits>
#include <deque>
#include <utility>
#include <complex>
#include <numeric>
#include <functional>
#include <stack>
#include <iomanip>
#include <ctime>
#include <valarray>
//#include "vc.h"
#ifdef _MSC_VER
#include <hash_set>
#include <hash_map>
using namespace stdext;
#else
#if __GNUC__ > 2
#include <ext/hash_set>
#include <ext/hash_map>
using namespace __gnu_cxx;
#else
#include <hash_set>
#include <hash_map>
#endif
#endif

template<class key>
struct hashtemp {

        enum {
                bucket_size = 4, min_buckets = 8
        };
        virtual size_t operator()(const key &p) const=0;

};

using namespace std;
#ifndef M_PI
const long double M_PI=acos((long double)-1);
#endif
#define rep(i,n) for(int  i=0;i<(int)(n);++i)
long double ZERO = 0;
const long double INF = 1 / ZERO, EPSILON = 1e-12;
#define all(c) (c).begin(),(c).end()
#define rep2(i,a,b) for(int i=(a);i<=((int)b);++i)
#define foreach(it,c) for(__typeof((c).begin()) it=(c).begin();it!=(c).end();++it)

#define sz(v) ((int)((v).size()))
#define let(x,y) __typeof(y) x(y)

#define rrep(i,n) for(int  i=((int)n)-1;i>=0;--i)
#define rall(c) (c).rbegin(),(c).rend()
#define rrep2(i,a,b) for(int i=(a);i>=((int)b);--i)
#define rforeach(it,c) for(__typeof((c).rbegin()) it=(c).rbegin();it!=(c).rend();++it)
#define rep2d(i, j, v) rep(i, sz(v)) rep(j, sz(v[i]))
#define foreach2d(i, j, v) foreach(i,v) foreach(j,*i)
typedef long long ll;

struct machine {
        long long d, p, r, g;
        bool operator<(const machine &m) const {
                return d < m.d;
        }
} machines[(int) 1e5 + 1];
long long mem[(int) 1e5 + 1];
//pair<int,pair<long long,long long> >  before[(int)1e5+1];
set<pair<ll, ll> > active;

inline bool checkValid(const pair<ll, ll> & line1, const pair<ll, ll> & line2,
                const pair<ll, ll> & line3) {
//      d=(C2-C1)/(M1-M2)
        double i1 = (line3.second - line1.second)
                        / (line1.first - line3.first * 1.0);
        double i2 = (line3.second - line2.second)
                        / (line2.first - line3.first * 1.0);
        if (i1 > i2)
                return 0;
        return 1;
}

int main() {
        std::ios_base::sync_with_stdio(false);
#ifndef ONLINE_JUDGE

        freopen("works.in", "rt", stdin);
//      freopen("1.in", "rt", stdin);
        freopen("1.out", "wt", stdout);
#endif

        int n, tt = 1, d;
        long long c;
        while (cin >> n >> c >> d, n || c || d) {
                active.clear();
                rep(i,n) {
                        machine &m = machines[i];
                        cin >> m.d >> m.p >> m.r >> m.g;
                }
                machines[n].d = d + 1;
                sort(machines, machines + n);
                rep(i,n+1) {
                        long long &r = mem[i];
                        r = c;
                        machine &c = machines[i];
                        while (active.size() >= 2) {
                                pair<ll, ll> line1 = *active.begin();
                                pair<ll, ll> line2 = *(++active.begin());
                                ll v1 = line1.second + line1.first * c.d;
                                ll v2 = line2.second + line2.first * c.d;
                                if (v1 < v2) {
                                        active.erase(active.begin());
                                } else {
                                        r = max(r,v1);
                                        break;
                                }
                        }
                        r = max(r,(*active.begin()).second + (*active.begin()).first * c.d);
                        set<pair<ll, ll> >::iterator it1, it2, it3, cur;
                        if (r >= c.p) {
                                //fun(day)=mem[i]-cur.p+cur.r+cur.g*(Day-cur.d-1ll);
                                //fun(day)=mem[i]-cur.p+cur.r+cur.g*Day-cur.g*cur.d-cur.g;
                                //fun(day)=M1day+C1
                                long long C1 = mem[i] - c.p + c.r - c.g * c.d - c.g;
                                long long M1 = c.g;
                                cur = it1 = it2 = it3 = active.insert(make_pair(M1, C1)).first;
                                --it1;
                                ++it3;

                                // If a parallel line exist before, remove
                                if (it2 != active.begin() && it1->first == it2->first) {
                                        active.erase(it1);
                                }
                                it1 = it2;
                                --it1;

                                // If a parallel line exist after, remove
                                if (it3 != active.end() && it2->first == it3->first) {
                                        active.erase(it2);
                                        continue;
                                }

                                if (it2 != active.begin() && it3 != active.end()) {
                                        if (!checkValid(*it1, *it2, *it3)) {
                                                active.erase(it2);
                                                continue;
                                        }
                                }

                                while (1) {
                                        it1 = it2 = it3 = cur;
                                        ++it2;
                                        if (it2 == active.end())
                                                break;
                                        ++(++it3);
                                        if (it3 == active.end())
                                                break;
                                        if (!checkValid(*it1, *it2, *it3)) {
                                                active.erase(it2);
                                        } else {
                                                break;
                                        }
                                }

                                while (1) {
                                        it1 = it2 = it3 = cur;
                                        if (it3 == active.begin())
                                                break;
                                        --it2;
                                        if (it2 == active.begin())
                                                break;
                                        --(--it1);
                                        if (!checkValid(*it1, *it2, *it3)) {
                                                active.erase(it2);
                                        } else {
                                                break;
                                        }
                                }
                        }

                }
                cout << ("Case ") << tt++ << ": " << mem[n] << endl;
        }
        return 0;
}
