#include <bits/stdc++.h>
using namespace std;

#define ll long long
#define MAX (5 * (int) 1e4 + 4)
#define F first
#define S second

ll mem[MAX];
pair<ll, ll> input[MAX];
set<pair<ll, ll>, greater<pair<ll, ll> > > act;
vector<pair<ll, ll> > fix;

bool checkValid(const pair<ll, ll> &L1, const pair<ll, ll> &L2, const pair<ll, ll> &L3) {
    double inters1 = (L3.S - L1.S) / (L1.F - L3.F * 1.0);
    double inters2 = (L3.S - L2.S) / (L2.F - L3.F * 1.0);
    if(inters1 > inters2) {
        return true;
    }
    return false;
}

int main () {
    for(int i = 0; i < MAX; i++) {
        mem[i] = (ll)1e18;
    }
    int n;
    scanf("%d", &n);
    for(int i = 0; i < n; i++) {
        scanf("%lld %lld", &input[i].F, &input[i].S);
    }
    fix.push_back({0, 0});
    sort(input, input + n);
    for(int i = n - 1; i >= 0; i--) {
        if(!input[i].F) {
            continue;
        }
        fix.push_back({input[i].F, input[i].S});
        for(int j = i - 1; j >= 0; j--) {
            if(input[j].S > input[i].S) {
                i = j + 1;
                break;
            }
            input[j].F = 0;
            input[j].S = 0;
        }
    }

    sort(fix.begin(), fix.end());

    for(int i = 0; i < fix.size(); i++) {
        ll &r = mem[i];
        set<pair<ll, ll>, greater<pair<ll, ll> > >::iterator it1, it2, it3, cur;
        while(act.size() >= 2) {
            it1 = it2 = it3 = cur = act.begin();
            (++it2);
            ll line1 = it1->F * fix[i].F + it1->S;
            ll line2 = it2->F * fix[i].F + it2->S;
            if(line1 >= line2) {
                act.erase(it1);
            }else {
                break;
            }
        }
        if(act.size()) {
            it1 = it2 = it3 = cur = act.begin();
            r = min(r, it1->F * fix[i].F + it1->S);
        }else {
            r = 0;
        }

        ll m = fix[i + 1].S;
        ll c = mem[i];
        it1 = it2 = it3 = cur = act.insert({m, c}).F;

        (++it2);
        if(it2 != act.end() && it1->F == it2->F) {
            act.erase(it1);
            continue;
        }

        it1 = it2 = cur;
        (--it2);
        if(it1 != act.begin() && it1->F == it2->F) {
            act.erase(it2);
        }

        it2 = cur;
        (--it1);
        (++it3);
        if(it2 != act.begin() && it3 != act.end()) {
            if(!checkValid(*it1, *it2, *it3)) {
                act.erase(it2);
                continue;
            }
        }

        while(1) {
            it1 = it2 = it3 = cur;
            (++it2);
            if(it2 == act.end()) {
                break;
            }
            ++(++it3);
            if(it3 == act.end()) {
                break;
            }
            if(!checkValid(*it1, *it2, *it3)) {
                act.erase(it2);
            }else {
                break;
            }
        }

        while(1) {
            it1 = it2 = it3 = cur;
            if(it1 == act.begin()) {
                break;
            }
            (--it2);
            if(it2 == act.begin()) {
                break;
            }
            --(--it3);
            if(!checkValid(*it1, *it2, *it3)) {
                act.erase(it2);
            }else {
                break;
            }
        }
    }
    printf("%lld", mem[fix.size() - 1]);
}
        /*for(int j = i; j >= 0; j--) {
            mem[i] = min(mem[i], mem[j - 1] + input[i].F * input[j].S);
        }*/
