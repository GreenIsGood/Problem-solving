#include <bits/stdc++.h>
using namespace std;

void computeLPSArray(string pat, int M, int *lps) {
    int len = 0;
    lps[0] = 0;
    int i = 1;
    while(i < M) {
        if(pat[i] == pat[len]) {
            len++;
            lps[i] = len;
            i++;
        }else {
            if(len != 0) {
                len = lps[len-1];
            }else {
                lps[i] = 0;
                i++;
            }
        }
    }
}

void KMPSearch(string pat,string txt,int *lps) {
    int M=pat.size();
    int N=txt.size();
    bool flg=0;
    int i=0;
    int j=0;
    while(i<N) {
        if(pat[j]==txt[i]) {
            j++;
            i++;
        }
        if(j==M) {
            flg=1;
            printf("%d\n", i-j);
            j = lps[j-1];
        }else if(i<N&&pat[j]!=txt[i]) {
            if(j!=0) {
                j=lps[j-1];
            }else {
                i++;
            }
        }
    }
    if(!flg) {
        puts("");
    }
}

int main() {
    int n;
    string pat;
    string txt;
    while(cin>>n) {
        cin>>pat>>txt;
        int lps[pat.size()];
        memset(lps,0,sizeof lps);
        computeLPSArray(pat,pat.size(),lps);
        KMPSearch(pat,txt,lps);
    }
}
