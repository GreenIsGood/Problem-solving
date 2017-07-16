#include <bits/stdc++.h>

using namespace std;

const int N = 1e5 + 5;

vector<int> adj[N];
int n, DFSN[N], LOW[N], id, root, rootCalls;
set<int> Arti;
vector<pair<int, int>> Bridges;

void init(){
  for(int i=0; i<N; i++)
    adj[i].clear();
  id = 0;
  memset(DFSN, 0, sizeof DFSN);
  memset(LOW, 0, sizeof LOW);
  Arti.clear();
  Bridges.clear();
}

void Tarjan(int u, int p = -1){
  DFSN[u] = LOW[u] = ++id;
  for(int v: adj[u]){
    if(!DFSN[v]){
      Tarjan(v, u);
      rootCalls += u == root;
      LOW[u] = min(LOW[u], LOW[v]);
      if(DFSN[u] < LOW[v]) Bridges.push_back({min(u, v), max(u, v)});
      if(DFSN[u] <= LOW[v] && root != u) Arti.insert(u);
    }else if(p != v)
      LOW[u] = min(LOW[u], DFSN[v]);
  }
}

int main()
{
  while(~scanf("%d", &n)){
    init();
    for(int i=0; i<n; i++){
      int u, k, v;
      scanf("%d (%d)", &u, &k);
      while(k--){
        scanf("%d", &v);
        adj[u].push_back(v);
      }
    }
    for(int i=0; i<n; i++){
      if(!DFSN[i]){
        root = i;
        rootCalls = 0;
        Tarjan(i);
        if(rootCalls > 1)
          Arti.insert(root);
      }
    }
    sort(Bridges.begin(), Bridges.end());
    printf("%d critical links\n", Bridges.size());
    for(auto i: Bridges)
      printf("%d - %d\n", i.first, i.second);
    puts("");
  }
  return 0;
}
