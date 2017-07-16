#include <bits/stdc++.h>

using namespace std;

const int N = 1e5 + 5;

vector<int> adj[N], comp[N], adjComp[N];
int n, DFSN[N], LOW[N], inStack[N], id, f[N], cnt, inDeg[N], outDeg[N];
stack<int> st;

void init(){
  for(int i=0; i<N; i++){
    comp[i].clear();
    adj[i].clear();
    adjComp[i].clear();
  }
  cnt = id = 0;
  memset(inDeg, 0, sizeof inDeg);
  memset(outDeg, 0, sizeof outDeg);
  memset(DFSN, 0, sizeof DFSN);
  memset(LOW, 0, sizeof LOW);
}


void Extract(int u){
  int node = -1;
  while(node != u){
    node = st.top();
    st.pop();
    comp[cnt].push_back(node);
    f[node] = cnt;
    inStack[node] = 0;
  }
  ++cnt;
}

void Tarjan(int u){
  DFSN[u] = LOW[u] = ++id;
  inStack[u] = 1;
  st.push(u);
  for(int v: adj[u]){
    if(!DFSN[v]){
      Tarjan(v);
      LOW[u] = min(LOW[u], LOW[v]);
    }else if(inStack[v]){
      LOW[u] = min(LOW[u], DFSN[v]);
    }
  }
  if(LOW[u] == DFSN[u]) Extract(u);
}

int T, m;


int main()
{
  scanf("%d", &T);
  while(T--){
    init();
    scanf("%d %d", &n, &m);
    for(int i=0; i<m; i++){
      int u, v;
      scanf("%d %d", &u, &v);
      adj[u].push_back(v);
    }
    for(int i=1; i<=n; i++){
      if(!DFSN[i]){
        Tarjan(i);
      }
    }
    for(int u=1; u<=n; u++){
      for(int v: adj[u]){
        int cu = f[u];
        int cv = f[v];
        if(cu == cv) continue;
        adjComp[cu].push_back(cv);
        ++inDeg[cv];
        ++outDeg[cu];
      }
    }
    int ans = 0;
    for(int i=0; i<cnt; i++)
      ans += !inDeg[i];
    printf("%d\n", ans);
  }
  return 0;
}
