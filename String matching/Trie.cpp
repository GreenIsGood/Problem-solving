#include<bits/stdc++.h>

using namespace std;

const int mx = 1e5 + 10;

int child[mx][128];
bool isEnd[mx];
int childCount[mx];
int nodeCount;

int addNode() {
	memset(child[nodeCount], -1, sizeof child[nodeCount]);
	isEnd[nodeCount] = false;
	childCount[nodeCount] = 0;
	return nodeCount++;
}

void init() {
	nodeCount = 0;
	addNode();
}

void insert(char *s) {
	int cur = 0;
	for (; *s; s++) {
		int &nxt = child[cur][(int) *s];
		if (nxt == -1)
			nxt = addNode(), childCount[cur]++;
		cur = nxt;
	}
	isEnd[cur] = true;
}

int main() {
	int tc;
	scanf("%d", &tc);
	while (tc--) {
		int n;
		init();
		scanf("%d", &n);
		while (n--) {
			char s[15];
			scanf("%s", s);
			insert(s);
		}
		for (int i = 0; i < nodeCount; i++) {
			if (isEnd[i] && childCount[i]) {
				puts("NO");
				goto barra;
			}
		}
		puts("YES");
		barra: ;
	}
}

