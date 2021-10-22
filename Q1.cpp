#include <bits/stdc++.h>
#define int long long
#define MP make_pair
#define PA pair<int,int>
#define fi first
#define se second
using namespace std;
const int maxl=200005,maxg=22;

int read(){
	int s=0,w=1;char ch=getchar();
	while (ch<'0'||ch>'9'){if (ch=='-')  w=-w;ch=getchar();}
	while (ch>='0'&&ch<='9'){s=(s<<1)+(s<<3)+(ch^'0');ch=getchar();}
	return s*w;
}
int n,q,cnt;
int head[maxl],w[maxl],fa[maxl];PA ans[maxl];
int maxpos[maxl],maxv[maxl],maxdis[maxl];

struct edge{int nxt,to,dis;}e[maxl<<1];
struct node{int u,v,c,t;}b[maxl];
struct Query{int u,x,id;}lis[maxl];

bool cmpx(Query tmpx,Query tmpy){return tmpx.x>tmpy.x;}
bool cmpt(node tmpx,node tmpy){return tmpx.t>tmpy.t;} 
void chkmax(int x,int &y){y=max(x,y);}
//Tree works
namespace Tree{
	int lg[maxl],dep[maxl],fa[maxl][maxg],f[maxl][maxg];
	void init(){
		for (int i=1;i<=n;i++)  lg[i]=log2(i)+1;
	}
	void dfs(int now,int fath){
		dep[now]=dep[fath]+1,fa[now][0]=fath;
		for (int i=1;i<=lg[dep[now]];i++){
			fa[now][i]=fa[fa[now][i-1]][i-1];
			f[now][i]=max(f[now][i-1],f[fa[now][i-1]][i-1]);
		}
		for (int i=head[now];i;i=e[i].nxt){
			int y=e[i].to;
			if (y!=fath)  continue;
			f[y][0]=e[i].dis,dfs(y,now);
		}
	}
	int LCA(int x,int y){
		if (dep[x]<dep[y])  swap(x,y);
		while (dep[x]>dep[y])  x=fa[x][lg[dep[x]-dep[y]]-1];
		if (x==y)  return x;
		
		for (int k=lg[dep[x]];k>=0;k--){
			if (fa[x][k]^fa[y][k])  x=fa[x][k],y=fa[y][k];
		}
		return fa[x][0];
	}
	int query_lian(int u,int v){
		int d=dep[u]-dep[v],res=0;
		for (int i=lg[dep[u]];i>=0;i--){
			if (d&(1<<i)){
				res=max(res,f[u][i]);
				u=fa[u][i];
			}
		}
		return res;
	}
	int query(int x,int y){
		int l=LCA(x,y);
		return max(query_lian(x,l),query_lian(y,l));
	}
	void work(){init(),dfs(1,0);}
}
//DSU


int query(int u, int v) {
        if (dep[u] > dep[v])
            swap(u, v);
        for (int i = LOG - 1; i >= 0; i--) {
            if (((dep[v] - dep[u]) >> i) & 1)
                v = table[i][v];
        }
        if (u == v)
            return u;
        for (int i = LOG - 1; i >= 0; i--) {
            if (table[i][u] != table[i][v]) {
                u = table[i][u];
                v = table[i][v];
            }
        }
        return table[0][u];
    }

    int ancestor(int u, int k) {
        if (k < 0)
            return u;
        int res = u;
        for (int i = 0; i < LOG; i++)
            if ((k >> i) & 1)
                res = table[i][res];
        return res;
    }

    int distance(int u, int v) {
        return dep[u] + dep[v] - 2 * dep[query(u, v)];
    }

    int depth(int u) {
        return dep[u];
    }

    int move(int u, int v, int k) {
        int anc = query(u, v);
        if (anc == u)
            return ancestor(v, dep[v] - dep[u] - k);
        else if (dep[u] - dep[anc] >= k)
            return ancestor(u, k);
        else
            return ancestor(v, dep[u] + dep[v] - 2 * dep[anc] - k);
    }

    int par(int u) {
        return table[0][u];
    }
};
using LCA = DoublingLowestCommonAncestor<UnWeightedGraph>;
// ----- library -------




//main work
namespace ducati{
	void add_edge(int u,int v,int w){
		cnt++;
		e[cnt].to=v,e[cnt].dis=w,e[cnt].nxt=head[u],head[u]=cnt;
	}
	void get_all_in(){
		n=read(),q=read();
		for (int i=1;i<=n;i++)  w[i]=read();
		for (int i=1;i<n;i++){
			b[i].u=read(),b[i].v=read(),b[i].t=read(),b[i].c=read();
			add_edge(b[i].u,b[i].v,b[i].c);
			add_edge(b[i].v,b[i].u,b[i].c);
		}
		for (int i=1;i<=q;i++)  lis[i].id=i,lis[i].x=read(),lis[i].u=read();
	}
	void work(){
		sort(lis+1,lis+q+1,cmpx),sort(b+1,b+n,cmpt);
		int r=1;
		int i=1;
		while(i<=q){
			while (r<=n&&b[r].t>=lis[i].x){
				DSU::Merge(b[r].u,b[r].v);
				r++;
			}
			int rt=DSU::Find(lis[i].u),res1,res2;
			res1=maxv[rt];
			res2=max(Tree::query(lis[i].u,maxpos[rt]),maxdis[rt]);
			ans[lis[i].id]=MP(res1,res2);
			i++;
		}
	}
	void print_ans(){
		for (int i=1;i<=q;i++)  printf("%d %d\n",ans[i].fi,ans[i].se);
	}
	void solve(){
		get_all_in(),Tree::work(),DSU::init();
		work(),print_ans();
	}
}

signed main(){
	ducati::solve();
	return 0;
}
