#pragma once

#include<graph_v_of_v_idealID_DPBF_only_ec.cuh>
/*this is the DPBF algorithm in Ding, Bolin, et al. "Finding top-k min-cost connected trees in databases." 2007 IEEE 23rd International Conference on Data Engineering. IEEE, 2007.

time complexity: O( 4^|Gamma| + 3^|Gamma||V|+ 2^|Gamma|* (|E| + |V|*(|Gamma| + log V)) )*/


const int N=2e3+10;
const double inf=1e18;
struct node{
    int statu,pos;double cost;
    node(int _statu=0,int _pos=0,double _cost=0):statu(_statu),pos(_pos),cost(_cost){}
    bool operator<(const node &k)const{
        return cost>k.cost;
    }
};
vector<vector<node> >pre;
vector<int>num;
vector<vector<double> >dp;
vector<vector<bool> >used;
vector<vector<int> >idp;

void getpath(int statu,int x,graph_hash_of_mixed_weighted& solu){
    if(pre[statu][x].statu==-1)return;
    int _statu=pre[statu][x].statu;
    int _x=pre[statu][x].pos;
    if(_statu==statu){
        graph_hash_of_mixed_weighted_add_edge(solu, x, _x, pre[statu][x].cost);
        getpath(_statu,_x,solu);
    }
    else{
        getpath(_statu,x,solu);
        getpath(statu-_statu,x,solu);
    }
}
__global__ void grow(int x,int statu,graph_v_of_v_idealID& input_graph,priority_queue<node>&q){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i<input_graph[x].size()){
        int y=input_graph[x][i].first;
        double w=input_graph[x][i].second;
        if(dp[statu][y]>dp[statu][x]+w){
            dp[statu][y]=dp[statu][x]+w;
            pre[statu][y]=node(statu,x,w);
            q.push(node(statu,y,dp[statu][y]));
        }
    }
}

__global__ void merge(int x,int statu,priority_queue<node>&q){
    int pos = blockIdx.x * blockDim.x + threadIdx.x;
    if(pos<idp[statu].size()){
        int i=idp[statu][pos];
        if(dp[statu|i][x]>dp[statu][x]+dp[i][x]){
            dp[statu|i][x]=dp[statu][x]+dp[i][x];
            pre[statu|i][x]=node(i,x,0);
            q.push(node(statu|i,x,dp[statu|i][x]));
        }
    }
}

graph_hash_of_mixed_weighted graph_v_of_v_idealID_DPBF_only_ec(
	graph_v_of_v_idealID& input_graph, graph_v_of_v_idealID& group_graph, std::unordered_set<int>& cumpulsory_group_vertices, double& RAM_MB) {
	int cnt=0;
    unordered_map<int,int>id;
    const int n=(int)input_graph.size();
    const int m=(1<<cnt)-1;

	for(auto it=cumpulsory_group_vertices.begin();it!=cumpulsory_group_vertices.end();it++){
		id[*it]=cnt++;
	}
    cudaMallocManaged((void**)&pre, (m+1)*n*sizeof(node));
    cudaMallocManaged((void**)&num, n*sizeof(int));
    cudaMallocManaged((void**)&dp, (m+1)*n*sizeof(double));
    cudaMallocManaged((void**)&used, (m+1)*n*sizeof(bool));
    cudaMallocManaged((void**)&idp, (m+1)*(m+1)*sizeof(int));
	pre.resize(m+1);
    num.resize(n);
    dp.resize(m+1);
    used.resize(m+1);
    for(int i=0;i<=m;i++){
        pre[i].resize(n);
        dp[i].resize(n);
        used[i].resize(n);
    }
    idp.resize(m+1);
    for(int i=0;i<=m;i++){
        int lim=m-i,z=0;
        for(int j=lim;j;j=(j-1)&lim)
            idp[i].push_back(j);
    }
	for(int x=0;x<n;x++){
		for(int i=0;i<group_graph[x].size();i++){
			int y=group_graph[x][i].first;
			num[x]|=1<<id[y];
		}
	}
	
	for(int i=0;i<=m;i++)
		for(int j=0;j<n;j++)
			dp[i][j]=inf,used[i][j]=0,pre[i][j]=node(-1,0,0);
    for(int i=0;i<n;i++)dp[0][i]=0;
    priority_queue<node>q;
    for(int i=0;i<n;i++)if(num[i]){
        dp[num[i]][i]=0;
        q.push(node(num[i],i,0));
    }
    int threadsPerBlock = 1024;
    int numBlocks;
	while(q.size()){
        node now=q.top();q.pop();
        if(used[now.statu][now.pos])continue;
        used[now.statu][now.pos]=1;
        int x=now.pos,statu=now.statu;
        if(statu==m){
			graph_hash_of_mixed_weighted solu;
            getpath(statu,x,solu);
            cudaFree(pre);
            cudaFree(num);
            cudaFree(dp);
            cudaFree(used);
			return solu;
        }
        //grow
        numBlocks = ((int)input_graph[x].size() + threadsPerBlock - 1) / threadsPerBlock;
        grow<<< numBlocks, threadsPerBlock >>>(x,statu,input_graph,q);
        cudaDeviceSynchronize();
        //merge
		numBlocks = ((int)idp[statu].size() + threadsPerBlock - 1) / threadsPerBlock;
        merge<<< numBlocks, threadsPerBlock >>>(x,statu,q);
        cudaDeviceSynchronize();
    }
}
