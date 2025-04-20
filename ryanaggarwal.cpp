#include <bits/stdc++.h>
using namespace std;
struct Edge {
    int from ,to, cost;
    bool operator<(const Edge& other) const {
        return cost < other.cost;
    }
};
const int INF = 1e9;

int numNodes;
map<int,int>map_ifhub;
map<int,int>map_ifhouse;
map<int,int>map_house_hub;
vector<int>visited_hubs;
vector<int>visited_houses;

vector<bool> visited;
vector<int>parent;
vector<int>path;
vector<int>vec, vec2;

struct DisjointSetUnion {
    vector<int > parent, rank;
    DisjointSetUnion(int  n) {
        parent.resize(n);
        rank.resize(n, 0);
        for(int i = 0; i <n; i++)
        {
            parent[i] = i;
        }
    }
    int  find(int  u) {
        if (parent[u] != u)
            parent[u] = find(parent[u]);
        return parent[u];
    }
    bool unite(int  u, int  v) {
        int  pu = find(u), pv = find(v);
        if (pu == pv) return false;
        if (rank[pu] < rank[pv])
            parent[pu] = pv;
        else {
            parent[pv] = pu;
            if (rank[pu] == rank[pv])
                rank[pu]++;
        }
        return true;
    }
};

void dijkstra(vector<int>&fuelStations, vector<vector<pair<int, int>>> &mstGraph, vector<pair<int, int>> &dist) {
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> pq;

    for(int i = 0; i <fuelStations.size(); i++)
    {
        pq.push({0,fuelStations[i]});
        dist[fuelStations[i]].second = fuelStations[i];
        dist[fuelStations[i]].first = 0;
    }
    
    while (!pq.empty()) {
        int u = pq.top().second;
        int d = pq.top().first;
        pq.pop();
        
        if (d > dist[u].first) continue;
        
        for (auto &element : mstGraph[u]) {

            int v = element.first, weight = element.second;
            if (dist[u].first + weight < dist[v].first) {
                dist[v].first = dist[u].first + weight;
                dist[v].second = dist[u].second;
                pq.push({dist[v].first, v});
            }
        }
    }
}

vector<vector<int>> distiii;
void floydWarshall(int V, vector<vector<int>>& distiii) {
    for (int k = 0; k < V; ++k) {
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) {
                if (distiii[i][k] < INF && distiii[k][j] < INF)
                    distiii[i][j] = min(distiii[i][j], distiii[i][k] + distiii[k][j]);
            }
        }
    }
}


void dfsCheck(int vertex,vector<bool>&visited, vector<bool>& finall,vector<vector<pair<int,int>>> &minimumSpanningTree
    ,int fuelstation,vector<pair<int,int>>&mindist, set<int> &houseHub, vector<int>&subtree){

    visited[vertex] = true;
    for( auto& pr : minimumSpanningTree[vertex])
    {
        int child=pr.first;
        if(visited[child])continue;
        int fuelstation_child=mindist[child].second;
        if(fuelstation_child!=fuelstation)continue;
        dfsCheck(child,visited, finall,minimumSpanningTree,fuelstation,mindist, houseHub, subtree);
        subtree[vertex]+= subtree[child];
    }
    if (subtree[vertex] == 0)
    {
        finall[vertex] =true;
    }
}

bool dfs(int vertex,vector<bool>&visited,vector<vector<pair<int,int>>> &minimumSpanningTree
    ,int &ans,int fuelstation,vector<pair<int,int>>&mindist,int &currentfuel){
    
    vec.push_back(vertex);
    visited[vertex]=1;
    int cnt=0;
    if(map_ifhub[vertex]==1){
        visited_hubs[vertex]=1;
    }
    if(map_ifhouse[vertex]==1 and visited_hubs[map_house_hub[vertex]]==1){
        visited_houses[vertex]=1;
    }
    bool done=1;
    for(auto pr:minimumSpanningTree[vertex]){
        int child=pr.first;
        int wt=pr.second;
        if(visited[child])continue;

        int fuelstation_child=mindist[child].second;
        if(fuelstation_child!=fuelstation)continue;

        if(currentfuel-wt<mindist[child].first){
            done=0;
        }
        else{
            ans+=2*wt;
            currentfuel -= wt;
            done=((done)&dfs(child,visited, minimumSpanningTree,ans,fuelstation,mindist,(currentfuel)));
            cnt++;
        }
    }
    if(done==0){
        visited[vertex]=0;
    }
    if (cnt>0)
    {
        vec.push_back(vertex);
    }

    

    return done;
}

void dfs2(int vertex,vector<bool>&visited2,vector<vector<pair<int,int>>> &minimumSpanningTree,int &ans, int fuelCapacity){
    vec2.push_back(vertex);
    visited2[vertex]=1;
    bool done=1;
    int cnt=0;
    for(auto pr:minimumSpanningTree[vertex]){
        int child=pr.first;
        int wt=pr.second;
        if(visited2[child])continue;
        if(wt>fuelCapacity){
            done=0;
            continue;
        }
        else{
            ans+=2*wt;
            cnt++;
            dfs2(child,visited2, minimumSpanningTree,ans,fuelCapacity);
        }
    }
    if (cnt>0)
    {
        vec2.push_back(vertex);
    }
}
void dfs3(int vertex,vector<vector<pair<int,int>>> &minimumSpanningTree,vector<int>&visited,vector<int>&parent){
    visited[vertex]=1;
    for(auto pr:minimumSpanningTree[vertex]){
        int child=pr.first;
        if(visited[child])continue;
        parent[child]=vertex;
        dfs3(child,minimumSpanningTree,visited,parent);
    }
    return;
}
void dfs4(int vertex,vector<bool>&visited,vector<vector<pair<int , int >>> &minimumSpanningTree,vector<int>&fuelStations,int &sum){
    visited[vertex]=1;
    vec.push_back(vertex);
    for(auto pr:minimumSpanningTree[vertex]){
        if(visited[pr.first]==1)continue;
        if(pr.first==fuelStations[0]){
            sum+=pr.second*2;
            return;
        }

        dfs4(pr.first,visited,minimumSpanningTree,fuelStations,sum);
        sum+=pr.second;

    }
    return;
 }

int  main() {
    ios::sync_with_stdio(0);
    cin.tie(0);

    // freopen("input3.txt", "r", stdin);
    // freopen("output1.txt", "w", stdout);

    int  numDeliveries, numRoads, numFuelStations, fuelCapacity;
    cin >> numDeliveries >> numNodes >> numRoads >> numFuelStations >> fuelCapacity;
    parent.assign(numNodes,-1);
    visited_hubs.assign(numNodes,0);
    visited_houses.assign(numNodes,0);
    visited.assign(numNodes,false);
    distiii.assign(numNodes, vector<int> (numNodes,INF));
    set<int> s, houseAndHubs;
    vector<int > deliveryHubs(numDeliveries);
    for (int  i = 0; i < numDeliveries; ++i){ 
        cin >> deliveryHubs[i];
        houseAndHubs.insert(deliveryHubs[i]);
        s.insert(deliveryHubs[i]);
    }

    vector<int > deliveryHouses(numDeliveries);
    for (int  i = 0; i < numDeliveries; ++i){ 
        cin >> deliveryHouses[i];
        houseAndHubs.insert(deliveryHouses[i]);
        s.insert(deliveryHouses[i]);
    }

    vector<int > fuelStations(numFuelStations);
    for (int  i = 0; i < numFuelStations; ++i){ 
        cin >> fuelStations[i];
        s.insert(fuelStations[i]);
    }

    vector<Edge> edges;
    for (int  i = 0; i < numRoads; ++i) {
        int  from, to, cost;
        cin >> from >> to >> cost;
        distiii[from][to] = cost;
        distiii[to][from] = cost;

        edges.push_back({from, to, cost});
    }
        for (int i = 0; i < numNodes; ++i) distiii[i][i] = 0;
    floydWarshall(numNodes,distiii);

///
    vector<Edge> edges2;

    for (int i = 0; i < numFuelStations; i++)
    {
        for(int j = i+1; j < numFuelStations; j++)
        {
            if (i!= j && distiii[fuelStations[j]][fuelStations[i]] != INF)
            {
                edges2.push_back({i, j, distiii[fuelStations[j]][fuelStations[i]]});
            }
        }
    }
    
    sort(edges2.begin(), edges2.end());
    DisjointSetUnion dsu2(numFuelStations);
    vector<vector<pair<int , int >>> minimumSpanningTree2(numFuelStations);
    int sum2 = 0;
    //cout<<edges2.size()<<endl;
    for (const auto &edge : edges2) {
        if (dsu2.unite(edge.from, edge.to)) {
            // sum2+= edge.cost;
            // if (edge.cost>fuelCapacity)
            // {
               // cout<<edge.from<<" "<<edge.to<<endl;
            // }
            minimumSpanningTree2[edge.from].push_back({edge.to, edge.cost});
            minimumSpanningTree2[edge.to].push_back({edge.from, edge.cost});
        }
    }
        vector<bool> visited2(numFuelStations,false);
    for(int i = 0; i < numFuelStations; i++)
    {
        if (visited2[i])
        {
            continue;
        }
        int capacity_forone = fuelCapacity;
        dfs2(i,visited2, minimumSpanningTree2,sum2, capacity_forone);
    }

/////


    sort(edges.begin(), edges.end());
    DisjointSetUnion dsu(numNodes);
    vector<vector<pair<int , int >>> minimumSpanningTree(numNodes);

    for (const auto &edge : edges) {
        if (dsu.unite(edge.from, edge.to)) {
            minimumSpanningTree[edge.from].push_back({edge.to, edge.cost});
            minimumSpanningTree[edge.to].push_back({edge.from, edge.cost});
        }
    }
    vector<pair<int,int>> minDist(numNodes);

    for(int i = 0; i <numNodes; i++)
    {
        minDist[i].first = INT_MAX;
        minDist[i].second = 0;
    }
    dijkstra(fuelStations,minimumSpanningTree,minDist);

    for(int i=0;i<deliveryHubs.size();i++){
        map_ifhub[deliveryHubs[i]]=1;
    }
    for(int i=0;i<deliveryHouses.size();i++){
        map_ifhouse[deliveryHouses[i]]=1;
    }
    for(int i=0;i<deliveryHouses.size();i++){
        map_house_hub[deliveryHouses[i]]=deliveryHubs[i];
    }

if(numNodes<=8){
        int start;
        for (int i = 0; i < numDeliveries; ++i) {
            if(minDist[deliveryHubs[i]].first>fuelCapacity/2) start=deliveryHubs[i];
        }
        vector<bool>visited(numNodes,0);
        
        int sum=0;

        dfs4(start,visited,minimumSpanningTree,fuelStations,sum);
        start=vec[1];

        bool temp=dfs(fuelStations[0],visited, minimumSpanningTree,sum,fuelStations[0],minDist, fuelCapacity);
        visited[start]=0;
        dfs4(start,visited,minimumSpanningTree,fuelStations,sum);
        cout<<vec.size()<<endl;
        for(auto e:vec){
            cout<<e<<" ";
        }
        cout<<endl;
        return 0;
    }
    for(int i = 0; i < numNodes; i++)
    {
        if (s.count(i) == 0 && minDist[i].first > fuelCapacity/2)
        {
            visited[i] = true;
        }
    }
    int sum=0;

    for(int i = 0; i <fuelStations.size(); i++)
    {
        parent[fuelStations[i]] = -1; 
    }

    vec.clear();

    for(int i = 0; i < vec2.size(); i++)
    {
        vector<int>subtree(numNodes,1);
        vector<bool> visitedDum(numNodes , false);

        for(int j = 0; j < numNodes; j++)
        {
            if (houseAndHubs.count(j) == 0)
            {
                subtree[j] = 0;
            }
        }

        dfsCheck(fuelStations[vec2[i]], visitedDum,visited,minimumSpanningTree,fuelStations[vec2[i]],minDist,houseAndHubs, subtree);

        while(1){
            int capacity_forone = fuelCapacity;
            bool temp=dfs(fuelStations[vec2[i]],visited, minimumSpanningTree,sum,fuelStations[vec2[i]],minDist, capacity_forone);
            if(visited[fuelStations[vec2[i]]]==1){ 
                break;
            }
        }
        for(auto e:vec){
            path.push_back(e);
        }
        vec.clear();
        if(i+1==vec2.size())break;
        //dijskta:fuelstation[vec2[i]]--fuelstation[vec2[i+1]],both exclusice 
        int start=fuelStations[vec2[i]];
        int end=fuelStations[vec2[i+1]];

        vector<int>visited2(numNodes,0);
        vector<int>parent2(numNodes,0);
        dfs3(start,minimumSpanningTree,visited2,parent2);

        vector<int>temp;
        int current=end;
        while(current!=start){
            if(current!=end){
                temp.push_back(current);
            }
            current=parent2[current];
        }
        reverse(temp.begin(),temp.end());
        for(auto e:temp){
            path.push_back(e);
        }

    }

/*    for(int i=0;i<numNodes;i++){
        if(map_ifhouse[i]==1 and visited_houses[i]==0){
            sum+=minDist[i].first*2;
        }
    }*/

    //cout<<"sum:"<<2*sum<<endl;

////////////////////////////////////////////////////////
    cout<<(path.size()*2) -1<<endl;
    for(auto e:path){
        cout<<e<<" ";
    }
    for(int i=1;i<path.size();i++){
        cout<<path[i]<<" ";
    }
    cout<<endl;


    map<pair<int,int>,int>mp;
    for(int i=0;i<minimumSpanningTree.size();i++){
        vector<pair<int,int>>v=minimumSpanningTree[i];
        for(auto pr:v){
            mp[{i,pr.first}]=1;
        }
    }
    int n;cin>>n;
    int a[n+1]={0};
    for(int i=1;i<=n;i++){
        cin>>a[i];
    }
    int cnt=0;
    for(int i=1;i<=n-1;i++){
        if(mp[{a[i],a[i+1]}]==0){
            cnt++;
        }
    }
    cout<<cnt<<endl;


    

    return 0;
}