#include <bits/stdc++.h>
using namespace std;
using pii = pair<long long , long long >;

struct Edge {
    int from ,to, cost;
    bool operator<(const Edge& other) const {
        return cost < other.cost;
    }
};
vector<int>vec,vec2,path;
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
    priority_queue<pii, vector<pii>, greater<pii>> pq;

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


bool dfs(int vertex,vector<bool>&visited,vector<vector<pair<int , int >>> &minimumSpanningTree,
    vector<int>&parent,int &ans,int fuelstation,vector<pair<int,int>>&mindist,int &currentfuel,map<int,int>&map_ifhub,vector<int>&visited_hubs,map<int,int>&map_ifhouse,vector<int>&visited_houses,map<int,int>&map_house_hub){
    visited[vertex]=1;
    vec.push_back(vertex);
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
            //
            done=0;
        }
        else{
            parent[child] = vertex;
            ans+=2*wt;
            currentfuel -= wt;
            done=((done)&dfs(child,visited, minimumSpanningTree,parent,ans,fuelstation,mindist,(currentfuel),map_ifhub,visited_hubs,map_ifhouse,visited_houses,map_house_hub));
            vec.push_back(vertex);
        }
    }
    if(done==0){
        visited[vertex]=0;
    }
    return done;
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

void dfs2(int vertex,vector<bool>&visited,vector<vector<pair<int , int >>> &minimumSpanningTree,vector<int>&path,vector<int>&fuelStations,int &sum){
    visited[vertex]=1;
    path.push_back(vertex);
    for(auto pr:minimumSpanningTree[vertex]){
        if(visited[pr.first]==1)continue;
        if(pr.first==fuelStations[0]){
            sum+=pr.second*2;
            return;
        }

        dfs2(pr.first,visited,minimumSpanningTree,path,fuelStations,sum);
        sum+=pr.second;

    }
    return;
 }


int  main() {  

    int  numDeliveries, numNodes, numRoads, numFuelStations, fuelCapacity;
    cin >> numDeliveries >> numNodes >> numRoads >> numFuelStations >> fuelCapacity;
    vector<int > deliveryHubs(numDeliveries);
    set<int> s, houseAndHubs;

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
        edges.push_back({from, to, cost});
    }
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
    map<int,int>map_ifhub;
    map<int,int>map_ifhouse;
    map<int,int>map_house_hub;

    for(int i=0;i<deliveryHubs.size();i++){
        map_ifhub[deliveryHubs[i]]=1;
    }
    for(int i=0;i<deliveryHouses.size();i++){
        map_ifhouse[deliveryHouses[i]]=1;
    }
    for(int i=0;i<deliveryHouses.size();i++){
        map_house_hub[deliveryHouses[i]]=deliveryHubs[i];
    }
    vector<int>visited_hubs(numNodes,0);
    vector<int>visited_houses(numNodes,0);
    vector<int>parent(numNodes,0);
    if(numNodes<=8){
            int start;
            for (int i = 0; i < numDeliveries; ++i) {
                if(minDist[deliveryHubs[i]].first>fuelCapacity/2) start=deliveryHubs[i];
            }
            vector<bool>visited(numNodes,0);
            vector<int>path;
            int sum=0;

            dfs2(start,visited,minimumSpanningTree,path,fuelStations,sum);
            start=path[1];
            bool temp=dfs(fuelStations[0],visited, minimumSpanningTree,parent,sum,fuelStations[0],minDist, fuelCapacity,map_ifhub,visited_hubs,map_ifhouse,visited_houses,map_house_hub);
            visited[start]=0;
            for(auto e:vec){
               path.push_back(e);
            }
            dfs2(start,visited,minimumSpanningTree,path,fuelStations,sum);
            cout<<path.size()<<endl;
            for(auto e:path){
                cout<<e<<" ";
            }
            cout<<endl;
            return 0;
        }

    int sum=0;
    for(int i = 0; i <fuelStations.size(); i++)
    {
        parent[fuelStations[i]] = -1; 
    }

    vector<int>petrol=fuelStations;//same vector hai


    for(int i = 0; i < petrol.size(); i++)
    {
        vector<bool> visited(numNodes,false);


        vector<int>subtree(numNodes,1);
        vector<bool> visitedDum(numNodes , false);

        for(int j = 0; j < numNodes; j++)
        {
            if (houseAndHubs.count(j) == 0)
            {
                subtree[j] = 0;
            }
        }

        dfsCheck(fuelStations[i], visitedDum,visited,minimumSpanningTree,fuelStations[i],minDist,houseAndHubs, subtree);

        while(1){

            for(int k = 0; k < numNodes; k++)
            {
                if (minDist[k].first > fuelCapacity/2)
                {
                    visited[k] = true;
                }
            }

            int cp=fuelCapacity;
            bool tep=dfs(petrol[i],visited, minimumSpanningTree,parent,sum,petrol[i],minDist, cp,map_ifhub,visited_hubs,map_ifhouse,visited_houses,map_house_hub);
            if(visited[petrol[i]]==1){
                break;
            }
        }
            for(auto e:vec){
                path.push_back(e);
            }
            vec.clear();
            if(i+1==petrol.size())break;

            int start=petrol[i];
            int end=petrol[i+1];

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
    cout<<(2*path.size())-1<<endl;
    for(int i=0;i<path.size();i++){
        cout<<path[i]<<" ";
    }
    path.pop_back();
    reverse(path.begin(),path.end());
    for(int i=0;i<path.size();i++){
        cout<<path[i]<<" ";
    }


    
    return 0;
}
