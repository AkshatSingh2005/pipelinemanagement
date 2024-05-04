#include <iostream>
#include <vector>
#include <queue>
#include <climits>
#include <functional>
#include <algorithm>
#include<cmath>
#include <set>
#include <fstream>
#include <cstdlib>
#ifdef _WIN32
#define CLEAR_SCREEN "cls"
#else
#define CLEAR_SCREEN "clear"
#endif

using namespace std;

// Structure to represent an edge
struct Edge {
    int to, capacity, flow, age; 
    float replacementScore;
    float probabilityBlocked; // Added probabilityBlocked
    Edge(int t, int c, int a) : to(t), capacity(c), flow(0), age(a) {
        
        replacementScore = (static_cast<float>(age)) * capacity;
        
        probabilityBlocked = calculateProbabilityBlocked();
    }

    float calculateProbabilityBlocked() const {
        
        float k = 0.1; 
        float x0 = 0; 
        float L = 1; 
        return L / (1 + exp(-k * (replacementScore - x0)));
    }

};

// Function to add an edge to the graph
void addEdge(vector<vector<Edge>>& graph, int u, int v, int capacity, int age) {
    graph[u].push_back(Edge(v, capacity, age));
    graph[v].push_back(Edge(u, 0, age)); 
}


bool compareEdges(const Edge& edge1, const Edge& edge2) {
    
    return edge1.replacementScore > edge2.replacementScore;
}

void pipeReplacement(vector<vector<Edge>>& graph , int n){
    for (int i = 0; i < n; ++i) {
        sort(graph[i].begin(), graph[i].end(), compareEdges);
    }
    cout << "Edges to be replaced:" << endl;
    for (int u = 0; u < n; ++u) {
        for (const Edge& e : graph[u]) {
            cout << "(" << u << ", " << e.to << ")" << " Replacement Score: " << e.replacementScore << endl;
        }
    }
}

// Function to adjust flow rates using the MFIP algorithm
void adjustFlowRates(vector<vector<Edge>>& graph, vector<vector<int>>& currentFlowRates) {

    int n = graph.size();
    for (int u = 0; u < n; ++u) {
        for (Edge& e : graph[u]) {
            currentFlowRates[u][e.to] = rand() % (e.capacity + 1); // Set flow rate randomly up to the capacity
        }
    }
}



// Function to identify bottleneck edges
void findBottleneckEdges(vector<vector<Edge>>& graph, int source, int sink, set<pair<int, int>>& bottlenecks) {
    for (Edge& e : graph[source]) {
        if (e.flow == e.capacity) {
            bottlenecks.insert({source, e.to});
        }
    }
    for (int u = 0; u < graph.size(); ++u) {
        for (Edge& e : graph[u]) {
            if (e.flow == e.capacity && u != source && e.to != sink) {
                for (Edge& backEdge : graph[e.to]) {
                    if (backEdge.to == u && backEdge.flow == backEdge.capacity) {
                        bottlenecks.insert({u, e.to});
                        break;
                    }
                }
            }
        }
    }
}

void increaseBottleneckCapacities(vector<vector<Edge>>& graph, const set<pair<int, int>>& bottlenecks) {
    for (const auto& edge : bottlenecks) {
        int u = edge.first;
        int v = edge.second;
        for (Edge& e : graph[u]) {
            if (e.to == v) {
                e.capacity += (e.capacity - e.flow);
                break;
            }
        }
    }
}

// Ford-Fulkerson algorithm with dynamic capacity adjustment and printing augmented path and increased capacity
int fordFulkersonDynamicCapacity(vector<vector<Edge>>& graph, int source, int sink) {
    int n = graph.size();
    int maxFlow = 0;

    while (true) {

        vector<int> parent(n, -1);
        queue<int> q;
        q.push(source);
        parent[source] = source;

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (const Edge& e : graph[u]) {
                if (parent[e.to] == -1 && e.capacity > e.flow) {
                    parent[e.to] = u;
                    q.push(e.to);
                }
            }
        }

        if (parent[sink] == -1)
            break;

        int bottleneck = INT_MAX;
        vector<int> path;
        for (int v = sink; v != source; v = parent[v]) {
            int u = parent[v];
            for (Edge& e : graph[u]) {
                if (e.to == v) {
                    bottleneck = min(bottleneck, e.capacity - e.flow);
                    path.push_back(v); 
                    break;
                }
            }
        }

        reverse(path.begin(), path.end());

        // Print augmented path
        cout << "Augmented path: ";
        for (int v : path) {
            cout << v << " -> ";
        }
        cout << sink << ", Flow: " << bottleneck << endl;

        for (size_t i = 0; i < path.size() - 1; ++i) {
            int u = path[i];
            int v = path[i + 1];
            for (Edge& e : graph[u]) {
                if (e.to == v) {
                    e.flow += bottleneck;
                    break;
                }
            }
        }

        maxFlow += bottleneck;


        for (int u : path) {
            for (Edge& e : graph[u]) {
                if (e.to == parent[u]) {
                    e.capacity += bottleneck;
                    // Print edge
                    cout << "Increased capacity of edge (" << u << ", " << e.to << ") to " << e.capacity << endl;
                    break;
                }
            }
        }
    }

    return maxFlow;

    
}


void generateJSON(const vector<vector<Edge>>& graph) {
    ofstream jsonFile("graph.json");
    if (!jsonFile.is_open()) {
        cerr << "Failed to open graph.json for writing." << endl;
        return;
    }

    jsonFile << "{\n";
    jsonFile << "  \"nodes\": [\n";
    for (int i = 0; i < graph.size(); ++i) {
        if (i > 0) {
            jsonFile << ",\n"; // Add comma before each item except for the first one
        }
        jsonFile << "    {\"id\": " << i << "}";
    }
    jsonFile << "\n";
    jsonFile << "  ],\n";
    jsonFile << "  \"links\": [\n";
    for (int u = 0; u < graph.size(); ++u) {
        for (size_t j = 0; j < graph[u].size(); ++j) {
            const Edge& e = graph[u][j];
            if (u > 0 || j > 0) {
                jsonFile << ",\n"; // Add comma before each item except for the first one
            }
            jsonFile << "    {\"source\": " << u << ", \"target\": " << e.to << ", \"capacity\": " << e.capacity << "}";
        }
    }
    jsonFile << "\n";
    jsonFile << "  ]\n";
    jsonFile << "}\n";

    jsonFile.close();
}


void generateJSONWithAugmentedPaths(const vector<vector<Edge>>& graph, const set<vector<int>>& augmentedPaths) {
    ofstream jsonFile("graph1.json");
    if (!jsonFile.is_open()) {
        cerr << "Failed to open graph1.json for writing." << endl;
        return;
    }

    jsonFile << "{\n";
    jsonFile << "  \"nodes\": [\n";
    for (int i = 0; i < graph.size(); ++i) {
        if (i > 0) {
            jsonFile << ",\n"; // Add comma before each item except for the first one
        }
        jsonFile << "    {\"id\": " << i << "}";
    }
    jsonFile << "\n";
    jsonFile << "  ],\n";
    jsonFile << "  \"links\": [\n";
    for (int u = 0; u < graph.size(); ++u) {
        for (size_t j = 0; j < graph[u].size(); ++j) {
            const Edge& e = graph[u][j];
            if (u > 0 || j > 0) {
                jsonFile << ",\n"; 
            }
            jsonFile << "    {\"source\": " << u << ", \"target\": " << e.to << ", \"capacity\": " << e.capacity;
            bool isAugmented = false;
            for (const vector<int>& path : augmentedPaths) {
                auto it = find(path.begin(), path.end(), u);
                if (it != path.end() && it + 1 != path.end() && *(it + 1) == e.to) {
                    isAugmented = true;
                    break;
                }
            }
            jsonFile << ", \"augmented\": " << (isAugmented ? "true" : "false");
            jsonFile << "}";
        }
    }
    jsonFile << "\n";
    jsonFile << "  ]\n";
    jsonFile << "}\n";

    jsonFile.close();
}

void generateJSONWithAugmentedPaths1(const vector<vector<Edge>>& graph, const set<vector<int>>& augmentedPaths) {
    ofstream jsonFile("graph2.json");
    if (!jsonFile.is_open()) {
        cerr << "Failed to open graph1.json for writing." << endl;
        return;
    }

    jsonFile << "{\n";
    jsonFile << "  \"nodes\": [\n";
    for (int i = 0; i < graph.size(); ++i) {
        if (i > 0) {
            jsonFile << ",\n";
        }
        jsonFile << "    {\"id\": " << i << "}";
    }
    jsonFile << "\n";
    jsonFile << "  ],\n";
    jsonFile << "  \"links\": [\n";
    for (int u = 0; u < graph.size(); ++u) {
        for (size_t j = 0; j < graph[u].size(); ++j) {
            const Edge& e = graph[u][j];
            if (u > 0 || j > 0) {
                jsonFile << ",\n"; // Add comma before each item except for the first one
            }
            jsonFile << "    {\"source\": " << u << ", \"target\": " << e.to << ", \"capacity\": " << e.capacity;
            bool isAugmented = false;

            for (const vector<int>& path : augmentedPaths) {
                auto it = find(path.begin(), path.end(), u);
                if (it != path.end() && it + 1 != path.end() && *(it + 1) == e.to) {
                    isAugmented = true;
                    break;
                }
            }
            jsonFile << ", \"augmented\": " << (isAugmented ? "true" : "false");
            jsonFile << "}";
        }
    }
    jsonFile << "\n";
    jsonFile << "  ]\n";
    jsonFile << "}\n";

    jsonFile.close();
}





// Function to print the residual graph
void printResidualGraph(vector<vector<Edge>>& graph) {
    cout << "Residual Graph:" << endl;
    for (int u = 0; u < graph.size(); ++u) {
        cout << "Node " << u << ": ";
        for (const Edge& e : graph[u]) {
            cout << "(" << e.to << ", " << e.capacity - e.flow << ") ";
        }
        cout << endl;
    }
    cout << endl;
}

// Ford-Fulkerson algorithm with printing of residual graph and tracking augmented paths
int fordFulkerson(vector<vector<Edge>> graph, int source, int sink, set<vector<int>>& augmentedPaths) {
    int n = graph.size();
    vector<bool> visited(n, false);
    int maxFlow = 0;

    // DFS function for finding augmenting path
    function<int(int, int, vector<int>&)> dfs = [&](int u, int minCapacity, vector<int>& path) {
        if (u == sink) {
            path.push_back(u);
            return minCapacity;
        }
        visited[u] = true;
        for (Edge& e : graph[u]) {
            if (!visited[e.to] && e.capacity > e.flow) {
                int newMinCapacity = min(minCapacity, e.capacity - e.flow);
                int flow = dfs(e.to, newMinCapacity, path);
                if (flow > 0) {
                    e.flow += flow;
                    for (Edge& backEdge : graph[e.to]) {
                        if (backEdge.to == u) {
                            backEdge.flow -= flow;
                            break;
                        }
                    }
                    path.push_back(u);
                    return flow;
                }
            }
        }
        return 0;
    };

    // Repeat DFS until no augmenting path exists
    vector<int> path;
    while (true) {
        fill(visited.begin(), visited.end(), false);
        path.clear();
        int flow = dfs(source, INT_MAX, path);
        if (flow == 0) break;
        maxFlow += flow;

        // Store augmented path
        augmentedPaths.insert(path);

        // Print augmented path
        cout << "Augmented path (Ford-Fulkerson): ";
        for (int i = path.size() - 1; i >= 0; --i) {
            cout << path[i];
            if (i > 0) cout << " -> ";
        }
        cout << ", Flow: " << flow << endl;

        // Print residual graph
        printResidualGraph(graph);
    }

    return maxFlow;
}


int edmondsKarp(vector<vector<Edge>>& graph, int source, int sink, set<vector<int>>& augmentedPaths) {
    int V = graph.size(); // Number of vertices
    vector<vector<int>> residualGraph(V, vector<int>(V, 0));
    for (int u = 0; u < V; ++u) {
        for (const Edge& e : graph[u]) {
            residualGraph[u][e.to] += e.capacity - e.flow;
        }
    }

    vector<int> parent(V);
    int maxFlow = 0;

    while (true) {
        fill(parent.begin(), parent.end(), -1);
        queue<int> q;
        q.push(source);
        parent[source] = source;

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            if (u == sink)
                break;

            for (int v = 0; v < V; ++v) {
                if (parent[v] == -1 && residualGraph[u][v] > 0) {
                    parent[v] = u;
                    q.push(v);
                }
            }
        }

        if (parent[sink] == -1)
            break;

        int pathFlow = INT_MAX;
        for (int v = sink; v != source; v = parent[v]) {
            int u = parent[v];
            pathFlow = min(pathFlow, residualGraph[u][v]);
        }

        vector<int> path;
        int curr = sink;
        while (curr != source) {
            path.push_back(curr);
            curr = parent[curr];
        }
        path.push_back(source);

        // Update augmented paths
        augmentedPaths.insert(path);

        // Print augmented path
        cout << "Augmented path (Edmonds-Karp): ";
        for (int i = path.size() - 1; i >= 0; --i) {
            cout << path[i];
            if (i > 0) cout << " -> ";
        }
        cout << ", Flow: " << pathFlow << endl;

        for (int v = sink; v != source; v = parent[v]) {
            int u = parent[v];
            residualGraph[u][v] -= pathFlow;
            residualGraph[v][u] += pathFlow;
        }

        maxFlow += pathFlow;

        // Print residual graph
        //printResidualGraph(residualGraph);
    }

    return maxFlow;
}


// Dinic's algorithm for maximum flow
int dinic(vector<vector<Edge>>& graph, int source, int sink) {
    int n = graph.size();
    vector<int> dist(n), ptr(n);
    int maxFlow = 0;

    while (true) {
        fill(dist.begin(), dist.end(), -1);
        queue<int> q;
        q.push(source);
        dist[source] = 0;

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (Edge& e : graph[u]) {
                if (dist[e.to] == -1 && e.capacity > e.flow) {
                    dist[e.to] = dist[u] + 1;
                    q.push(e.to);
                }
            }
        }

        if (dist[sink] == -1)
            break;

        fill(ptr.begin(), ptr.end(), 0);
        function<int(int, int)> dfs = [&](int u, int minCapacity) {
            if (u == sink)
                return minCapacity;
            for (int& i = ptr[u]; i < graph[u].size(); ++i) {
                Edge& e = graph[u][i];
                if (dist[e.to] == dist[u] + 1 && e.capacity > e.flow) {
                    int flow = dfs(e.to, min(minCapacity, e.capacity - e.flow));
                    if (flow > 0) {
                        e.flow += flow;
                        for (Edge& backEdge : graph[e.to]) {
                            if (backEdge.to == u) {
                                backEdge.flow -= flow;
                                break;
                            }
                        }
                        return flow;
                    }
                }
            }
            return 0;
        };

        int flow;
        while ((flow = dfs(source, INT_MAX)) > 0) {
            maxFlow += flow;
        }
    }

    return maxFlow;
}

void callFordFulkerson(vector<vector<Edge>>& graph, int source, int sink, set<vector<int>>& augmentedPaths){
    int maxFlow = fordFulkerson(graph, source, sink, augmentedPaths);

    set<vector<int>> reversedAugmentedPaths;

    // Iterate over each vector in augmentedPaths, reverse it, and insert into reversedAugmentedPaths
    for (const vector<int>& path : augmentedPaths) {
        vector<int> reversedPath = path;
        reverse(reversedPath.begin(), reversedPath.end());
        reversedAugmentedPaths.insert(reversedPath);
    }

    generateJSONWithAugmentedPaths(graph, reversedAugmentedPaths);

    cout << "Maximum flow: " << maxFlow << endl;

    for (const vector<int>& path : augmentedPaths) {
        cout << "Path: ";
        for (int node : path) {
            cout << node << " -> ";
        }
        cout << "Sink" << endl;
    }
}

void callEdmondKarp(vector<vector<Edge>>& graph, int source, int sink, set<vector<int>>& augmentedPaths){
    int maxFlow = edmondsKarp(graph, source, sink, augmentedPaths);

    set<vector<int>> reversedAugmentedPaths;

    // Iterate over each vector in augmentedPaths, reverse it, and insert into reversedAugmentedPaths
    for (const vector<int>& path : augmentedPaths) {
        vector<int> reversedPath = path;
        reverse(reversedPath.begin(), reversedPath.end());
        reversedAugmentedPaths.insert(reversedPath);
    }

    generateJSONWithAugmentedPaths1(graph, reversedAugmentedPaths);

    cout << "Maximum flow: " << maxFlow << endl;

    for (const vector<int>& path : augmentedPaths) {
        cout << "Path: ";
        for (int node : path) {
            cout << node << " -> ";
        }
        cout << "Sink" << endl;
    }
}

void printProbablity(vector<vector<Edge>> graph){
    for (int u = 0; u < graph.size(); ++u) {
        for (const auto& edge : graph[u]) {
            if(edge.probabilityBlocked != 0.5)
            cout << "Edge from " << u << " to " << edge.to << ", Probability of being blocked: " << edge.probabilityBlocked << endl;
        }
    }
}

int main() {

    // Display welcome screen
    cout << "Welcome to Graph Operations!" << endl;
    cout << "Press any key to continue..." << endl;
    cin.get(); // Wait for user to press any key
    system(CLEAR_SCREEN); // Clear terminal


    int n, m; // Number of vertices and edges respectively
    cout << "Enter the number of vertices and edges: ";
    cin >> n >> m;
    vector<vector<Edge>> graph(n);
    vector<vector<int>> currentFlowRates(n, vector<int>(n, 0));
    set<vector<int>> augmentedPaths; // Set to store augmented paths
    set<vector<int>> augmentedPaths1; 


    
    // Read graph information from user input
    for (int i = 0; i < m; ++i) {
        int u, v, capacity, age;
        
        cout << "Enter the source vertex for edge " << i+1 << ": ";
        cin >> u;
        cout << "Enter the destination vertex for edge " << i+1 << ": ";
        cin >> v;
        cout << "Enter the capacity for edge " << i+1 << ": ";
        cin >> capacity;
        cout << "Enter the age for edge " << i+1 << ": ";
        cin >> age;
        addEdge(graph, u, v, capacity, age);
    }

    int source, sink; // Source and sink vertices
    cout << "Enter the source and sink vertices: ";
    cin >> source >> sink;

    // Generate JSON representation of the graph
    generateJSON(graph);

    // Option selection menu
    int option;
    while (true) {
        cout << "Graph Operations Menu:" << endl;
        cout << "1. Find Max Flow" << endl;
        cout << "2. Find Pipes Needing Replacement" << endl;
        cout << "3. Find Apparent Blocked Pipes" << endl;
        cout << "4. MFIP" << endl;
        cout << "5. Exit" << endl;
        cout << "Enter your choice: ";
        cin >> option;

        system(CLEAR_SCREEN); 

        if (option == 1) {
            if (option == 1) {
                // Handle finding max flow
                int algorithmChoice;
                cout << "Choose algorithm for finding max flow:" << endl;
                cout << "1. Ford-Fulkerson" << endl;
                cout << "2. Edmonds-Karp" << endl;
                cout << "3. Dinic" << endl;
                cout << "Enter your choice: ";
                cin >> algorithmChoice;



                if (algorithmChoice == 1) {
                    // Call Ford-Fulkerson function
                    callFordFulkerson(graph, source, sink, augmentedPaths);
                } else if (algorithmChoice == 2) {
                    // Call Edmonds-Karp function
                    callEdmondKarp(graph, source, sink, augmentedPaths1);
                } else if (algorithmChoice == 3) {
                    // Call Dinic function
                    cout << "Maximum flow using Dinic's algorithm: " << dinic(graph, source, sink) << endl;
                } else {
                    cout << "Invalid algorithm choice. Please choose again." << endl;
                }
            }

        } else if (option == 2) {
            // Handle finding pipes needing replacement
            pipeReplacement(graph, n);
        } else if (option == 3) {
            // Handle finding apparent blocked pipes

            printProbablity(graph);
        } else if (option == 4) {
            // Handle MFIP

            cout << fordFulkersonDynamicCapacity(graph,source,sink);
        } else if (option == 5) {
            // Exit the program
            cout << "Exiting..." << endl;
            break;
        } else {
            cout << "Invalid option. Please choose again." << endl;
        }

        // Clear the terminal
        system("pause"); // Pause for the user to press a key
        system("clear || cls");

    }





/*

    // Input edges and capacities
    cout << "Enter the edges, their capacities, and age (from_vertex to_vertex capacity age):" << endl;
    for (int i = 0; i < m; ++i) {
        int u, v, capacity, age;
        
        cout << "Enter the source vertex for edge " << i+1 << ": ";
        cin >> u;
        cout << "Enter the destination vertex for edge " << i+1 << ": ";
        cin >> v;
        cout << "Enter the capacity for edge " << i+1 << ": ";
        cin >> capacity;
        //cout << "Enter the age for edge " << i+1 << ": ";
        //cin >> age;
        addEdge(graph, u, v, capacity, age);
    }

    int source, sink; // Source and sink vertices
    cout << "Enter the source and sink vertices: ";
    cin >> source >> sink;

    // Generate JSON representation of the graph
    //generateJSON(graph);
    
    callEdmondKarp(graph , source , sink , augmentedPaths1);

    

    

    // Calculate maximum flow using Ford-Fulkerson algorithm
    //cout << "Maximum flow using Ford-Fulkerson algorithm: " << fordFulkerson(graph, source, sink) << endl;


    // Calculate maximum flow using Edmonds-Karp algorithm
    //cout << "Maximum flow using Edmonds-Karp algorithm: " << edmondsKarp(graph, source, sink) << endl;

    // Calculate maximum flow using Dinic's algorithm
    //cout << "Maximum flow using Dinic's algorithm: " << dinic(graph, source, sink) << endl;

    // Print edges to be replaced
    //pipeReplacement(graph, n);
*/
    return 0;

}
