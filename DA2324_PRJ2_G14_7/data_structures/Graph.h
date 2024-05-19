// Original code by Gonçalo Leão
// Updated by DA 2023/2024 Team

#ifndef FEUP_DA_2024_PRJ2_G07_T14_GRAPH_H
#define FEUP_DA_2024_PRJ2_G07_T14_GRAPH_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include "data_structures/MutablePriorityQueue.h"
#include "../Class/Edges.h"


using namespace std;

template <class T>
class Edge;

#define INF std::numeric_limits<double>::max()

/************************* Vertex  **************************/

template <class T>
class Vertex {
public:
    Vertex(T in);
    bool operator<(Vertex<T> & vertex) const; // // required by MutablePriorityQueue

    T getInfo() const;
    std::vector<Edge<T> *> getAdj() const;
    bool isVisited() const;
    bool isProcessing() const;
    unsigned int getIndegree() const;
    double getDist() const;
    Edge<T> *getPath() const;
    std::vector<Edge<T> *> getIncoming() const;
    void addChild(string s);
    void eraseChild();
    std::vector<string> getChild();

    void setInfo(T info);
    void setVisited(bool visited);
    void setProcesssing(bool processing);
    void setIndegree(unsigned int indegree);
    void setDist(double dist);
    void setPath(Edge<T> *path);
    Edge<T> * addEdge(Vertex<T> *dest, double w);
    bool removeEdge(T in);
    void removeOutgoingEdges();

    double getFlow() const;

    void setFlow(double flow);

    friend class MutablePriorityQueue<Vertex>;
protected:
    T info;                // info node
    std::vector<Edge<T> *> adj;  // outgoing edges
    std::vector<string> child;

    double flow = 0; // for flow-related problems

    // auxiliary fields
    bool visited = false; // used by DFS, BFS, Prim ...
    bool processing = false; // used by isDAG (in addition to the visited attribute)
    unsigned int indegree; // used by topsort
    double dist = 0;
    Edge<T> *path = nullptr;

    std::vector<Edge<T> *> incoming; // incoming edges

    int queueIndex = 0; 		// required by MutablePriorityQueue and UFDS

    void deleteEdge(Edge<T> *edge);
};

/********************** Edge  ****************************/

template <class T>
class Edge {
public:
    Edge(Vertex<T> *orig, Vertex<T> *dest, double w);

    Vertex<T> * getDest() const;
    double getWeight() const;
    bool isSelected() const;
    Vertex<T> * getOrig() const;
    Edge<T> *getReverse() const;
    double getFlow() const;

    void setSelected(bool selected);
    void setReverse(Edge<T> *reverse);
    void setFlow(double flow);
protected:
    Vertex<T> *dest; // destination vertex
    double weight; // edge weight, can also be used for capacity

    // auxiliary fields
    bool selected = false;

    // used for bidirectional edges
    Vertex<T> *orig;
    Edge<T> *reverse = nullptr;

    double flow; // for flow-related problems
};

/********************** Graph  ****************************/

template <class T>
class Graph {
public:
    //Graph(const unordered_map<string, Vertex<T> *> &vertexSet);

    ~Graph();
    /*
    * Auxiliary function to find a vertex with a given the content.
    */
    Vertex<T> *findVertex(const T &in) const;
    /*
     *  Adds a vertex with a given content or info (in) to a graph (this).
     *  Returns true if successful, and false if a vertex with that content already exists.
     */
    bool addVertex(const T &in);
    bool removeVertex(const T &in);

    //void MaxFlow(Graph<string> *g, const unordered_map<string, Reservoirs *> *wr, const unordered_map<string, Cities *> *ds);

    /*
     * Adds an edge to a graph (this), given the contents of the source and
     * destination vertices and the edge weight (w).
     * Returns true if successful, and false if the source or destination vertex does not exist.
     */
    bool addEdge(const T &sourc, const T &dest, double w);
    bool removeEdge(const T &source, const T &dest);
    bool addBidirectionalEdge(const T &sourc, const T &dest, double w);

    int getNumVertex() const;
    unordered_map<string, Vertex<T> *> getVertexSet() const;

    std:: vector<T> dfs() const;
    std:: vector<T> dfs(const T & source) const;
    void dfsVisit(Vertex<T> *v,  std::vector<T> & res) const;
    std::vector<T> bfs(const T & source) const;

    bool isDAG() const;
    bool dfsIsDAG(Vertex<T> *v) const;
    std::vector<T> topsort() const;

//    void Backtrack(Graph<T> *g, const unordered_map<string, Edges *> *edges);
//    void BTUtil(Graph<T> *g, vector<string> &path, double dist, vector<string> &shortestPath, double &shortestDist, Vertex<T> *checkVertex);
    vector<Vertex<T>*> nearestNeighborTSP();

//    void Backtracking(Graph<T> *g, const unordered_map<string, Edges *> *edges);
//    void TriangleApproximation(Graph<T> *g, const unordered_map<string, Edges *> *edges);
//    void OtherHeuristic(Graph<T> *g, const unordered_map<string, Edges *> *edges);
protected:
    unordered_map<string, Vertex<T> *> vertexSet;    // vertex set

    double ** distMatrix = nullptr;   // dist matrix for Floyd-Warshall
    int **pathMatrix = nullptr;   // path matrix for Floyd-Warshall

    /*
     * Finds the index of the vertex with a given content.
     */
    int findVertexIdx(const T &in) const;


};

void deleteMatrix(int **m, int n);
void deleteMatrix(double **m, int n);


/************************* Vertex  **************************/

template <class T>
Vertex<T>::Vertex(T in): info(in) {}
/*
 * Auxiliary function to add an outgoing edge to a vertex (this),
 * with a given destination vertex (d) and edge weight (w).
 */
template <class T>
Edge<T> * Vertex<T>::addEdge(Vertex<T> *d, double w) {
    auto newEdge = new Edge<T>(this, d, w);
    adj.push_back(newEdge);
    d->incoming.push_back(newEdge);
    return newEdge;
}

/*
 * Auxiliary function to remove an outgoing edge (with a given destination (d))
 * from a vertex (this).
 * Returns true if successful, and false if such edge does not exist.
 */
template <class T>
bool Vertex<T>::removeEdge(T in) {
    bool removedEdge = false;
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge<T> *edge = *it;
        Vertex<T> *dest = edge->getDest();
        if (dest->getInfo() == in) {
            it = adj.erase(it);
            deleteEdge(edge);
            removedEdge = true; // allows for multiple edges to connect the same pair of vertices (multigraph)
        }
        else {
            it++;
        }
    }
    return removedEdge;
}

/*
 * Auxiliary function to remove an outgoing edge of a vertex.
 */
template <class T>
void Vertex<T>::removeOutgoingEdges() {
    auto it = adj.begin();
    while (it != adj.end()) {
        Edge<T> *edge = *it;
        it = adj.erase(it);
        deleteEdge(edge);
    }
}

template <class T>
bool Vertex<T>::operator<(Vertex<T> & vertex) const {
    return this->dist < vertex.dist;
}

template <class T>
T Vertex<T>::getInfo() const {
    return this->info;
}

template <class T>
vector<Edge<T> *> Vertex<T>::getAdj() const {
    return this->adj;
}

template <class T>
bool Vertex<T>::isVisited() const {
    return this->visited;
}

template <class T>
bool Vertex<T>::isProcessing() const {
    return this->processing;
}

template <class T>
unsigned int Vertex<T>::getIndegree() const {
    return this->indegree;
}

template <class T>
double Vertex<T>::getDist() const {
    return this->dist;
}

template <class T>
Edge<T> *Vertex<T>::getPath() const {
    return this->path;
}

template <class T>
std::vector<Edge<T> *> Vertex<T>::getIncoming() const {
    return this->incoming;
}

template <class T>
void Vertex<T>::setInfo(T in) {
    this->info = in;
}

template <class T>
void Vertex<T>::setVisited(bool visited) {
    this->visited = visited;
}

template <class T>
void Vertex<T>::setProcesssing(bool processing) {
    this->processing = processing;
}

template <class T>
void Vertex<T>::setIndegree(unsigned int indegree) {
    this->indegree = indegree;
}

template <class T>
void Vertex<T>::addChild(std::string s) {
    this->child.push_back(s);
}

template <class T>
void Vertex<T>::eraseChild() {
    child.clear();
}

template <class T>
std::vector<string> Vertex<T>::getChild(){
    return this->child;
}

template <class T>
void Vertex<T>::setDist(double dist) {
    this->dist = dist;
}

template <class T>
void Vertex<T>::setPath(Edge<T> *path) {
    this->path = path;
}

template <class T>
void Vertex<T>::deleteEdge(Edge<T> *edge) {
    Vertex<T> *dest = edge->getDest();
    // Remove the corresponding edge from the incoming list
    auto it = dest->incoming.begin();
    while (it != dest->incoming.end()) {
        if ((*it)->getOrig()->getInfo() == info) {
            it = dest->incoming.erase(it);
        }
        else {
            it++;
        }
    }
    delete edge;
}

template<class T>
double Vertex<T>::getFlow() const {
    return flow;
}

template<class T>
void Vertex<T>::setFlow(double flow) {
    this->flow = flow;
}

/********************** Edge  ****************************/

template <class T>
Edge<T>::Edge(Vertex<T> *orig, Vertex<T> *dest, double w): orig(orig), dest(dest), weight(w) {}

template <class T>
Vertex<T> * Edge<T>::getDest() const {
    return this->dest;
}

template <class T>
double Edge<T>::getWeight() const {
    return this->weight;
}

template <class T>
Vertex<T> * Edge<T>::getOrig() const {
    return this->orig;
}

template <class T>
Edge<T> *Edge<T>::getReverse() const {
    return this->reverse;
}

template <class T>
bool Edge<T>::isSelected() const {
    return this->selected;
}

template <class T>
double Edge<T>::getFlow() const {
    return flow;
}

template <class T>
void Edge<T>::setSelected(bool selected) {
    this->selected = selected;
}

template <class T>
void Edge<T>::setReverse(Edge<T> *reverse) {
    this->reverse = reverse;
}

template <class T>
void Edge<T>::setFlow(double flow) {
    this->flow = flow;
}

/********************** Graph  ****************************/

template <class T>
int Graph<T>::getNumVertex() const {
    return vertexSet.size();
}

template <class T>
unordered_map<string, Vertex<T> *> Graph<T>::getVertexSet() const {
    return this->vertexSet;
}

/*
 * Auxiliary function to find a vertex with a given content.
 */
template <class T>
Vertex<T> * Graph<T>::findVertex(const T &in) const {
    for (auto &pair : vertexSet)
        if (pair.second->getInfo() == in)
            return pair.second;
    return nullptr;
}

/*
 * Finds the index of the vertex with a given content.
 */
template <class T>
int Graph<T>::findVertexIdx(const T &in) const {
    for (unsigned i = 0; i < vertexSet.size(); i++)
        if (vertexSet[i]->getInfo() == in)
            return i;
    return -1;
}
/*
 *  Adds a vertex with a given content or info (in) to a graph (this).
 *  Returns true if successful, and false if a vertex with that content already exists.
 */
template <class T>
bool Graph<T>::addVertex(const T &in) {
    if (findVertex(in) != nullptr)
        return false;
    Vertex<T> *vertexnew = new Vertex<T>(in);
    vertexSet.insert({in, vertexnew});
    return true;
}

/*
 *  Removes a vertex with a given content (in) from a graph (this), and
 *  all outgoing and incoming edges.
 *  Returns true if successful, and false if such vertex does not exist.
 */
template <class T>
bool Graph<T>::removeVertex(const T &in) {
//    for (auto it = vertexSet.begin(); it != vertexSet.end(); it++) {
////        for(auto &pair : vertexSet){
////            it = pair.second;
//        if ((*it)->getInfo() == in) {
//            auto v = *it;
//            v->removeOutgoingEdges();
//            for (auto u : vertexSet) {
//                u->removeEdge(v->getInfo());
//            }
//            vertexSet.erase(it);
//            delete v;
//            return true;
//        }
//    }
    for(auto &pair : vertexSet) {
        if(pair.second->getInfo() == in) {
            auto v = pair.second;
            v->removeOutgoingEdges();
            for (auto &pair : vertexSet) {
                pair.second->removeEdge(v->getInfo());
            }
            vertexSet.erase(v->getInfo());
            delete v;
            return true;
        }
    }

    return false;
}

/*
 * Adds an edge to a graph (this), given the contents of the source and
 * destination vertices and the edge weight (w).
 * Returns true if successful, and false if the source or destination vertex does not exist.
 */
template <class T>
bool Graph<T>::addEdge(const T &sourc, const T &dest, double w) {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr)
        return false;
    v1->addEdge(v2, w);
    return true;
}

/*
 * Removes an edge from a graph (this).
 * The edge is identified by the source (sourc) and destination (dest) contents.
 * Returns true if successful, and false if such edge does not exist.
 */
template <class T>
bool Graph<T>::removeEdge(const T &sourc, const T &dest) {
    Vertex<T> * srcVertex = findVertex(sourc);
    if (srcVertex == nullptr) {
        return false;
    }
    return srcVertex->removeEdge(dest);
}

template <class T>
bool Graph<T>::addBidirectionalEdge(const T &sourc, const T &dest, double w) {
    auto v1 = findVertex(sourc);
    auto v2 = findVertex(dest);
    if (v1 == nullptr || v2 == nullptr)
        return false;
    auto e1 = v1->addEdge(v2, w);
    auto e2 = v2->addEdge(v1, w);
    e1->setReverse(e2);
    e2->setReverse(e1);
    return true;
}

/****************** DFS ********************/

/*
 * Performs a depth-first search (dfs) traversal in a graph (this).
 * Returns a vector with the contents of the vertices by dfs order.
 */
template <class T>
std::vector<T> Graph<T>::dfs() const {
    std::vector<T> res;
    for (auto v : vertexSet)
        v->setVisited(false);
    for (auto v : vertexSet)
        if (!v->isVisited())
            dfsVisit(v, res);
    return res;
}

/*
 * Performs a depth-first search (dfs) in a graph (this) from the source node.
 * Returns a vector with the contents of the vertices by dfs order.
 */
template <class T>
std::vector<T> Graph<T>::dfs(const T & source) const {
    std::vector<int> res;
    // Get the source vertex
    auto s = findVertex(source);
    if (s == nullptr) {
        return res;
    }
    // Set that no vertex has been visited yet
    for (auto v : vertexSet) {
        v->setVisited(false);
    }
    // Perform the actual DFS using recursion
    dfsVisit(s, res);

    return res;
}

/*
 * Auxiliary function that visits a vertex (v) and its adjacent, recursively.
 * Updates a parameter with the list of visited node contents.
 */
template <class T>
void Graph<T>::dfsVisit(Vertex<T> *v, std::vector<T> & res) const {
    v->setVisited(true);
    res.push_back(v->getInfo());
    for (auto & e : v->getAdj()) {
        auto w = e->getDest();
        if (!w->isVisited()) {
            dfsVisit(w, res);
        }
    }
}

/****************** BFS ********************/
/*
 * Performs a breadth-first search (bfs) in a graph (this), starting
 * from the vertex with the given source contents (source).
 * Returns a vector with the contents of the vertices by bfs order.
 */
template <class T>
std::vector<T> Graph<T>::bfs(const T & source) const {
    std::vector<int> res;
    // Get the source vertex
    auto s = findVertex(source);
    if (s == nullptr) {
        return res;
    }

    // Set that no vertex has been visited yet
    for (auto v : vertexSet) {
        v->setVisited(false);
    }

    // Perform the actual BFS using a queue
    std::queue<Vertex<T> *> q;
    q.push(s);
    s->setVisited(true);
    while (!q.empty()) {
        auto v = q.front();
        q.pop();
        res.push_back(v->getInfo());
        for (auto & e : v->getAdj()) {
            auto w = e->getDest();
            if ( ! w->isVisited()) {
                q.push(w);
                w->setVisited(true);
            }
        }
    }
    return res;
}

/****************** isDAG  ********************/
/*
 * Performs a depth-first search in a graph (this), to determine if the graph
 * is acyclic (acyclic directed graph or DAG).
 * During the search, a cycle is found if an edge connects to a vertex
 * that is being processed in the stack of recursive calls (see theoretical classes).
 * Returns true if the graph is acyclic, and false otherwise.
 */

template <class T>
bool Graph<T>::isDAG() const {
    for (auto v : vertexSet) {
        v->setVisited(false);
        v->setProcesssing(false);
    }
    for (auto v : vertexSet) {
        if (! v->isVisited()) {
            if ( ! dfsIsDAG(v) ) return false;
        }
    }
    return true;
}

/**
 * Auxiliary function that visits a vertex (v) and its adjacent, recursively.
 * Returns false (not acyclic) if an edge to a vertex in the stack is found.
 */
template <class T>
bool Graph<T>::dfsIsDAG(Vertex<T> *v) const {
    v->setVisited(true);
    v->setProcesssing(true);
    for (auto e : v->getAdj()) {
        auto w = e->getDest();
        if (w->isProcessing()) return false;
        if (! w->isVisited()) {
            if (! dfsIsDAG(w)) return false;
        }
    }
    v->setProcesssing(false);
    return true;
}

/****************** toposort ********************/
//=============================================================================
// Exercise 1: Topological Sorting
//=============================================================================
//
/*
 * Performs a topological sorting of the vertices of a graph (this).
 * Returns a vector with the contents of the vertices by topological order.
 * If the graph has cycles, returns an empty vector.
 * Follows the algorithm described in theoretical classes.
 */

template<class T>
std::vector<T> Graph<T>::topsort() const {
    std::vector<int> res;

    for (auto v : vertexSet) {
        v->setIndegree(0);
    }
    for (auto v : vertexSet) {
        for (auto e : v->getAdj()) {
            unsigned int indegree = e->getDest()->getIndegree();
            e->getDest()->setIndegree(indegree + 1);
        }
    }

    std::queue<Vertex<T> *> q;
    for (auto v : vertexSet) {
        if (v->getIndegree() == 0) {
            q.push(v);
        }
    }

    while( !q.empty() ) {
        Vertex<T> * v = q.front();
        q.pop();
        res.push_back(v->getInfo());
        for(auto e : v->getAdj()) {
            auto w = e->getDest();
            w->setIndegree(w->getIndegree() - 1);
            if(w->getIndegree() == 0) {
                q.push(w);
            }
        }
    }

    if ( res.size() != vertexSet.size() ) {
        //std::cout << "Impossible topological ordering!" << std::endl;
        res.clear();
        return res;
    }

    return res;
}

inline void deleteMatrix(int **m, int n) {
    if (m != nullptr) {
        for (int i = 0; i < n; i++)
            if (m[i] != nullptr)
                delete [] m[i];
        delete [] m;
    }
}

inline void deleteMatrix(double **m, int n) {
    if (m != nullptr) {
        for (int i = 0; i < n; i++)
            if (m[i] != nullptr)
                delete [] m[i];
        delete [] m;
    }
}

template <class T>
Graph<T>::~Graph() {
    deleteMatrix(distMatrix, vertexSet.size());
    deleteMatrix(pathMatrix, vertexSet.size());
}


template <class T>
void BTUtil(Graph<T> *g, vector<string> &path, double dist, vector<string> &shortestPath, double &shortestDist, Vertex<T> *checkVertex, int &numPathsTested, int &numRecursiveCalls , int &numBackTracks) {
    //cout << "teste" << endl;
    numRecursiveCalls++;

    if(path.size() == g->getVertexSet().size()) {
        for(auto & e : checkVertex->getAdj()) {
            auto w = e->getDest();
            if(w == g->findVertex("0")) {
                //cout << "teste2" << endl;
                numPathsTested++;
                dist += e->getWeight();
                path.push_back(w->getInfo());
                if(dist < shortestDist) {
                    shortestPath = path;
                    shortestDist = dist;
                }
                dist -= e->getWeight();
                path.pop_back();
            }
        }
        numBackTracks++;
        return;
    }

    for(auto & e : checkVertex->getAdj()) {
        auto w = e->getDest();
        if (!(w->isVisited())) {
            //cout << "testeRecursion" << endl;
            path.push_back(w->getInfo());
            dist += e->getWeight();
            w->setVisited(true);

            BTUtil(g, path, dist, shortestPath, shortestDist, w, numPathsTested, numRecursiveCalls ,numBackTracks);

            path.pop_back();
            dist -= e->getWeight();
            w->setVisited(false);

            numBackTracks++;
        }
    }
}

template<class T>
void Backtracking(Graph<T> *g, const unordered_map<string, Edges *> *edges) {
    vector<string> path;
    vector<string> shortestPath;
    double dist = 0;
    double shortestDist = INF;
    int numPathsTested = 0;
    int numRecursiveCalls = 0;
    int numBackTracks = 0;
    string initial = "0";
    Vertex<T>* initialVertex = g->findVertex("0");
//    path.push_back(initial);

    for (auto v : g->getVertexSet()) {
        if(v.second->getInfo() != initial)
            v.second->setVisited(false);
//        else v.second->setVisited(true);
    }

    auto start = chrono::high_resolution_clock::now();

    if (initialVertex) {
        initialVertex->setVisited(true);
        path.push_back(initialVertex->getInfo());
        BTUtil(g, path, dist, shortestPath, shortestDist, initialVertex, numPathsTested, numRecursiveCalls ,numBackTracks);
        initialVertex->setVisited(false);
    }

//    BTUtil(g, path, dist, shortestPath, shortestDist, initialVertex, numPathsTested, numRecursiveCalls ,numBackTracks);

    auto end = chrono::high_resolution_clock::now();
    auto durationMicroseconds = chrono::duration_cast<chrono::microseconds>(end - start).count();
    auto durationMilliseconds = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    auto duration = chrono::duration_cast<chrono::seconds >(end - start).count();

    cout << "Shortest Path Found: ";
    for (unsigned int i = 0; i < shortestPath.size(); ++i) {
        cout << shortestPath.at(i);
        if (i != shortestPath.size()-1) {
            cout << " -> ";
        }
    }
    cout << endl;
    cout << "Total Distance Traveled in path: " << shortestDist << endl;
    cout << "Total Number of Paths Tested: " << numPathsTested << endl;
    cout << "Total Number of Recursive Calls made: " << numRecursiveCalls << endl;
    cout << "Total Number of Backtracks done: " << numBackTracks << endl;
    if(duration != 0)
        cout << "Execution time: " << duration << " seconds" << endl;
    else if(durationMilliseconds != 0)
        cout << "Execution time: " << durationMilliseconds << " milliseconds" << endl;
    else
        cout << "Execution time: " << durationMicroseconds << " microseconds" << endl;
}




template<class T>
void TriangleApproximation(Graph<T> *g, const unordered_map<string, Edges *> *edges) {

    unordered_map<string , Vertex<T> *> vertix = g->getVertexSet();
    for(auto e: vertix){
        e.second->setVisited(false);
        e.second->setDist(INF);
        e.second->eraseChild();
    }

    auto start = chrono::high_resolution_clock::now();
    Vertex<T> *r = g->findVertex("0");
    r->setDist(0);
    MutablePriorityQueue<Vertex<T>> pq;
    pq.insert(r);
    while(!pq.empty()){
        auto i = pq.extractMin();
        i->setVisited(true);
        if(i->getInfo()!= "0"){
            i->getPath()->getOrig()->addChild(i->getInfo());
        }
        for(auto & q: i->getAdj()){
            Vertex<T> *w = q->getDest();
            if (!w->isVisited()){
                auto dist2 = w->getDist();
                if(q->getWeight() < dist2){
                    w->setDist(q->getWeight());
                    w->setPath(q);
                    if(dist2 == INF){
                        pq.insert(w);
                    }else{
                        pq.decreaseKey(w);
                    }
                }
            }
        }
    }
    vector<string> result;

    stack<Vertex<T> *> stack;
    stack.push(r);

    while (!stack.empty()) {
        Vertex<T> *v = stack.top();
        stack.pop();

        for (const auto &child : v->getChild()) {
            stack.push(g->findVertex(child));
        }

        result.push_back(v->getInfo());
    }

    result.push_back("0");
    double totalDistance = 0;
    for (int l = 0; l < result.size() - 1; ++l) {
        auto j = g->findVertex(result[l]);
        for (auto & k: j->getAdj()){
            if(k->getDest()->getInfo()==result[l+1]){
                totalDistance += k->getWeight();
            }
        }
    }

    auto end = chrono::high_resolution_clock::now();
    auto durationMicroseconds = chrono::duration_cast<chrono::microseconds>(end - start).count();
    auto durationMilliseconds = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    auto duration = chrono::duration_cast<chrono::seconds >(end - start).count();
    
    cout << "Approximate Path Found: ";
    for (size_t i = 0; i < result.size(); ++i) {
        cout << result[i];
        if (i != result.size() - 1) {
            cout << " -> ";
        }
    }
    cout << endl;
    cout << "Total Distance Traveled in path: " << totalDistance << endl;
    if(duration != 0)
        cout << "Execution time: " << duration << " seconds" << endl;
    else if(durationMilliseconds != 0)
        cout << "Execution time: " << durationMilliseconds << " milliseconds" << endl;
    else
        cout << "Execution time: " << durationMicroseconds << " microseconds" << endl;
}


template <class T>
vector<Vertex<T>*> Graph<T>::nearestNeighborTSP() {
    vector<Vertex<T>*> tour;

    // Initialize a set to keep track of visited vertices
    unordered_set<Vertex<T>*> visited;

    // Choose an arbitrary starting vertex (in this case, the first vertex)
    if (vertexSet.empty()) {
        return tour; // Return an empty tour if the graph is empty
    }
    string initial = "0";
    Vertex<T>* currentVertex = findVertex("0");

    tour.push_back(currentVertex);
    visited.insert(currentVertex);

    // Repeat until all vertices are visited
    while (visited.size() < vertexSet.size()) {
        double minDistance = std::numeric_limits<double>::max();
        Vertex<T>* nearestNeighbor = nullptr;

        // Find the nearest unvisited neighbor of the current vertex
        for (Edge<T>* edge : currentVertex->getAdj()) {
            Vertex<T>* neighbor = edge->getDest();
            if (visited.find(neighbor) == visited.end() && edge->getWeight() < minDistance) {
                minDistance = edge->getWeight();
                nearestNeighbor = neighbor;
            }
        }

        if (nearestNeighbor) {
            // Mark the nearest neighbor as visited
            visited.insert(nearestNeighbor);
            // Add the nearest neighbor to the tour
            tour.push_back(nearestNeighbor);
            // Update the current vertex to be the nearest neighbor
            currentVertex = nearestNeighbor;
        } else {
            // If no unvisited neighbors found, break the loop
            break;
        }
    }

    // Add the edge from the last vertex to the starting vertex to complete the tour
    if (!tour.empty()) {
        tour.push_back(tour.front());
    }

    return tour;
}

/*
 * template<class T>
double calculateDistance(Vertex<T>* v1, Vertex<T>* v2) {
    // Example using Euclidean distance
    double dx = v1->x - v2->x;
    double dy = v1->y - v2->y;
    return sqrt(dx * dx + dy * dy);
}
*/

template<class T>
double getEdgeWeight(Vertex<T>* v1, Vertex<T>* v2) {
    for (auto & e : v1->getAdj()) {
        if (e->getDest() == v2) {
            return e->getWeight();
        }
    }
    return numeric_limits<double>::infinity();
}

template<class T>
void NearestNeighbor(Graph<T> *g, const unordered_map<string, Edges *> *edges) {
    auto start = chrono::high_resolution_clock::now();

    vector<Vertex<T>*> tspTour = g->nearestNeighborTSP();

    auto end = chrono::high_resolution_clock::now();
    auto durationMicroseconds = chrono::duration_cast<chrono::microseconds>(end - start).count();
    auto durationMilliseconds = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    auto duration = chrono::duration_cast<chrono::seconds >(end - start).count();

    double totalDistance = 0.0;
    for (size_t i = 0; i < tspTour.size() - 1; i++) {
        totalDistance+= getEdgeWeight(tspTour[i], tspTour[i + 1]);
    }

    cout << "Shortest Path Found: ";
    for (unsigned int i = 0; i < tspTour.size(); ++i) {
        cout << tspTour.at(i)->getInfo();
        if (i != tspTour.size()-1) {
            cout << " -> ";
        }
    }
    cout << endl;

    cout << "Total distance:" << totalDistance << endl;

    if(duration != 0)
        cout << "Execution time: " << duration << " seconds" << endl;
    else if(durationMilliseconds != 0)
        cout << "Execution time: " << durationMilliseconds << " milliseconds" << endl;
    else
        cout << "Execution time: " << durationMicroseconds << " microseconds" << endl;
}

/*
static inline double convertToRadians(double degree) {
    return degree * M_PI / 180.0;
}

// Function to calculate Haversine distance between two coordinates in meters
static inline double haversineDistance(double lat1, double lon1, double lat2, double lon2) {
    double rad_lat1 = convertToRadians(lat1);
    double rad_lon1 = convertToRadians(lon1);
    double rad_lat2 = convertToRadians(lat2);
    double rad_lon2 = convertToRadians(lon2);

    double delta_lat = rad_lat2 - rad_lat1;
    double delta_lon = rad_lon2 - rad_lon1;

    double a = sin(delta_lat / 2.0) * sin(delta_lat / 2.0) +
               cos(rad_lat1) * cos(rad_lat2) *
               sin(delta_lon / 2.0) * sin(delta_lon / 2.0);

    double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

    double earth_radius = 6371000.0; // Earth radius in meters
    return earth_radius * c;
}


template<class T>
void NearestNeighbor(Graph<T> *g, const unordered_map<string, Edges *> *edges) {
    auto start = chrono::high_resolution_clock::now();

    vector<Vertex<T>*> tspTour = g->nearestNeighborTSP();

    auto end = chrono::high_resolution_clock::now();
    auto durationMicroseconds = chrono::duration_cast<chrono::microseconds>(end - start).count();
    auto durationMilliseconds = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    auto duration = chrono::duration_cast<chrono::seconds >(end - start).count();

    double totalDistance = 0.0;
    for (size_t i = 0; i < tspTour.size() - 1; i++) {
        // Extract coordinates of two vertices
        double lat1 = tspTour[i]->getLatitude();
        double lon1 = tspTour[i]->getLongitude();
        double lat2 = tspTour[i + 1]->getLatitude();
        double lon2 = tspTour[i + 1]->getLongitude();

        // Calculate distance between them using Haversine formula
        totalDistance += haversineDistance(lat1, lon1, lat2, lon2);
    }

    cout << "Shortest Path Found: ";
    for (unsigned int i = 0; i < tspTour.size(); ++i) {
        cout << tspTour.at(i)->getInfo();
        if (i != tspTour.size()-1) {
            cout << " -> ";
        }
    }
    cout << endl;

    cout << "Total distance:" << totalDistance << " meters" << endl;

    if(duration != 0)
        cout << "Execution time: " << duration << " seconds" << endl;
    else if(durationMilliseconds != 0)
        cout << "Execution time: " << durationMilliseconds << " milliseconds" << endl;
    else
        cout << "Execution time: " << durationMicroseconds << " microseconds" << endl;
}
*/

#endif // FEUP_DA_2024_PRJ2_G07_T14_GRAPH_H
