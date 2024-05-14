// Original code by Gonçalo Leão
// Updated by DA 2023/2024 Team

#ifndef FEUP_DA_2024_PRJ2_G07_T14_GRAPH_H
#define FEUP_DA_2024_PRJ2_G07_T14_GRAPH_H

#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <iomanip>
#include "data_structures/MutablePriorityQueue.h"


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
//template<class T>
//Graph<T>::Graph(const unordered_map<string, Vertex<T> *> & vertexSet):vertexSet(vertexSet){}


/******************************* EdmondsKarp *******************************/

/**
 * @brief Function to test the given vertex 'w' and visit it if conditions are met
 * @tparam T Type of vertex data
 * @param q Queue of vertices
 * @param e Edge connecting current vertex to 'w'
 * @param w Vertex to test and visit
 * @param residual Residual capacity of the edge
 */

// Time Complexity: O(1) - The function performs simple operations such as checking conditions and enqueueing the vertex 'w'.

// Function to test the given vertex 'w' and visit it if conditions are met
    template <class T>
    void testAndVisit(std::queue< Vertex<T>*> &q, Edge<T> *e, Vertex<T> *w, double residual) {
        // Check if the vertex 'w' is not visited and there is residual capacity
        if (! w->isVisited() && residual > 0) {
            // Mark 'w' as visited, set the path through which it was reached, and enqueue it
            w->setVisited(true);
            w->setPath(e);
            q.push(w);
        }
    }

/**
 * @brief Function to find an augmenting path using Breadth-First Search
 * @tparam T Type of vertex data
 * @param g Pointer to the graph
 * @param s Source vertex
 * @param t Target vertex
 * @return true if an augmenting path is found, false otherwise
 */

// Time Complexity: O(V + E), where V is the number of vertices and E is the number of edges in the graph. This is because the BFS algorithm traverses each vertex and each edge once.

// Function to find an augmenting path using Breadth-First Search
template <class T>
bool findAugmentingPath(Graph<T> *g, Vertex<T> *s, Vertex<T> *t) {
    // Mark all vertices as not visited
    for(auto v : g->getVertexSet()) {
        v.second->setVisited(false);
    }

    // Mark the source vertex as visited and enqueue it
    s->setVisited(true);
    std::queue<Vertex<T> *> q;
    q.push(s);

    // BFS to find an augmenting path
    while( ! q.empty() && ! t->isVisited()) {
        //std::cout << "queue: " << q.front() << " " << endl; // Assuming getData() returns the vertex data
        auto v = q.front();
        q.pop();

        // Process outgoing edges
        for(auto e: v->getAdj()) {
            testAndVisit(q, e, e->getDest(), e->getWeight() - e->getFlow());
        }

        // Process incoming edges
//        for(auto e: v->getIncoming()) {
//            testAndVisit(q, e, e->getOrig(), e->getFlow());
//            //cout << "OH O FLOW: " << e->getWeight() << endl;
//        }
    }

    // Return true if a path to the target is found, false otherwise
    return t->isVisited();
}


/**
 * @brief Function to find the minimum residual capacity along the augmenting path
 * @tparam T Type of vertex data
 * @param s Source vertex of the path
 * @param t Target vertex of the path
 * @return Minimum residual capacity along the augmenting path
 */

// Time Complexity: O(V), where V is the number of vertices in the augmenting path. This is because the function traverses the augmenting path once.

// Function to find the minimum residual capacity along the augmenting path
template <class T>
double findMinResidualAlongPath(Vertex<T> *s, Vertex<T> *t) {
    double f = INF;
    // Traverse the augmenting path to find the minimum residual capacity
    for (auto v = t; v != s; ) {
        auto e = v->getPath();
        if (e->getDest() == v) {
            f = std::min(f, e->getWeight() - e->getFlow());
            v = e->getOrig();
        }
        else {
            f = std::min(f, e->getFlow());
            v = e->getDest();
        }
    }

    // Return the minimum residual capacity
    return f;
}

/**
 * @brief Function to augment flow along the augmenting path with the given flow value
 * @tparam T Type of vertex data
 * @param s Source vertex of the path
 * @param t Target vertex of the path
 * @param f Flow value to augment
 */

// Time Complexity: O(V), where V is the number of vertices in the augmenting path. This is because the function traverses the augmenting path once.

// Function to augment flow along the augmenting path with the given flow value
template <class T>
void augmentFlowAlongPath(Vertex<T> *s, Vertex<T> *t, double f) {
    // Traverse the augmenting path and update the flow values accordingly
    for (auto v = t; v != s; ) {
        auto e = v->getPath();
        double flow = e->getFlow();
        if (e->getDest() == v) {
            e->setFlow(flow + f);
            v = e->getOrig();
        } else {
            e->setFlow(flow - f);
            v = e->getDest();
        }
    }
}


/**
 * @brief Main function implementing the Edmonds-Karp algorithm
 * @tparam T Type of vertex data
 * @param g Pointer to the graph
 * @param source Source vertex data
 * @param sink Sink vertex data
 */

// Time Complexity: O(VE^2), where V is the number of vertices and E is the number of edges in the graph. In each iteration, the BFS takes O(V + E) time, and there can be at most O(E) iterations.

// Main function implementing the Edmonds-Karp algorithm
template <class T>
void edmondsKarp(Graph<T> *g, string source, string sink) {
    // Find source and target vertices in the graph
    Vertex<T>* s = g->findVertex(source);
    Vertex<T>* t = g->findVertex(sink);

    // Validate source and target vertices
    if (s == nullptr || t == nullptr || s == t)
        throw std::logic_error("Invalid source and/or target vertex");


    // Initialize flow on all edges to 0
    for (auto v : g->getVertexSet()) {
        for (auto e: v.second->getAdj()) {
            e->setFlow(0);
        }
    }

    // While there is an augmenting path, augment the flow along the path
    while( findAugmentingPath(g, s, t) ) {
        double f = findMinResidualAlongPath(s, t);
        augmentFlowAlongPath(s, t, f);
    }


    for (auto & pair : g->getVertexSet()) {
        Vertex<T> *vertexFlow = pair.second;
        double incomingFlow = 0;
        for (auto e : vertexFlow->getIncoming()) {
            incomingFlow += e->getFlow();
        }
        vertexFlow->setFlow(incomingFlow);
    }
}

///**
// * @brief Function to calculate maximum flow in the network
// * @tparam T Type of vertex data
// * @param g Pointer to the graph
// * @param wr Map of reservoirs
// * @param ds Map of cities
// * @param i Parameter to choose output format
// */
//
//// Time Complexity: O(VE^2), where V is the number of vertices and E is the number of edges in the graph. This is because the function calls edmondsKarp, which has this time complexity.
//
//template<class T>
//void MaxFlow(Graph<T> *g, const unordered_map<string, Reservoirs *> *wr, const unordered_map<string, Cities *> *ds, int i) {
//    string source = "MasterSource";
//    string sink = "MasterSink";
//    g->addVertex(source);
//    g->addVertex(sink);
//    double maxFlow = 0;
//    double totalUndersupply= 0;
//    double undersupply= 0;
//    double demand2 = 0;
//    double flow = 0;
//
//    for (auto &pair: *wr) {
//        string rcode = pair.first;
//        Reservoirs *r = pair.second;
//        double weight = r->getMaxDelivery();
//
//        g->addEdge(source, rcode, pair.second->getMaxDelivery());
//    }
//
//    for (auto &pair: *ds) {
//        string dcode = pair.first;
//        Cities *c = pair.second;
//        double demand = c->getDemand();
//        double flow = 0;
//        g->addEdge(dcode, sink, pair.second->getDemand());
//    }
//
//    edmondsKarp(g, source, sink);
//
//    ofstream outputFile("output.txt");
//
//    if (outputFile.is_open()) { // check if the file was opened successfully
//        if(i == 1) {
//            outputFile << setw(20) << left << "City" << setw(20) << left << "Code" << setw(20) << left << "Value"
//                       << endl;
//            cout << setw(20) << left << "City" << setw(20) << left << "Code" << setw(20) << left << "Value" << endl;
//        }else if(i==2) {
//            cout << setw(20) << left << "City" << setw(20) << left << "Code" << setw(20) << left << "Demand" << setw(20) << left << "Value" << setw(20) << left << "Undersupply" << endl;
//        }
//        for (auto &pair: g->getVertexSet()) {
//            Vertex<T> *vertexDS = pair.second;
//            double flowTotal = 0;
//
//            for (auto e: vertexDS->getIncoming()) {
//                flowTotal += e->getFlow();
//            }
//
//        }
//        for (auto &pair: *ds) {
//            string code = pair.first;
//            Cities *cities = pair.second;
//
//            string name = cities->getCity();
//            double demand2 = cities->getDemand();
//            double flow = g->findVertex(code)->getFlow();
//            undersupply = demand2 - flow;
//            totalUndersupply += undersupply;
//            maxFlow += flow;
//            if(i==2) {
//                cout << setw(20) << left << name << setw(20) << left << code << setw(20) << left << demand2 << setw(20) << left << flow << setw(20) << left << undersupply << endl;
//
//            }else if(i==1){
//                outputFile << setw(20) << left << name << setw(20) << left << code << setw(20) << left << flow << endl;
//                cout << setw(20) << left << name << setw(20) << left << code << setw(20) << left << flow<< endl;
//            }
//        }
//
//        cout << endl << "Max Flow:" << maxFlow << endl;
//        if(i== 2) {
//            cout << "Total Undersupply:" << totalUndersupply << endl;
//        }
//        cout << "Data was written to output.txt\n";
//        outputFile.close();// close the file when done
//    } else {
//        std::cerr << "Error opening file\n";
//    }
//}
//
///**
// * @brief Function to simulate the impact of reservoir failure on the network
// * @tparam T Type of vertex data
// * @param g Pointer to the graph before failure
// * @param g2 Pointer to the graph after failure
// * @param wr Map of reservoirs
// * @param pipes Map of pipes
// * @param ds Map of cities
// * @param code Code of the reservoir to simulate failure
// */
//
//// Time Complexity: O(V + E), where V is the number of vertices and E is the number of edges in the graph. This is because the function calls edmondsKarp twice, which has this time complexity.
//
//template<class T>
//void reservoirFailure(Graph<T> *g, Graph<T> *g2, const unordered_map<string, Reservoirs *> *wr, const unordered_map<string, Pipes *> *pipes, const unordered_map<string, Cities *> *ds, const string& code) {
//    string source = "MasterSource";
//    string sink = "MasterSink";
//
//    g->addVertex(source);
//    g->addVertex(sink);
////    g2->addVertex(source);
////    g2->addVertex(sink);
//
//    for (const auto &pair : *wr) {
//        string rcode = pair.first;
//        Reservoirs *r = pair.second;
//        double weight = r->getMaxDelivery();
//        g->addEdge(source, rcode, pair.second->getMaxDelivery());
//        //g2->addEdge(source, rcode, pair.second->getMaxDelivery());
//    }
//    for (const auto &pair : *ds) {
//        string dcode = pair.first;
//        Cities *c = pair.second;
//        double demand = c->getDemand();
//        double flow = 0;
//        g->addEdge(dcode, sink, pair.second->getDemand());
//        //g2->addEdge(dcode, sink, pair.second->getDemand());
//    }
//    /*for (auto &pair: *pipes) {
//        string code2 = pair.first;
//        Pipes *pipe = pair.second;
//        if (pipe->getServicePointA() == code) {
//            cout << "Removing Pipeline " << code << " from " << pipe->getServicePointA() << " to " << pipe->getServicePointB() << endl;
//            g2->removeEdge(pipe->getServicePointA(), pipe->getServicePointB());
//        }
//    }*/
//    //g2->removeVertex(code);
//    edmondsKarp(g, source, sink);
//    //edmondsKarp(g2, source, sink);
//
//    // Store the initial flow for each city
//    unordered_map<string, double> initialFlow;
//    for (auto &pair: *ds) {
//        string code = pair.first;
//        Cities *cities = pair.second;
//        double flow = g->findVertex(code)->getFlow();
//        initialFlow[code] = flow;
//    }
//
//    g->removeVertex(code);
//
////    Pipes *pump = pipes->at(code);
////    //cout << "Removing Pumping Station " << pumpToRemove << " from " << pump->getServicePointA() << " to " << pump->getServicePointB() << endl;
////    g->removeEdge(pump->getServicePointA(), pump->getServicePointB());
//
//    edmondsKarp(g, source, sink);
////    int i = 0;
////    for (const auto &pair : *ds) {
////        string code3 = pair.first;
////        Cities *cities = pair.second;
////        string name = cities->getCity();
////        double flow1 = g->findVertex(code3)->getFlow();
////        double flow2 = g2->findVertex(code3)->getFlow();
////
////        if (flow1 != flow2) {
////            if (i == 0) {
////                cout << setw(20) << left << "City" << setw(20) << left << "Code" << setw(20) << left << "OldFlow" << setw(20) << left << "NewFlow" << endl;
////            }
////            i++;
////            cout << setw(20) << left << name << setw(20) << left << code3 << setw(20) << left << flow1 << setw(20) << left << flow2 << endl;
////        }
////    }
////
////    if (i == 0) {
////        cout << "No impact on the flow of the cities";
////    }
//
//
//    // Output the impact of the Reservoir failure
//    bool anyCityAffected = false;
//    for (auto &pair: *ds) {
//
//        string code = pair.first;
//        Cities *cities = pair.second;
//        double newFlow = g->findVertex(code)->getFlow();
//        double oldFlow = initialFlow[code];
//
//        if (newFlow != oldFlow) {
//            anyCityAffected = true;
//            cout << "City " << code << " is affected:" << endl;
//            cout << "Old Flow: " << oldFlow << endl;
//            cout << "New Flow: " << newFlow << endl;
//        }
//    }
//
//    if (!anyCityAffected) {
//        cout << "No city is affected." << endl;
//    }
//}
//
///**
// * @brief Function to simulate the impact of pumping station failure on the network
// * @tparam T Type of vertex data
// * @param g Pointer to the graph
// * @param wr Map of reservoirs
// * @param ds Map of cities
// * @param pipes Map of pipes
// * @param pumpToRemove Code of the pumping station to remove
// */
//
//// Time Complexity: O(V + E), where V is the number of vertices and E is the number of edges in the graph. This is because the function calls edmondsKarp twice, which has this time complexity.
//
//template<class T>
//void pumpingStationFailureImpact(Graph<T> *g, const unordered_map<string, Reservoirs *> *wr, const unordered_map<string, Cities *> *ds, const unordered_map<string, Pipes *> *pipes, const string& pumpToRemove) {
//    string source = "MasterSource";
//    string sink = "MasterSink";
//    g->addVertex(source);
//    g->addVertex(sink);
//    double totalUndersupply = 0;
//    double undersupply = 0;
//
//    // Add all reservoirs and their outgoing edges to the graph
//    for (auto &pair: *wr) {
//        string rcode = pair.first;
//        Reservoirs *r = pair.second;
//        double weight = r->getMaxDelivery();
//        g->addEdge(source, rcode, pair.second->getMaxDelivery());
//    }
//
//    // Add all cities and their incoming edges to the graph
//    for (auto &pair: *ds) {
//        string dcode = pair.first;
//        Cities *c = pair.second;
//        double demand = c->getDemand();
//        g->addEdge(dcode, sink, pair.second->getDemand());
//    }
//
//    // Perform max flow calculation
//    edmondsKarp(g, source, sink);
//
//    // Store the initial flow for each city
//    unordered_map<string, double> initialFlow;
//    for (auto &pair: *ds) {
//        string code = pair.first;
//        Cities *cities = pair.second;
//        double flow = g->findVertex(code)->getFlow();
//        initialFlow[code] = flow;
//    }
//
//    // Check if the selected pump exists
//    if (pipes->find(pumpToRemove) == pipes->end()) {
//        cout << "Error: Pumping Station " << pumpToRemove << " not found." << endl;
//        return;
//    }
//
//    // Remove the selected pump and recalculate the max flow
//    Pipes *pump = pipes->at(pumpToRemove);
//    //cout << "Removing Pumping Station " << pumpToRemove << " from " << pump->getServicePointA() << " to " << pump->getServicePointB() << endl;
//    g->removeEdge(pump->getServicePointA(), pump->getServicePointB());
//
//    edmondsKarp(g, source, sink);
//
//    // Output the impact of the pumping station failure
//    bool anyCityAffected = false;
//    for (auto &pair: *ds) {
//        string code = pair.first;
//        Cities *cities = pair.second;
//        double newFlow = g->findVertex(code)->getFlow();
//        double oldFlow = initialFlow[code];
//
//        if (newFlow != oldFlow) {
//            anyCityAffected = true;
//            cout << "City " << code << " is affected:" << endl;
//            cout << "Old Flow: " << oldFlow << endl;
//            cout << "New Flow: " << newFlow << endl;
//        }
//    }
//
//    if (!anyCityAffected) {
//        cout << "No city is affected." << endl;
//    }
//}
//
///**
// * @brief Function to simulate the impact of pipeline failure on the network
// * @tparam T Type of vertex data
// * @param g Pointer to the graph
// * @param wr Map of reservoirs
// * @param ds Map of cities
// * @param pipes Map of pipes
// * @param pipeToRemove Vector of pipeline codes to remove
// */
//
//// Time Complexity: O(P * (V + E)), where P is the number of pipelines to be removed, V is the number of vertices, and E is the number of edges in the graph.
//// This is because the function iterates over each pipeline to be removed and calls edmondsKarp once for each removal, which has a time complexity of O(V + E).
//
//template<class T>
//void pipelineFailureImpact(Graph<T> *g, const unordered_map<string, Reservoirs *> *wr, const unordered_map<string, Cities *> *ds, const unordered_map<string, Pipes *> *pipes, const vector<string>& pipeToRemove) {
//    string source = "MasterSource";
//    string sink = "MasterSink";
//    g->addVertex(source);
//    g->addVertex(sink);
//    double totalUndersupply = 0;
//    double undersupply = 0;
//
//    // Add all reservoirs and their outgoing edges to the graph
//    for (auto &pair: *wr) {
//        string rcode = pair.first;
//        Reservoirs *r = pair.second;
//        double weight = r->getMaxDelivery();
//        g->addEdge(source, rcode, pair.second->getMaxDelivery());
//    }
//
//    // Add all cities and their incoming edges to the graph
//    for (auto &pair: *ds) {
//        string dcode = pair.first;
//        Cities *c = pair.second;
//        double demand = c->getDemand();
//        g->addEdge(dcode, sink, pair.second->getDemand());
//    }
//
//    // Perform max flow calculation
//    edmondsKarp(g, source, sink);
//
//    // Store the initial flow for each city
//    unordered_map<string, double> initialFlow;
//    for (auto &pair: *ds) {
//        string code = pair.first;
//        Cities *cities = pair.second;
//        double flow = g->findVertex(code)->getFlow();
//        initialFlow[code] = flow;
//    }
//
//    // Loop through each pipeline to be removed
//    for (const string& codePipe : pipeToRemove) {
//        // Check if the selected pipe exists
//        if (pipes->find(codePipe) == pipes->end()) {
//            cout << "Error: Pipeline " << codePipe << " not found." << endl;
//            continue; // Move to the next pipeline if this one doesn't exist
//        }
//
//        // Remove the selected pipeline and recalculate the max flow
//        bool pipeRemoved = false;
//        for (auto &pair: *pipes) {
//            string code = pair.first;
//            Pipes *pipe = pair.second;
//            if (code == codePipe) {
//                cout << "Removing Pipeline " << codePipe << " from " << pipe->getServicePointA() << " to " << pipe->getServicePointB() << endl;
//                g->removeEdge(pipe->getServicePointA(), pipe->getServicePointB());
//                //g->removeEdge(pipe->getServicePointB(), pipe->getServicePointA());
//                pipeRemoved = true;
//                break; // Exit the loop after removing the pipeline
//            }
//        }
//
//        if (!pipeRemoved) {
//            // This shouldn't happen if the pipeline exists, but it's good to check
//            cout << "Error: Failed to remove pipeline " << codePipe << endl;
//        }
//    }
//
//    edmondsKarp(g, source, sink);
//
//    // Output the impact of the pipeline failure
//    bool anyCityAffected = false;
//    for (auto &pair: *ds) {
//
//        string code = pair.first;
//        Cities *cities = pair.second;
//        double newFlow = g->findVertex(code)->getFlow();
//        double oldFlow = initialFlow[code];
//
//        if (newFlow != oldFlow) {
//            anyCityAffected = true;
//            cout << "City " << code << " is affected:" << endl;
//            cout << "Old Flow: " << oldFlow << endl;
//            cout << "New Flow: " << newFlow << endl;
//        }
//    }
//
//    if (!anyCityAffected) {
//        cout << "No city is affected." << endl;
//    }
//}
//template<class T>
//void pipelineFailureImpact(Graph<T> *g, const unordered_map<string, Reservoirs *> *wr, const unordered_map<string, Cities *> *ds, const unordered_map<string, Pipes *> *pipes, const vector<string>& pipeToRemove) {
//    string source = "MasterSource";
//    string sink = "MasterSink";
//    g->addVertex(source);
//    g->addVertex(sink);
//    double totalUndersupply = 0;
//    double undersupply = 0;
//
//    // Add all reservoirs and their outgoing edges to the graph
//    for (auto &pair: *wr) {
//        string rcode = pair.first;
//        Reservoirs *r = pair.second;
//        double weight = r->getMaxDelivery();
//        g->addEdge(source, rcode, pair.second->getMaxDelivery());
//    }
//
//    // Add all cities and their incoming edges to the graph
//    for (auto &pair: *ds) {
//        string dcode = pair.first;
//        Cities *c = pair.second;
//        double demand = c->getDemand();
//        g->addEdge(dcode, sink, pair.second->getDemand());
//    }
//
//    // Perform max flow calculation
//    edmondsKarp(g, source, sink);
//
//    // Store the initial flow for each city
//    unordered_map<string, double> initialFlow;
//    for (auto &pair: *ds) {
//        string code = pair.first;
//        Cities *cities = pair.second;
//        double flow = g->findVertex(code)->getFlow();
//        initialFlow[code] = flow;
//    }
//
//    // Loop through each pipeline to be removed
//    for (const string& codePipe : pipeToRemove) {
//        // Check if the selected pipe exists
//        if (pipes->find(codePipe) == pipes->end()) {
//            cout << "Error: Pipeline " << codePipe << " not found." << endl;
//            continue; // Move to the next pipeline if this one doesn't exist
//        }
//
//        // Remove the selected pipeline and recalculate the max flow
//        for (auto &pair: *pipes) {
//            string code = pair.first;
//            Pipes *pipe = pair.second;
//            if (code == codePipe) {
//                cout << "From Point:" << pipe->getServicePointA() << " to "  << "Point:" << pipe->getServicePointB() << endl;
//                g->removeEdge(pipe->getServicePointA(), pipe->getServicePointB());
//            }
//        }
//    }
//
//    edmondsKarp(g, source, sink);
//
//    // Check if any city is affected by the pipeline failure
//    bool cityAffected = false;
//    for (auto &pair: *ds) {
//        string code = pair.first;
//        Cities *cities = pair.second;
//        double newFlow = g->findVertex(code)->getFlow();
//        double oldFlow = initialFlow[code];
//
//        if (newFlow != oldFlow) {
//            cityAffected = true;
//        }
//    }
//        // Output the impact of the pipeline failure
//        if (cityAffected) {
//            for (const string& codePipe1 : pipeToRemove) {
//                cout << "Case: Pipeline " << codePipe1 << " is removed. Cities affected:" << endl;
//                cout << setw(20) << left << "City" << setw(20) << left << "Old Flow" << setw(20) << left << "New Flow"
//                     << endl;
//            }
//            for (auto &pair: *ds) {
//                string code = pair.first;
//                Cities *cities = pair.second;
//                double newFlow = g->findVertex(code)->getFlow();
//                double oldFlow = initialFlow[code];
//                if (newFlow != oldFlow) {
//                    cout << setw(20) << left << code << setw(20) << left << oldFlow << setw(20) << left << newFlow
//                         << endl;
//                }
//            }
//            cout << endl;
//        } else {
//            for (const string& codePipe1 : pipeToRemove) {
//                cout << "Case: Pipeline " << codePipe1 << " is removed. No city is affected." << endl << endl;
//            }
//        }
//
//}

/*
template<class T>
void pipelineFailureImpact(Graph<T> *g, const unordered_map<string, Reservoirs *> *wr, const unordered_map<string, Cities *> *ds, const unordered_map<string, Pipes *> *pipes, const string& pipeToRemove) {
    string source = "MasterSource";
    string sink = "MasterSink";
    g->addVertex(source);
    g->addVertex(sink);
    double totalUndersupply = 0;
    double undersupply = 0;

    // Add all reservoirs and their outgoing edges to the graph
    for (auto &pair: *wr) {
        string rcode = pair.first;
        Reservoirs *r = pair.second;
        double weight = r->getMaxDelivery();

        g->addEdge(source, rcode, pair.second->getMaxDelivery());
    }

    // Add all cities and their incoming edges to the graph
    for (auto &pair: *ds) {
        string dcode = pair.first;
        Cities *c = pair.second;
        double demand = c->getDemand();
        g->addEdge(dcode, sink, pair.second->getDemand());
    }

    // Perform max flow calculation
    edmondsKarp(g, source, sink);

    // Store the initial flow for each city
    unordered_map<string, double> initialFlow;
    for (auto &pair: *ds) {
        string code = pair.first;
        Cities *cities = pair.second;
        double flow = g->findVertex(code)->getFlow();
        initialFlow[code] = flow;
    }

    // Check if the selected pipe exists
    if (pipes->find(pipeToRemove) == pipes->end()) {
        cout << "Error: Pipeline " << pipeToRemove << " not found." << endl;
        return;
    }

    for(auto &pair : *pipes) {
        string code = pair.first;
        Pipes *pipe = pair.second;
        if(code == pipeToRemove){
            g->removeEdge(pipe->getServicePointA(), pipe->getServicePointB());
        }
    }
    // Remove the selected pipeline and recalculate the max flow
    edmondsKarp(g, source, sink);

    // Check if any city is affected by the pipeline failure
    bool cityAffected = false;
    for (auto &pair: *ds) {
        string code = pair.first;
        Cities *cities = pair.second;
        double demand = cities->getDemand();
        double newFlow = g->findVertex(code)->getFlow();
        double oldFlow = initialFlow[code];
        undersupply = demand - newFlow;

        // Check if the new flow is less than the demand
//        if (undersupply > 0) {
//            cityAffected = true;
//            totalUndersupply += undersupply;
//        }
        if(newFlow != oldFlow){
            cityAffected = true;
        }
    }

    // Output the impact of the pipeline failure
    if (cityAffected) {
        cout << "Case: Pipeline " << pipeToRemove << " is removed. Cities affected:" << endl;
        cout << setw(20) << left << "City" << setw(20) << left << "Old Flow" << setw(20) << left << "New Flow" << endl;
        for (auto &pair: *ds) {
            string code = pair.first;
            Cities *cities = pair.second;
            double newFlow = g->findVertex(code)->getFlow();
            double oldFlow = initialFlow[code];
            undersupply = cities->getDemand() - newFlow;
            if (newFlow != oldFlow) {
                cout << setw(20) << left << code << setw(20) << left << oldFlow << setw(20) << left << newFlow << endl;
            }
        }
        cout << endl;
    } else {
        cout << "Case: Pipeline " << pipeToRemove << " is removed. No city is affected." << endl << endl;
    }

// Output the total undersupply caused by the pipeline failure
// cout << "Total Undersupply:" << totalUndersupply << endl;

//    Restore the pipeline for the next iteration
//    g->addEdge(source, pipeToRemove, (*wr)[pipeToRemove]->getMaxDelivery());
}
*/

/*
template<class T>
void WaterNetworkUndersupply(Graph<T> *g, const unordered_map<string, Reservoirs *> *wr, const unordered_map<string, Cities *> *ds) {
    string source = "MasterSource";
    string sink = "MasterSink";
    g->addVertex(source);
    g->addVertex(sink);
    double maxFlow = 0;

    for (auto &pair: *wr) {
        string rcode = pair.first;
        Reservoirs *r = pair.second;
        double weight = r->getMaxDelivery();

        g->addEdge(source, rcode, pair.second->getMaxDelivery());
    }

    for (auto &pair: *ds) {
        string dcode = pair.first;
        Cities *c = pair.second;
        double demand = c->getDemand();

        g->addEdge(dcode, sink, pair.second->getDemand());
    }

    edmondsKarp(g, source, sink);


        for (auto &pair: g->getVertexSet()) {
            Vertex<T> *vertexDS = pair.second;
            double flowTotal = 0;

            for (auto e: vertexDS->getIncoming()) {
                flowTotal += e->getFlow();
            }

        }
        for (auto &pair: *ds) {
            string code = pair.first;
            Cities *cities = pair.second;

            string name = cities->getCity();
            double demand = cities->getDemand();
            double flow = g->findVertex(code)->getFlow();
            //double teste = g->findVertex(code)->getFlow();
            maxFlow += flow;

            cout << name << ", " << code << ", " << flow << endl;
        }

        cout << endl << "Max Flow:" << maxFlow << endl;
        cout << "Data was written to output.txt\n";


}
*/




#endif // FEUP_DA_2024_PRJ2_G07_T14_GRAPH_H