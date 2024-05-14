#include "Edges.h"


Edges::Edges(const string &initialNode, const string &finalNode, double distancia) : initial_Node(initialNode),
                                                                                     final_Node(finalNode),
                                                                                     distancia(distancia) {}

const string &Edges::getInitialNode() const {
    return initial_Node;
}

void Edges::setInitialNode(const string &initialNode) {
    initial_Node = initialNode;
}

const string &Edges::getFinalNode() const {
    return final_Node;
}

void Edges::setFinalNode(const string &finalNode) {
    final_Node = finalNode;
}

double Edges::getDistancia() const {
    return distancia;
}

void Edges::setDistancia(double distancia) {
    Edges::distancia = distancia;
}
