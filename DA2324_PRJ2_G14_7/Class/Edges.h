#ifndef DA2324_PRJ2_G14_7_EDGES_H
#define DA2324_PRJ2_G14_7_EDGES_H

#include <iostream>
#include <string>

using namespace std;

class Edges {
private:
    string initial_Node;
    string final_Node;
    double distancia;
public:
    Edges(const string &initialNode, const string &finalNode, double distancia);

    const string &getInitialNode() const;
    void setInitialNode(const string &initialNode);
    const string &getFinalNode() const;
    void setFinalNode(const string &finalNode);
    double getDistancia() const;
    void setDistancia(double distancia);
};


#endif //DA2324_PRJ2_G14_7_EDGES_H
