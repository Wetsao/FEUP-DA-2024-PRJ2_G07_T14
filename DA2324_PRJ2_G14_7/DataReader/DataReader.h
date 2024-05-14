#ifndef FEUP_DA_2024_PRJ2_G07_T14_DATAREADER_H
#define FEUP_DA_2024_PRJ2_G07_T14_DATAREADER_H

#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <limits>
#include <sstream>
#include "data_structures/Graph.h"
#include "Class/Edges.h"


using namespace std;

class DataReader {
private:
    unordered_map<string, Edges *> edges;
    Graph<string> g;
public:
    DataReader();

    void getEdgesFile(const string &filename);

    const unordered_map<string, Edges *> &getEdges() const;

    void cenario41();
    void cenario42();
    void cenario43();
    void cenario44();
};



#endif //FEUP_DA_2024_PRJ2_G07_T14_DATAREADER_H
