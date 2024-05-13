#ifndef FEUP_DA_2024_PRJ2_G07_T14_DATAREADER_H
#define FEUP_DA_2024_PRJ2_G07_T14_DATAREADER_H

#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <limits>
#include <sstream>
#include "data_structures/Graph.h"

using namespace std;

class DataReader {
private:
    unordered_map<string, Reservoirs *> WR;
    unordered_map<string, Stations *> PS;
    unordered_map<string, Cities *> DS;
    unordered_map<string, Pipes *> pipes;
    Graph<string> g;
public:
    DataReader();


    void getCities(const string &filename);
    void getReservoirs(const string &filename);
    void getStations(const string &filename);
    void getPipes(const string &filename);

    const unordered_map<string, Reservoirs *> &getWR() const;
    const unordered_map<string, Stations *> &getPS() const;
    const unordered_map<string, Cities *> &getDS() const;
    const unordered_map<string, Pipes *> &getPipesgetter() const;

    void MaxFlowAllCities(int i);
    void MaxFlowSpecificCity(const string &code);

    void pipelinefailureImpact(const vector<string>& pipeToRemove);
    void reservoirfailure(const string& code);

    void pumpingstationfailure(const string &code);
};



#endif //FEUP_DA_2024_PRJ2_G07_T14_DATAREADER_H
