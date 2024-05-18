#include "DataReader.h"

DataReader::DataReader() {}

void DataReader::getEdgesFile(const string &filename) {
    int i= 0;
    string origem, destino, distancia;

    ifstream file(filename);

    if (!file)
    {
        cerr << "Error: file " << filename << " not found" << endl;
    }
    string line;
    getline(file, line); // Ignore header line
    while (getline(file, line)) {
        stringstream ss(line);

        getline(ss, origem, ',');
        getline(ss, destino, ',');
        getline(ss, distancia, ',');

        try {
//            int idorigem = stoi(origem);
//            int iddestino = stoi(destino);
            double dist = stod(distancia);

            auto* edges = new Edges(origem, destino, dist);
            this->edges.insert({to_string(i),edges});
            auto* edges2 = new Edges(destino, origem, dist);
            this->edges.insert({to_string(++i), edges2});
            g.addVertex(origem);
            g.addVertex(destino);
            g.addBidirectionalEdge(origem, destino, dist);
            //g.addEdge(origem, destino, dist);
        } catch (const std::invalid_argument& e) {
            cerr << "Invalid argument: " << e.what() << endl;
        } catch (const std::out_of_range& e) {
            cerr << "Out of range: " << e.what() << endl;
        }
        i++;
    }
    file.close();
}

const unordered_map<string, Edges *> &DataReader::getEdges() const {
    return edges;
}

void DataReader::cenario41() {
    Backtracking(&g, &getEdges());
}

void DataReader::cenario42() {
    //TriangleApproximation(&g, &getEdges());
}

void DataReader::cenario43() {
    //OtherHeuristic(&g, &getEdges());
}

void DataReader::cenario44() {
}
