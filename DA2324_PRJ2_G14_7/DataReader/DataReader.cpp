#include "DataReader.h"

DataReader::DataReader() {};

void DataReader::getCities(const string &filename) {
    string origen, destino, distancia;


    ifstream file(filename);


    if (!file.is_open())
    {
        cerr << "Error: file " << filename << " not found" << endl;
        return;
    }

    string line;
    getline(file, line); // Ignore header line
    while (getline(file, line)) {
        stringstream ss(line);

        getline(ss, City, ',');
        getline(ss, Id, ',');
        getline(ss, Code, ',');
        getline(ss, Demand, ',');
        getline(ss, Population, ',');

        try {
            int idi = stoi(Id);
            double dem = stod(Demand);
            double pop = stod(Population);

            auto* cities = new Cities(City, idi, Code, dem, pop);
            this->DS.insert({Code, cities});
            g.addVertex(Code);
        } catch (const std::invalid_argument& e) {
            cerr << "Invalid argument: " << e.what() << endl;
        } catch (const std::out_of_range& e) {
            cerr << "Out of range: " << e.what() << endl;
        }
    }
    file.close();
}