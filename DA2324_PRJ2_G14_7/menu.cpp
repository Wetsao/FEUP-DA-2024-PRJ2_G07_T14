#include "menu.h"

void Menu::getDisplay(){
    string state = "start";
    Graph<string> Portugal;
    DataReader reader;
    while(!(state =="finished")) {
        int instruction;

        if(state == "start") {
            cout << "Analysis for Water Supply Management" << endl;
            cout << "1: Madeira's Network" << endl;
            cout << "2: Portugal's Network" << endl;
            cin >> instruction;
            if (instruction == 1) {
                reader.getCities("../Project1DataSetSmall/Cities_Madeira.csv");
                reader.getReservoirs("../Project1DataSetSmall/Reservoirs_Madeira.csv");
                reader.getStations("../Project1DataSetSmall/Stations_Madeira.csv");
                reader.getPipes("../Project1DataSetSmall/Pipes_Madeira.csv");

                state = "dataLoaded";
            }
            if (instruction == 2) {
                reader.getCities("../Project1LargeDataSet/Cities.csv");
                reader.getReservoirs("../Project1LargeDataSet/Reservoir.csv");
                reader.getStations("../Project1LargeDataSet/Stations.csv");
                reader.getPipes("../Project1LargeDataSet/Pipes.csv");

                state = "dataLoaded";
            }

        }
        if(state == "dataLoaded"){
            cout << "1: Get Max Flow" << endl;
            cout << "2: Water Supply Undersupply" << endl;
            cout << "3: Pipeline Failure Impact" << endl;
            cout << "4: Reservoirs Failure Impact" << endl;
            cout << "5: Pumping Station Failure Impact" << endl;
            cin >> instruction;
            if(instruction == 1){
                state = "selectMode";

            }
            if(instruction == 2){
                reader.MaxFlowAllCities(2);
                state="finished";
            }
            if (instruction == 3) {
                vector<string> codes;
                cout << "Select The Pipe To Remove (type 'done' when finished):" << endl;
                string pipeCode;
                while (true) {
                    cout << "Pipe code: ";
                    cin >> pipeCode;
                    if (pipeCode == "done") {
                        break; // Exit the loop if the user inputs 'done'
                    }
                    codes.push_back(pipeCode); // Add the input pipe code to the vector
                }

                // Now you have all the pipe codes in the 'codes' vector
                // You can pass this vector to the function pipelineFailureImpact
                reader.pipelinefailureImpact(codes);


//            if(instruction == 3){
//                vector<string> code;
//                cout << "Select The Pipe To Remove:";
//                cin >> code;
//                reader.pipelinefailureImpact(code);
//                state="finished";
            }
            if(instruction== 4){
                string code;
                cout << "Select the Reservoir:";
                cin >> code;
                reader.reservoirfailure(code);
                state="finished";
            }
            if(instruction == 5){
                string code;
                cout << "Select the Pumping Station:";
                cin >>code;
                reader.pumpingstationfailure(code);
                state="finished";
            }
        }
        if(state == "selectMode"){
            cout << "1: All Cities" << endl;
            cout << "2: Specific City" << endl;
            cout << "0: Back" << endl;
            cin >> instruction;
            if(instruction==0){
                state="dataLoaded";
            }
            if(instruction==1){
                reader.MaxFlowAllCities(1);
                state="finished";
            }
            if(instruction==2){
                string name;
                cout << "Please Introduce the city code" << endl;
                cin >> name;
                reader.MaxFlowSpecificCity(name);
                state="finished";
            }
        }
        if(state == "finished"){
            cout << endl << "Execution Finished" << endl;
        }
    }
}
