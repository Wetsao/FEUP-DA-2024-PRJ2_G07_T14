#include "menu.h"

void Menu::getDisplay(){
    int state = 1;
    // State = 1 Means Start
    // State = 2 Means Data Loaded
    // State = 3 Means Finished Program
    DataReader reader;
    string datasetPath;
    while(!(state ==3)) {
            switch(state){
                case 1:
                    cout << "You Can Select The Dataset By Typing The Path To File Of The Dataset (EX: C:\\Users\\User\\Desktop\\e\\2324\\2ano\\2sem\\DA\\FEUP-DA-2024-PRJ2_G07_T14\\DA2324_PRJ2_G14_7\\DataSet\\Toy-Graphs\\shipping.csv)" << endl;
                    cout << "Select Dataset:";
                    cin >> datasetPath;
                    reader.getEdgesFile(datasetPath);
                    state = 2;
                    break;
                case 2:
                    int instruction2;
                    cout << "Select Approach:" << endl;
                    cout << "1: Bracktracking" << endl;
                    cout << "2: Triangle Approximation" << endl;
                    cout << "3: Nearest Neighbor" << endl;
                    cout << "4: TSP in the Real World" << endl;
                    cin >> instruction2;
                    switch(instruction2){
                        case 1:
                            reader.cenario41();
                            state++;
                            break;
                        case 2:
                            state++;
                            break;
                        case 3:
                            reader.cenario43();
                            state++;
                            break;
                        case 4:
                            state++;
                            break;
                        default:
                            cout << "Error Selecting The Aprroach Please Try Again" << endl;
                            break;
                    }
                    break;
            }
            cout << "Program Finished Executing" << endl;
    }
}
