#ifndef FEUP_DA_2024_PRJ2_G07_T14_MENU_H
#define FEUP_DA_2024_PRJ2_G07_T14_MENU_H

#include <sstream>
#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include "data_structures/Graph.h"
#include "DataReader/DataReader.h"
#include "data_structures/Graph.h"

using namespace std;

class Menu {

public:
    static int MainMenu(); //main menu
    static bool is_number(const string &s);
    static void  returnToMainMenu();

    void getDisplay();
};

#endif //FEUP_DA_2024_PRJ2_G07_T14_MENU_H
