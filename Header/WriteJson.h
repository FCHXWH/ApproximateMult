#pragma once
#ifndef __WRITE_JSON__
#define __WRITE_JSON__
#include <iostream>
#include "wireConnect.h"
#include <json/json.h>
#include <fstream>
#include <string>
using namespace std;
using namespace Json;

void WriteToJson(vector<vector<column_sol_t>> allSols, String file);
#endif // !__WRITE_JSON__
                               