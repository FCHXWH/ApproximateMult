#pragma once
#ifndef MULTIPLIER_TREE
#define MULTIPLIER_TREE
#include <iostream>
#include <vector>
#include <sstream>
#include "gurobi_c++.h"
using namespace std;
extern int MULT_SIZE;
extern int stages_num;
extern string filename;
extern vector<int> input_patterns;
void generate_input_patterns();
void generate_variables_of_multiplier(vector<vector<GRBVar>>& variables_fa, vector<vector<GRBVar>>& variables_ha,
	vector<vector<GRBVar>>& variables_V, vector<vector<GRBVar>>& variables_ac32,
	vector<vector<GRBVar>>& variables_ac42, vector<vector<GRBVar>>& variables_error, GRBModel& model);
void Initial_constraints_of_V(vector<vector<GRBVar>> variables_V, GRBModel& model);
void generate_constraints_of_h_f_V(vector<vector<GRBVar>> variables_fa, vector<vector<GRBVar>> variables_ha, vector<vector<GRBVar>> variables_V,
	vector<vector<GRBVar>> variables_ac32, vector<vector<GRBVar>> variables_ac42, vector<vector<GRBVar>> variables_error, vector<double> MED_ac42,
	vector<double> MED_ac32, GRBModel& model);
void generate_constraints_of_last_stage(vector<vector<GRBVar>> variables_V, GRBModel& model);
void generate_cost_of_adders(vector<vector<GRBVar>> variables_fa, vector<vector<GRBVar>> variables_ha, GRBVar& fa_num, GRBVar& ha_num,
	vector<vector<GRBVar>> variables_ac32, vector<vector<GRBVar>> variables_ac42, GRBVar& ac32_num, GRBVar& ac42_num, GRBModel& model);
void generate_total_error(vector<vector<GRBVar>> variables_error, GRBVar& E, GRBModel& model);


#endif // !MULTIPLIER_TREE
