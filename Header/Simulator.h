#pragma once
#ifndef __SIMULATOR_H__
#define __SIMULATOR_H__
#include <iostream>
#include <vector>
#include <string>
using namespace std;
int final_mult(int mult1, int mult2);
void sim_h();
void sim_f();
vector<double> max_error_AC42_propagation(vector<double> probs);
vector<double> max_error_AC32_propagation(vector<double> probs);
void ConstructCT(vector<vector<int>>& F, vector<vector<int>>& H, vector<vector<int>>& AC32,
	vector<vector<int>>& AC42, vector<vector<double>>& E, vector<vector<int>>& V);
void generate_MEDs(vector<double>& MED_ac42, vector<double>& MED_ac32);
#endif // !__SIMULATOR_H__

