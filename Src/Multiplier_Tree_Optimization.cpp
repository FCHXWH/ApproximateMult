/* Copyright 2018, Gurobi Optimization, LLC */

#include "gurobi_c++.h"
#include "Multiplier_Tree_Optimization.h"

void generate_input_patterns()
{
	//Construct input pattern
	for (int i = 1; i <= 2 * MULT_SIZE - 1; i++)
	{
		input_patterns.push_back(MULT_SIZE - abs(i - MULT_SIZE));
		//cout << *(input_patterns.end() - 1) << ",";
	}
}

void generate_variables_of_multiplier(vector<vector<GRBVar>>& variables_fa, vector<vector<GRBVar>>& variables_ha,
	vector<vector<GRBVar>>& variables_V, vector<vector<GRBVar>>& variables_ac32,
	vector<vector<GRBVar>>& variables_ac42, vector<vector<GRBVar>>& variables_error, GRBModel& model)
{
	for (int i = 0; i < stages_num + 1; i++)
	{
		vector<GRBVar> variables_fa_of_each_stage;
		vector<GRBVar> variables_ha_of_each_stage;
		vector<GRBVar> variables_ac32_of_each_stage;
		vector<GRBVar> variables_ac42_of_each_stage;
		vector<GRBVar> variables_error_of_each_stage;
		vector<GRBVar> variables_V_of_each_stage;
		stringstream ss;
		string tmp_s;
		for (int j = 0; j < input_patterns.size() + stages_num; j++)
		{
			if (i != 0)
			{
				ss << "f" << i << "_" << j;
				ss >> tmp_s;
				ss.clear(); ss.str("");
				GRBVar variable_of_fa = model.addVar(0, INFINITY, 0, GRB_INTEGER, tmp_s);
				variables_fa_of_each_stage.push_back(variable_of_fa);

				ss << "h" << i << "_" << j;
				ss >> tmp_s;
				ss.clear(); ss.str("");
				GRBVar variable_of_ha = model.addVar(0, INFINITY, 0, GRB_INTEGER, tmp_s);
				variables_ha_of_each_stage.push_back(variable_of_ha);

				ss << "32ac" << i << "_" << j;
				ss >> tmp_s;
				ss.clear(); ss.str("");
				GRBVar variable_of_32ac = model.addVar(0, INFINITY, 0, GRB_INTEGER, tmp_s);
				variables_ac32_of_each_stage.push_back(variable_of_32ac);

				ss << "42ac" << i << "_" << j;
				ss >> tmp_s;
				ss.clear(); ss.str("");
				GRBVar variable_of_42ac = model.addVar(0, INFINITY, 0, GRB_INTEGER, tmp_s);
				variables_ac42_of_each_stage.push_back(variable_of_42ac);

				ss << "e" << i << "_" << j;
				ss >> tmp_s;
				ss.clear(); ss.str("");
				GRBVar variable_of_e = model.addVar(0, INFINITY, 0, GRB_CONTINUOUS, tmp_s);
				variables_error_of_each_stage.push_back(variable_of_e);
			}

			ss << "V" << i << "_" << j;
			ss >> tmp_s;
			ss.clear();
			GRBVar variable_of_V = model.addVar(0, INFINITY, 0, GRB_INTEGER, tmp_s);
			variables_V_of_each_stage.push_back(variable_of_V);
		}
		if (i != 0)
		{
			variables_fa.push_back(variables_fa_of_each_stage);
			variables_ha.push_back(variables_ha_of_each_stage);
			variables_ac32.push_back(variables_ac32_of_each_stage);
			variables_ac42.push_back(variables_ac42_of_each_stage);
			variables_error.push_back(variables_error_of_each_stage);
		}
		variables_V.push_back(variables_V_of_each_stage);
		
	}

}

void Initial_constraints_of_V(vector<vector<GRBVar>> variables_V, GRBModel& model)
{
	for (int i = 0; i < input_patterns.size() + stages_num; i++)
	{
		if (i < stages_num)
		{
			model.addConstr(variables_V[0][i] == 0);
		}
		else
		{
			model.addConstr(variables_V[0][i] == input_patterns[i - stages_num]);
		}
		
	}
}

void generate_constraints_of_h_f_V(vector<vector<GRBVar>> variables_fa, vector<vector<GRBVar>> variables_ha, vector<vector<GRBVar>> variables_V,
	vector<vector<GRBVar>> variables_ac32, vector<vector<GRBVar>> variables_ac42, vector<vector<GRBVar>> variables_error, vector<double> MED_ac42, 
	vector<double> MED_ac32, GRBModel& model)
{
	for (int i = 1; i < variables_V.size(); i++)
	{
		for (int j = 0; j < variables_V[i].size(); j++)
		{
			if (j == variables_V[i].size() - 1)
			{
				model.addConstr(variables_V[i][j] == 1);
				model.addConstr(variables_fa[i - 1][j] == 0);
				model.addConstr(variables_ha[i - 1][j] == 0);
				model.addConstr(variables_ac32[i - 1][j] == 0);
				model.addConstr(variables_ac42[i - 1][j] == 0);
				model.addConstr(variables_error[i - 1][j] == 0);
			}
			else
			{
				model.addConstr(variables_error[i - 1][j] == MED_ac32[i - 1] * variables_ac32[i - 1][j] + MED_ac42[i - 1] * variables_ac42[i - 1][j]);
				model.addConstr(3 * variables_fa[i - 1][j] + 2 * variables_ha[i - 1][j] + 3 * variables_ac32[i - 1][j] + 4 * variables_ac42[i - 1][j] <= variables_V[i - 1][j]);
				model.addConstr(variables_V[i][j] == variables_V[i - 1][j] - variables_ha[i - 1][j] - 2 * variables_fa[i - 1][j] + (variables_ha[i - 1][j + 1] + variables_fa[i - 1][j + 1]) - variables_ac32[i - 1][j] - 2 * variables_ac42[i - 1][j]);
			}


		}
	}
}

//void generate_constraints_of_h_f_V(vector<vector<GRBVar>> variables_fa, vector<vector<GRBVar>> variables_ha, vector<vector<GRBVar>> variables_V,
//	vector<vector<GRBVar>> variables_ac32,vector<vector<GRBVar>> variables_ac42, vector<vector<GRBVar>> variables_error, GRBModel& model)
//{
//	for (int i = 1; i < variables_V.size(); i++)
//	{
//		for (int j = 0; j < variables_V[i].size(); j++)
//		{
//			if (j == variables_V[i].size() - 1)
//			{
//				model.addConstr(variables_V[i][j] == 1);
//				model.addConstr(variables_fa[i - 1][j] == 0);
//				model.addConstr(variables_ha[i - 1][j] == 0);
//				model.addConstr(variables_ac32[i - 1][j] == 0);
//				model.addConstr(variables_ac42[i - 1][j] == 0);
//				model.addConstr(variables_error[i - 1][j] == 0);
//			}
//			else
//			{
//				model.addConstr(variables_error[i - 1][j] == (1 / double(64)) * variables_ac32[i - 1][j] + (14 / double(256)) * variables_ac42[i - 1][j]);
//				model.addConstr(3 * variables_fa[i - 1][j] + 2 * variables_ha[i - 1][j] + 3 * variables_ac32[i - 1][j] + 4 * variables_ac42[i - 1][j] <= variables_V[i - 1][j]);
//				model.addConstr(variables_V[i][j] == variables_V[i - 1][j] - variables_ha[i - 1][j] - 2 * variables_fa[i - 1][j] + (variables_ha[i - 1][j + 1] + variables_fa[i - 1][j + 1]) - variables_ac32[i - 1][j] - 2 * variables_ac42[i - 1][j]);
//			}
//			
//
//		}
//	}
//}

void generate_constraints_of_last_stage(vector<vector<GRBVar>> variables_V, GRBModel& model)
{
	for (int i = 0; i < variables_V[0].size() - 1; i++)
	{
		model.addConstr(variables_V[stages_num][i] <= 2);
		if (i >= stages_num)
		{
			model.addConstr(variables_V[stages_num][i] >= 1);
		}
	}
	//model.addConstr(variables_V[stages_num][variables_V[0].size() - 1] == 1);
}

void generate_cost_of_adders(vector<vector<GRBVar>> variables_fa, vector<vector<GRBVar>> variables_ha, GRBVar& fa_num, GRBVar& ha_num,
	vector<vector<GRBVar>> variables_ac32, vector<vector<GRBVar>> variables_ac42, GRBVar& ac32_num, GRBVar& ac42_num, GRBModel& model)
{
	GRBLinExpr fs = 0;
	GRBLinExpr hs = 0;
	GRBLinExpr ac32s = 0;
	GRBLinExpr ac42s = 0;
	for (int i = 0; i < variables_fa.size(); i++)
	{
		for (int j = 0; j < variables_fa[i].size(); j++)
		{
			fs += variables_fa[i][j];
			hs += variables_ha[i][j];
			ac32s += variables_ac32[i][j];
			ac42s += variables_ac42[i][j];
			//obj += (3 * variables_fa[i][j] + 2 * variables_ha[i][j]);
		}
	}
	ha_num = model.addVar(0, INFINITY, 0, GRB_INTEGER, "hs");
	fa_num = model.addVar(0, INFINITY, 0, GRB_INTEGER, "fs");
	ac32_num= model.addVar(0, INFINITY, 0, GRB_INTEGER, "ac32s");
	ac42_num = model.addVar(0, INFINITY, 0, GRB_INTEGER, "ac42s");
	model.addConstr(fs == fa_num);
	model.addConstr(hs == ha_num);
	model.addConstr(ac32s == ac32_num);
	model.addConstr(ac42s == ac42_num);
	//return obj;
}

void generate_total_error(vector<vector<GRBVar>> variables_error, GRBVar& E, GRBModel& model)
{
	vector<GRBVar> Es;
	for (int j = 0; j < variables_error[0].size(); j++)
	{
		Es.push_back(model.addVar(0, INFINITY, 0, GRB_CONTINUOUS, "E_" + to_string(j)));
	}
	for (int j = 0; j < variables_error[0].size(); j++)
	{
		GRBLinExpr columnerror = 0;
		for (int i = 0; i < variables_error.size(); i++)
		{
			columnerror += variables_error[i][j];
		}
		model.addConstr(Es[j] == columnerror);
	}
	E = model.addVar(0, INFINITY, 0, GRB_CONTINUOUS, "E");
	GRBLinExpr sumerror = 0;
	for (int j = 0; j < Es.size(); j++)
	{
		sumerror += pow(2, Es.size() - j - 1) * Es[j];
	}
	model.addConstr(E == sumerror);
}
