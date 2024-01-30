#include "gurobi_c++.h"
#include "Multiplier_Tree_Optimization.h"
#include "ConnectionOrderOpt.h"
#include "wireConnect.h"
#include <iostream>
#include <vector>
#include "Simulator.h"
#include <sstream>
#include <fstream>
#include <string>
#include "WriteJson.h"
using namespace std;

vector<vector<double>> allcoefs;
int EP_Approx;

vector<int> input_patterns;
int PARAMETER_W, PARAMETER_L, MULT_SIZE, stages_num;
//0->Constrain error & Optimize area
//1->Constrain area & Optimize error
int Opt_Mode;
float MED_bound, AREA_bound;
float Weight;
string filename;
int simulator_approach;
//vector<int> para_Ms = { 10,16,20,32,64 };
//vector<int> para_Ws = { 5,8,/*10*/13 };
//vector<int> para_Ls = { 10,9,8,7 };
//vector<int> para_Ss = { 5,6,7,8,10 };

int main()
{
	//test
	max_error_AC32_propagation({ 0.21,0.45,0.55 });
	
	//initialize error propagation parameters
	EP_Approx = 2;
	allcoefs.resize(4);
	Init_EP_Coefs(allcoefs[0], "model_0_0.sol");
	Init_EP_Coefs(allcoefs[1], "model_0_1.sol");
	Init_EP_Coefs(allcoefs[2], "model_0_2.sol");
	Init_EP_Coefs(allcoefs[3], "model_1_0.sol");
	//initial parameter m, stages_num
	MULT_SIZE = 8;
	stages_num = 2;

	//initial parameter optimization mode : 
	Opt_Mode = 0;
	MED_bound = 115;
	AREA_bound = 150;
	Weight = 0.01;
	vector<float> bounds;
	//define runtime upper bound
	int Time_Bound_s = 1200;
	//update the file path
	filename = "gurobi";
	//initialize simulation approach : 0->MyHDL, 1->Simulator.out
	simulator_approach = 1;

	bounds = { 115 };
	for (auto bound : bounds)
	{
		if (!Opt_Mode)
			MED_bound = bound;
		else if (Opt_Mode == 1)
			AREA_bound = bound;
		else
			Weight = bound;
		int nFA, nHA, nAC42, nAC32;
		float Opt1_Analyzed_Error, Opt2_Analyzed_Error;
		//optimization process
		try {
			/*******Set up environment******/
			GRBEnv env = GRBEnv();
			GRBModel model = GRBModel(env);
			model.set(GRB_StringParam_LogFile, filename + ".log");
			//model.set(GRB_INT_PAR_MIPFOCUS, "1");
			model.set(GRB_INT_PAR_PRESOLVE, "2");
			/*model.set(GRB_DBL_PAR_TIMELIMIT, to_string(Time_Bound_s));*/
			//model.set(GRB_DBL_PAR_HEURISTICS, "0.5");
			//model.set(GRB_DBL_PAR_MIPGAP, "0.005");
			/*******generate input patterns******/
			input_patterns = {};
			generate_input_patterns();

			/*******Creat e variables******/
			vector<double> MED_ac42(stages_num, double(14) / 256);
			vector<double> MED_ac32(stages_num, double(1) / 64);
			/*for (int i = 0; i < stages_num; i++)
			{
				MED_ac32[i] = pow(1.5, i) * MED_ac32[i];
				MED_ac42[i] = pow(2, i) * MED_ac42[i];
			}*/
			//generate_MEDs(MED_ac42, MED_ac32);

			/*******Create variables******/
			//Compressor tree
			vector<vector<GRBVar>> variables_fa, variables_ha;
			vector<vector<GRBVar>> variables_V;
			vector<vector<GRBVar>> variables_ac32, variables_ac42;
			vector<vector<GRBVar>> variables_error;
			GRBVar fa_num, ha_num, ac32_num, ac42_num, E;
			generate_variables_of_multiplier(variables_fa, variables_ha, variables_V, variables_ac32, variables_ac42, variables_error, model);

			/*Add general constraints*/
			//Compressor Tree
			Initial_constraints_of_V(variables_V, model);
			generate_constraints_of_h_f_V(variables_fa, variables_ha, variables_V, variables_ac32, variables_ac42, variables_error, MED_ac42, MED_ac32, model);
			generate_constraints_of_last_stage(variables_V, model);

			// Set objective
			GRBLinExpr obj = 0;
			generate_cost_of_adders(variables_fa, variables_ha, fa_num, ha_num, variables_ac32, variables_ac42, ac32_num, ac42_num, model);
			generate_total_error(variables_error, E, model);
			// mode 0: constrain error & minimize area
			if (!Opt_Mode)
			{
				model.addConstr(E <= MED_bound);
				obj += (5.05 * fa_num + 2.66 * ha_num + 2.66 * ac32_num + 4.78 * ac42_num);
				//obj += (5 * fa_num + 2 * ha_num + 3 * ac32_num + 5 * ac42_num);
			}
			//mode 1: constrain area & minimize error
			else if (Opt_Mode == 1)
			{
				model.addConstr((5.05 * fa_num + 2.66 * ha_num + 2.66 * ac32_num + 4.78 * ac42_num) <= AREA_bound);
				obj += E;
			}
			else
			{
				obj += (10 * fa_num + 5 * ha_num + 5 * ac32_num + 9 * ac42_num) + Weight * E;
			}
			model.setObjective(obj, GRB_MINIMIZE);

			// Save problem
			model.write(filename + ".lp");

			// Optimize model
			model.update();

			model.optimize();

			//Save solution

			model.write(filename + ".sol");

			cout << fa_num.get(GRB_StringAttr_VarName) << " "
				<< fa_num.get(GRB_DoubleAttr_X) << endl;
			nFA = fa_num.get(GRB_DoubleAttr_X);
			cout << ha_num.get(GRB_StringAttr_VarName) << " "
				<< ha_num.get(GRB_DoubleAttr_X) << endl;
			nHA = ha_num.get(GRB_DoubleAttr_X);
			cout << ac32_num.get(GRB_StringAttr_VarName) << " "
				<< ac32_num.get(GRB_DoubleAttr_X) << endl;
			nAC32 = ac32_num.get(GRB_DoubleAttr_X);
			cout << ac42_num.get(GRB_StringAttr_VarName) << " "
				<< ac42_num.get(GRB_DoubleAttr_X) << endl;
			nAC42 = int(ac42_num.get(GRB_DoubleAttr_X));
			cout << E.get(GRB_StringAttr_VarName) << " "
				<< E.get(GRB_DoubleAttr_X) << endl;
			cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
			Opt1_Analyzed_Error = E.get(GRB_DoubleAttr_X);
			/*repeatly*/
			vector<int> last_stage;
			for (int i = 0; i < variables_V[variables_V.size() - 1].size(); i++)
			{
				cout << variables_V[variables_V.size() - 1][i].get(GRB_DoubleAttr_X) << ", ";
			}

			for (int i = variables_V.size() - 1; i < variables_V[variables_V.size() - 1].size(); i++)
			{
				last_stage.push_back(variables_V[variables_V.size() - 1][i].get(GRB_DoubleAttr_X));
			}

			cout << endl << endl;

			//post verification
			int err_times = 0;
			int test_times = 0;
			cout << "Error distribution of approximate 4:2 compressor:" << endl;
			max_error_AC42_propagation({ 0.25,0.25,0.25,0.25 });
			cout << endl;
			cout << "Error distribution of approximate 3:2 compressor:" << endl;
			max_error_AC32_propagation({ 0.25,0.25,0.25 });
			cout << endl;
			cout << "Distribution of 3:2 compressor:" << endl;
			sim_f();
			cout << endl;
			cout << "Distribution of 2:2 compressor:" << endl;
			sim_h();
			cout << endl;
			int err_sum = 0;
			/*string err_filename = "heuristic_error.txt";
			ofstream f1;
			f1.open(err_filename.c_str());
			for (int i = pow(2, MULT_SIZE) - 1; i >= 0; i--)
			{
				for (int j = 0; j < pow(2, MULT_SIZE); j++)
				{
					int appro_mult = final_mult(i, j);
					int exact_mult = i * j;
					int error_distance = abs(appro_mult - exact_mult);
					err_sum += error_distance;
					test_times++;
					if (error_distance > 0)
					{
						err_times++;
						cout << "This is wrong and error distance is " << error_distance << "!" << endl;
						cout << "The total wrong time is " << err_times << endl;
						cout << "The error rate is " << double(err_times) / test_times << endl;

						f1 << "This is wrong and error distance is " << error_distance << "!" << endl;
						f1 << "The total wrong time is " << err_times << endl;
						f1 << "The error rate is " << double(err_times) / test_times << endl;
					}

				}
			}
			cout << endl;
			cout << "The error rate is " << double(err_times) / test_times << endl;
			cout << "The MED is " << double(err_sum) / test_times << endl;

			f1 << endl;
			f1 << "The error rate is " << double(err_times) / test_times << endl;
			f1 << "The MED is " << double(err_sum) / test_times << endl;
			f1.close();*/
		}
		catch (GRBException e) {
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Exception during optimization" << endl;
		}

		// optimize
		string filename1 = "gurobi_CO";
		try {
			/*******Set up environment******/
			GRBEnv env = GRBEnv();
			GRBModel model_CO = GRBModel(env);
			model_CO.set(GRB_StringParam_LogFile, filename + ".log");
			model_CO.set(GRB_INT_PAR_MIPFOCUS, "1");
			model_CO.set(GRB_INT_PAR_NUMERICFOCUS, "1");
			//model_CO.set(GRB_INT_PAR_CUTS, "2");
			//model_CO.set(GRB_INT_PAR_MIQCPMETHOD, "0");
			//model_CO.set(GRB_INT_PAR_NODEMETHOD, "1");
			//model_CO.set(GRB_INT_PAR_PREQLINEARIZE, "1");
			model_CO.set(GRB_INT_PAR_PRESOLVE, "2");
			model_CO.set(GRB_DBL_PAR_TIMELIMIT, to_string(Time_Bound_s));
			model_CO.set(GRB_DBL_PAR_HEURISTICS, "0.4");
			//model_CO.set(GRB_DBL_PAR_MIPGAP, "0.005");
			model_CO.set(GRB_INT_PAR_NONCONVEX, "2");

			vector<vector<Vars_Column>> allVars;
			vector<vector<GRBVar>> allMEDs;
			vector<int> initial_stage = input_patterns;
			vector<vector<int>> F(stages_num, vector<int>(2 * MULT_SIZE + stages_num - 1, 0));
			vector<vector<int>> H(stages_num, vector<int>(2 * MULT_SIZE + stages_num - 1, 0));
			vector<vector<int>> AC32(stages_num, vector<int>(2 * MULT_SIZE + stages_num - 1, 0));
			vector<vector<int>> AC42(stages_num, vector<int>(2 * MULT_SIZE + stages_num - 1, 0));
			vector<vector<double>> E(stages_num, vector<double>(2 * MULT_SIZE + stages_num - 1, 0.0));
			vector<vector<int>> V(stages_num, vector<int>(2 * MULT_SIZE + stages_num - 1, 0));
			ConstructCT(F, H, AC32, AC42, E, V);
			initial_stage.insert(initial_stage.begin(), stages_num, 0);
			V.insert(V.begin(), initial_stage);
			F.insert(F.end(), vector<int>(F[0].size(), 0));
			H.insert(H.end(), vector<int>(H[0].size(), 0));
			AC32.insert(AC32.end(), vector<int>(AC32[0].size(), 0));
			AC42.insert(AC42.end(), vector<int>(AC42[0].size(), 0));
			Sols sol(V, F, H, AC42, AC32);

			generate_variables_of_ConnectOrderOpt(allVars, sol, model_CO);
			Initial_probs_of_V(allVars[0], model_CO);
			allMEDs = generate_constraints_of_ConnectOrderOpt(allVars, EP_Approx, model_CO);
			GRBLinExpr MED = 0;
			for (int i = 0; i < allMEDs.size(); i++)
				for (int j = 0; j < allMEDs[i].size(); j++)
					MED += allMEDs[i][j] * pow(2, allMEDs[i].size() - j - 1);
			model_CO.setObjective(MED, GRB_MINIMIZE);

			// Save problem
			model_CO.write(filename1 + ".lp");

			// Optimize model

			model_CO.optimize();

			//Save solution

			model_CO.write(filename1 + ".sol");

			Opt2_Analyzed_Error = model_CO.get(GRB_DoubleAttr_ObjVal);
			// post verification
			vector<vector<column_sol_t>> allSols;
			allSols = Init_Sols(allVars);
			string path = "D:\\Microsoft VS Studio 2019\\Project\\ApproxMULT_MyHDL\\ApproxMULT_MyHDL\\";
			string file = (!Opt_Mode) ? "CompressorTree_CEOA" + to_string(int(MED_bound)) + ".json"
				: ((Opt_Mode == 1) ? ("CompressorTree_CAOE" + to_string(int(AREA_bound)) + ".json") :
					("CompressorTree_CoOpt" + to_string(Weight) + ".json"));
			WriteToJson(allSols, path + file);
			string err_filename;
			if (!Opt_Mode)
			{
				if (simulator_approach == 0)
					err_filename = "CEOA\\" + to_string(MULT_SIZE) + "\\" + to_string(int(MED_bound)) + "\\opt_error_" + to_string(MED_bound) + ".txt";
				else if (simulator_approach == 1)
					err_filename = "CEOA\\" + to_string(MULT_SIZE) + "\\" + to_string(int(MED_bound)) + "\\simulator_opt_error_" + to_string(MED_bound) + ".txt";
				system(("mkdir CEOA\\" + to_string(MULT_SIZE) + "\\" + to_string(int(MED_bound))).c_str());
			}
			else if (Opt_Mode == 1)
			{
				if (simulator_approach == 0)
					err_filename = "CAOE\\" + to_string(MULT_SIZE) + "\\" + to_string(int(AREA_bound)) + "\\opt_error_" + to_string(AREA_bound) + ".txt";
				else if (simulator_approach == 1)
					err_filename = "CAOE\\" + to_string(MULT_SIZE) + "\\" + to_string(int(AREA_bound)) + "\\simulator_opt_error_" + to_string(AREA_bound) + ".txt";
				system(("mkdir CAOE\\" + to_string(MULT_SIZE) + "\\" + to_string(int(AREA_bound))).c_str());
			}
			else
			{
				if (simulator_approach == 0)
					err_filename = "CoOpt\\" + to_string(MULT_SIZE) + "\\" + to_string(Weight) + "\\opt_error_" + to_string(Weight) + ".txt";
				else if (simulator_approach == 1)
					err_filename = "CoOpt\\" + to_string(MULT_SIZE) + "\\" + to_string(Weight) + "\\simulator_opt_error_" + to_string(Weight) + ".txt";
				system(("mkdir CoOpt\\" + to_string(MULT_SIZE) + "\\" + to_string(Weight)).c_str());
			}
			ofstream f2;
			f2.open(err_filename.c_str());
			if (!Opt_Mode)
				f2 << "Error bound : " << MED_bound << endl;
			else if (Opt_Mode == 1)
				f2 << "Area bound : " << AREA_bound << endl;
			else
				f2 << "Weight : " << Weight << endl;
			f2 << "FA number : " << nFA << "\t" << "HA number : " << nHA << "\t" << "AC32 number : " << nAC32 << "\t" << "AC42 number : " << nAC42 << endl;
			f2 << "Analytical error of opt1 : " << Opt1_Analyzed_Error << "\t" << "Analytical error of opt2 : " << Opt2_Analyzed_Error << endl;
			f2 << "============================================================" << endl;
			f2.close();
			if (!Opt_Mode)
			{
				if (simulator_approach == 0)
					err_filename = "opt_error_" + to_string(MED_bound) + ".txt";
				else if (simulator_approach == 1)
					err_filename = "simulator_opt_error_" + to_string(MED_bound) + ".txt";
				//system(("mkdir CEOA\\" + to_string(MULT_SIZE) + "\\" + to_string(int(MED_bound))).c_str());
				system(("cd CEOA\\" + to_string(MULT_SIZE) + "\\" + to_string(int(MED_bound)) + "& python \"" + path + "ApproxMULT_MyHDL.py\" \"" + path + file + "\" ILP_ApproxMult " + err_filename + " " + to_string(MULT_SIZE) + " " + to_string(simulator_approach)).c_str());
				
				if (simulator_approach == 1)
				{	
					string command_dos = "cd CEOA\\" + to_string(MULT_SIZE) + "\\" + to_string(int(MED_bound)) + "& ";
					string commmand_yosys = "yosys.exe -p \"read_verilog Multiplier_KoggeStone.v; synth; write_blif Multiplier_KoggeStone.blif;\"";
					system((command_dos + commmand_yosys).c_str());
					/*string exact_mult_filepath = "\/mnt\/d\/Microsoft\ VS\ Studio\ 2019\/Project\/Global_Optimization_of_Mutliplier\/ApproximateMult\/ExactMult\/" + to_string(MULT_SIZE) + ".blif";
					string approximate_mult_filepath = "\/mnt\/d\/Microsoft\ VS\ Studio\ 2019\/Project\/Global_Optimization_of_Mutliplier\/ApproximateMult\/CEOA\/" + to_string(MULT_SIZE) + "\/" + to_string(int(MED_bound)) + "\/Multiplier_KoggeStone.blif";
					string error_report_filepath = "\/mnt\/d\/Microsoft\ VS\ Studio\ 2019\/Project\/Global_Optimization_of_Mutliplier\/ApproximateMult\/CEOA\/" + to_string(MULT_SIZE) + "\/" + to_string(int(MED_bound)) + "\/" + err_filename;
					string command_linux = "wsl.exe \"simulation.out -a " + exact_mult_filepath + " -b " + approximate_mult_filepath + " -f " + error_report_filepath + "; exit;\"";
					system((command_dos + command_linux).c_str());*/
				}
			}
			else if (Opt_Mode == 1)
			{
				if (!Opt_Mode)
					err_filename = "opt_error_" + to_string(AREA_bound) + ".txt";
				else if (simulator_approach == 1)
					err_filename = "simulator_opt_error_" + to_string(AREA_bound) + ".txt";
				//system(("mkdir CAOE\\" + to_string(MULT_SIZE) + "\\" + to_string(int(AREA_bound))).c_str());
				//system("python --verision");
				system(("cd CAOE\\" + to_string(MULT_SIZE) + "\\" + to_string(int(AREA_bound)) + "& python \"" + path + "ApproxMULT_MyHDL.py\" \"" + path + file + "\" ILP_ApproxMult " + err_filename + " " + to_string(MULT_SIZE) + " " + to_string(simulator_approach)).c_str());
				if (simulator_approach == 1)
				{
					string command_dos = "cd CAOE\\" + to_string(MULT_SIZE) + "\\" + to_string(int(AREA_bound)) + "& ";
					string commmand_yosys = "yosys.exe -p \"read_verilog Multiplier_KoggeStone.v; synth; write_blif Multiplier_KoggeStone.blif;\"";
					system((command_dos + commmand_yosys).c_str());
					/*string command_linux = "wsl.exe -e \"simulation.out -a " + exact_mult_filepath + " -b " + "Multiplier_KoggeStone.blif -f " + err_filename + "; exit;\"";
					system((command_dos + command_linux).c_str());*/
				}
			}
			else
			{

				if (!Opt_Mode)
					err_filename = "opt_error_" + to_string(Weight) + ".txt";
				else if (simulator_approach == 1)
					err_filename = "simulator_opt_error_" + to_string(Weight) + ".txt";
				//system(("mkdir CoOpt\\" + to_string(MULT_SIZE) + "\\" + to_string(Weight)).c_str());
				system(("cd CoOpt\\" + to_string(MULT_SIZE) + "\\" + to_string(Weight) + "& python \"" + path + "ApproxMULT_MyHDL.py\" \"" + path + file + "\" ILP_ApproxMult " + err_filename + " " + to_string(MULT_SIZE) + " " + to_string(simulator_approach)).c_str());
				if (simulator_approach == 1)
				{
					string command_dos = "cd CoOpt\\" + to_string(MULT_SIZE) + "\\" + to_string(Weight) + "& ";
					string commmand_yosys = "yosys.exe -p \"read_verilog Multiplier_KoggeStone.v; synth; write_blif Multiplier_KoggeStone.blif;\"";
					system((command_dos + commmand_yosys).c_str());
					/*string command_linux = "wsl.exe \"simulation.out -a " + exact_mult_filepath + " -b " + "Multiplier_KoggeStone.blif -f " + err_filename + "; exit;\"";
					system((command_dos + command_linux).c_str());*/
				}
			}
			//int test_times = 0, err_times = 0;
			//long err_sum = 0;
			//for (int i = pow(2, MULT_SIZE) - 1; i >= 0; i--)
			//{
			//	for (int j = 0; j < pow(2, MULT_SIZE); j++)
			//	{
			//		int appro_mult = final_mult(i, j, allSols);
			//		int exact_mult = i * j;
			//		int error_distance = abs(appro_mult - exact_mult);
			//		err_sum += error_distance;
			//		test_times++;
			//		/*cout << "Inputs: " << i << ", " << j << endl;
			//		cout << "approximate multiplier's result is " << appro_mult << endl;
			//		cout << "exact multiplier's result is " << exact_mult << endl;*/
			//		if (error_distance > 0)
			//		{
			//			err_times++;
			//			/*cout << "This is wrong and error distance is " << error_distance << "!" << endl;
			//			cout << "The total wrong time is " << err_times << endl;
			//			cout << "The error rate is " << double(err_times) / test_times << endl;*/

			//			f2 << "This is wrong and error distance is " << error_distance << "!" << endl;
			//			f2 << "The total wrong time is " << err_times << endl;
			//			f2 << "The error rate is " << double(err_times) / test_times << endl;
			//		}

			//	} 
			//}
			//cout << endl;
			//cout << "The error rate is " << double(err_times) / test_times << endl;
			//cout << "The MED is " << double(err_sum) / test_times << endl;

			//f2 << endl;
			//f2 << "The error rate is " << double(err_times) / test_times << endl;
			//f2 << "The MED is " << double(err_sum) / test_times << endl;
			//f2.close();
		}
		catch (GRBException e) {
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Exception during optimization" << endl;
		}
	}
	// generate verilog files using .py script
	system("pause");
	return 0;
}