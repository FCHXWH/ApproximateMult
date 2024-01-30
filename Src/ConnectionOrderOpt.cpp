#include "ConnectionOrderOpt.h"
#include <sstream>

void Init_EP_Coefs(vector<double>& coefs, string solname)
{
	ifstream f_coef; string s; stringstream ss;
	f_coef.open(solname.c_str());
	while (getline(f_coef, s))
	{
		if (s[0] != '#')
		{
			auto it = s.find(" ");
			double value;
			string svalue = s.substr(it + 1, s.size() - 1);
			ss << svalue;
			ss >> value; ss.clear(); ss.str("");
			if (value <= 0.001)
				value = 0;
			coefs.push_back(value);
		}
	}
	f_coef.close();
}

//AC42_Input_Probs : pa4_0 pa4_1 pa4_2 pa4_3
//AC42_Input_QuadTerms : probs4[2] * probs4[3]  probs4[0] * probs4[3]  probs4[0] * probs4[1]  
//                       probs4[1] * probs4[3]  probs4[1] * probs4[2]
//AC42_Output_Probs: OUTPUT1 OUTPUT2 MED
vector<GRBQuadExpr> Exact_AC42_Output_Probs(vector<GRBVar> AC42_Input_Probs, vector<GRBVar> AC42_Input_QuadTerms)
{
	vector<GRBQuadExpr> AC42_Output_Probs;
	AC42_Output_Probs.resize(3, 0);

	//OUTPUT1: (pa4_1+pa4_0-1)(1-pa4_3 pa4_2 )+pa4_3 pa4_0(pa4_2 pa4_1+2)+pa4_1 pa4_0 (1-2*pa4_3 )+1
	//OUTPUT2: (pa4_2+pa4_3-1)(1-pa4_1 pa4_0 )+pa4_3 pa4_2 (pa4_1 pa4_0- 1) + 1
	//MED: pa4_0 pa4_1 (pa4_2+pa4_3)+pa4_1 pa4_3 (pa4_0+pa4_2)- 2*pa4_0 pa4_1 pa4_2 pa4_3
	AC42_Output_Probs[0] += (AC42_Input_Probs[1] + AC42_Input_Probs[0] - 1) * (1 - AC42_Input_QuadTerms[0])
		+ AC42_Input_QuadTerms[1] * (AC42_Input_QuadTerms[4] + 2) + AC42_Input_QuadTerms[2] * (1 - 2 * AC42_Input_Probs[3]) + 1;
	AC42_Output_Probs[1] += (AC42_Input_Probs[2] + AC42_Input_Probs[3] - AC42_Input_QuadTerms[0] - 1) * (1 - AC42_Input_QuadTerms[2]) + 1;

	AC42_Output_Probs[2] += AC42_Input_QuadTerms[2] * (AC42_Input_Probs[2] + AC42_Input_Probs[3]) + AC42_Input_QuadTerms[3] * (AC42_Input_Probs[0]
		+ AC42_Input_Probs[2]) - 2 * AC42_Input_QuadTerms[0] * AC42_Input_QuadTerms[2];
	return AC42_Output_Probs;
}

vector<GRBQuadExpr> Approx_AC42_Output_Probs(vector<GRBVar> AC42_Input_Probs, vector<GRBVar> AC42_Input_QuadTerms)
{
	vector<GRBQuadExpr> AC42_Output_Probs;
	AC42_Output_Probs.resize(3, 0);

	//OUTPUT1: 
	//OUTPUT2: 
	//MED: 

	AC42_Output_Probs[0] += (AC42_Input_Probs[1] + AC42_Input_Probs[0] - 1) * (1 - AC42_Input_QuadTerms[0]) + 2 * AC42_Input_QuadTerms[1]
		+ AC42_Input_QuadTerms[2] * (1 - 2 * AC42_Input_Probs[3]) + 1;
	/*AC42_Output_Probs[0] += (AC42_Input_Probs[1] + AC42_Input_Probs[0] - 1) + 2 * AC42_Input_QuadTerms[1]
		+ (1 - 2 * AC42_Input_Probs[3]) + 1;*/
	AC42_Output_Probs[1] += (AC42_Input_Probs[2] + AC42_Input_Probs[3] - 1) * (1 - AC42_Input_QuadTerms[2]) - AC42_Input_QuadTerms[0] + 1;
	//AC42_Output_Probs[1] += (AC42_Input_Probs[2] + AC42_Input_Probs[3] - 1) - AC42_Input_QuadTerms[0] + 1;

	AC42_Output_Probs[2] += AC42_Input_QuadTerms[2] * (AC42_Input_Probs[2] + AC42_Input_Probs[3]) + AC42_Input_QuadTerms[3] * (AC42_Input_Probs[0] + AC42_Input_Probs[2]);
	//AC42_Output_Probs[2] += (AC42_Input_Probs[2] + AC42_Input_Probs[3]) + AC42_Input_QuadTerms[3];
	return AC42_Output_Probs;
}

vector<GRBQuadExpr> Approx_AC42_Output_Probs_2(vector<GRBVar> AC42_Input_Probs, vector<GRBVar> AC42_Input_QuadTerms)
{
	vector<GRBQuadExpr> AC42_Output_Probs;
	AC42_Output_Probs.resize(3, 0);

	//OUTPUT1: 
	//OUTPUT2: 
	//MED: 
	int vars_num = 4;
	for (int i = 0; i < vars_num; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (allcoefs[j][i] != 0)
				AC42_Output_Probs[j] += AC42_Input_Probs[i] * AC42_Input_Probs[i] * allcoefs[j][i];
		}
	}
	int index = vars_num;
	for (int i = 0; i < vars_num - 1; i++)
	{
		for (int j = i + 1; j < vars_num; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				if (allcoefs[k][index] != 0)
					AC42_Output_Probs[k] += AC42_Input_Probs[i] * AC42_Input_Probs[j] * allcoefs[k][index];
			}
			index++;
		}
	}
	for (int i = 0; i < vars_num; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (allcoefs[j][index] != 0)
				AC42_Output_Probs[j] += AC42_Input_Probs[i] * allcoefs[j][index];
		}
		index++;
	}
	for (int i = 0; i < 3; i++)
	{
		if (allcoefs[i][index] != 0)
			AC42_Output_Probs[i] += allcoefs[i][index];
	}
	return AC42_Output_Probs;
}

//AC32_Input_Probs : pa3_0 pa3_1 pa3_2
//AC32_Input_QuadTerms : probs3[0] * probs3[1]
vector<GRBQuadExpr> Exact_AC32_Output_Probs(vector<GRBVar> AC32_Input_Probs, vector<GRBVar> AC32_Input_QuadTerms)
{
	vector<GRBQuadExpr> AC32_Output_Probs;
	AC32_Output_Probs.resize(3, 0);

	//OUTPUT1: pa3_0+pa3_1-pa3_1 pa3_0
	//OUTPUT2: pa3_2+pa3_1 pa3_0-pa3_2 pa3_1 pa3_0
	//MED: pa3_0 pa3_1 pa3_2
	AC32_Output_Probs[0] += AC32_Input_Probs[0] + AC32_Input_Probs[1] - AC32_Input_QuadTerms[0];
	AC32_Output_Probs[1] += AC32_Input_Probs[2] + AC32_Input_QuadTerms[0] * (1 - AC32_Input_Probs[2]);
	AC32_Output_Probs[2] += AC32_Input_QuadTerms[0] * AC32_Input_Probs[2];
	return AC32_Output_Probs;
}

vector<GRBQuadExpr> Approx_AC32_Output_Probs(vector<GRBVar> AC32_Input_Probs, vector<GRBVar> AC32_Input_QuadTerms)
{
	vector<GRBQuadExpr> AC32_Output_Probs;
	AC32_Output_Probs.resize(3, 0);

	//OUTPUT1: 
	//OUTPUT2: 
	//MED: 
	AC32_Output_Probs[0] += AC32_Input_Probs[0] + AC32_Input_Probs[1] - AC32_Input_QuadTerms[0];
	AC32_Output_Probs[1] += AC32_Input_Probs[2] + AC32_Input_QuadTerms[0];
	AC32_Output_Probs[2] += AC32_Input_QuadTerms[0] * AC32_Input_Probs[2];
	//AC32_Output_Probs[2] += 0.1 * (AC32_Input_Probs[1] * AC32_Input_Probs[2] + AC32_Input_Probs[0] * AC32_Input_Probs[2] + AC32_Input_Probs[0] * AC32_Input_Probs[1]);
	return AC32_Output_Probs;
}

vector<GRBQuadExpr> Approx_AC32_Output_Probs_2(vector<GRBVar> AC32_Input_Probs, vector<GRBVar> AC32_Input_QuadTerms)
{
	vector<GRBQuadExpr> AC32_Output_Probs;
	AC32_Output_Probs.resize(3, 0);

	//OUTPUT1: 
	//OUTPUT2: 
	//MED: 
	AC32_Output_Probs[0] += AC32_Input_Probs[0] + AC32_Input_Probs[1] - AC32_Input_QuadTerms[0];
	AC32_Output_Probs[1] += AC32_Input_Probs[2] + AC32_Input_QuadTerms[0];
	int vars_num = 3;
	for (int i = 0; i < vars_num; i++)
	{
		if (allcoefs[3][i] != 0)
			AC32_Output_Probs[2] += AC32_Input_Probs[i] * AC32_Input_Probs[i] * allcoefs[3][i];
	}
	int index = vars_num;
	for (int i = 0; i < vars_num - 1; i++)
	{
		for (int j = i + 1; j < vars_num; j++)
		{
			if (allcoefs[3][index] != 0)
				AC32_Output_Probs[2] += AC32_Input_Probs[i] * AC32_Input_Probs[j] * allcoefs[3][index];
			index++;
		}
	}
	for (int i = 0; i < vars_num; i++)
	{
		if (allcoefs[3][index] != 0)
			AC32_Output_Probs[2] += AC32_Input_Probs[i] * allcoefs[3][index];
		index++;
	}
	AC32_Output_Probs[2] += allcoefs[3][index];
	return AC32_Output_Probs;
}


//F_Input_Probs : pf_0 pf_1 pf_2
//F_Input_QuadTerms : probsf[0] * probsf[1]  probsf[0] * probsf[2]  probsf[1] * probsf[2]
vector<GRBQuadExpr> Exact_F_Output_Probs(vector<GRBVar> F_Input_Probs, vector<GRBVar> F_Input_QuadTerms)
{
	vector<GRBQuadExpr> F_Output_Probs;
	F_Output_Probs.resize(2);

	//OUTPUT1: (pf_0+pf_1+pf_2 )-2(pf_0 pf_1+pf_0 pf_2+pf_1 pf_2 )+4pf_0 pf_1 pf_2
	//OUTPUT2: pf_0 pf_1+pf_0 pf_2+pf_1 pf_2-2*pf_0 pf_1 pf_2
	F_Output_Probs[0] = (F_Input_Probs[0] + F_Input_Probs[1] + F_Input_Probs[2]) - 2 * (F_Input_QuadTerms[0] + F_Input_QuadTerms[1] + F_Input_QuadTerms[2])
		+ 4 * F_Input_QuadTerms[0] * F_Input_Probs[2];
	F_Output_Probs[1] = F_Input_QuadTerms[0] + F_Input_QuadTerms[1] + F_Input_QuadTerms[2] - 2 * F_Input_QuadTerms[0] * F_Input_Probs[2];
	return F_Output_Probs;
}

vector<GRBQuadExpr> Approx_F_Output_Probs(vector<GRBVar> F_Input_Probs, vector<GRBVar> F_Input_QuadTerms)
{
	vector<GRBQuadExpr> F_Output_Probs;
	F_Output_Probs.resize(2);

	//OUTPUT1: 
	//OUTPUT2: 
	F_Output_Probs[0] = (F_Input_Probs[0] + F_Input_Probs[1] + F_Input_Probs[2]) - (F_Input_QuadTerms[0] + F_Input_QuadTerms[1] + F_Input_QuadTerms[2]);
	F_Output_Probs[1] = F_Input_QuadTerms[0] + F_Input_QuadTerms[1] + F_Input_QuadTerms[2];
	return F_Output_Probs;
}

void generate_variables_of_ConnectOrderOpt(vector<vector<Vars_Column>>& allVars, Sols sol, GRBModel& model)
{
	vector<vector<int>> Vs;
	vector<vector<int>> Fs, Hs;
	vector<vector<int>> AC42, AC32;
	Vs = sol.Vs;
	Fs = sol.Fs;
	Hs = sol.Hs;
	AC42 = sol.AC42;
	AC32 = sol.AC32;
	allVars.resize(Vs.size());
	for (int i = 0; i < Vs.size(); i++)
	{
		vector<Vars_Column> allVars_eachstage;
		allVars_eachstage.resize(Vs[0].size());
		for (int j = 0; j < Vs[i].size(); j++)
		{
			Vars_Column col(Vs[i][j], AC42[i][j], AC32[i][j], Fs[i][j], Hs[i][j]);
			if (i > 0)
			{
				if (j == Vs[i].size() - 1)
				{
					col.nfc_ls = 0;
					col.nhc_ls = 0;
				}
				else
				{
					col.nfc_ls = allVars[i - 1][j + 1].nF;
					col.nhc_ls = allVars[i - 1][j + 1].nH;
				}
				col.nac42_ls = allVars[i - 1][j].nac42 * 2;
				col.nac32_ls = allVars[i - 1][j].nac32 * 2;
				col.nfs_ls = allVars[i - 1][j].nF;
				col.nhs_ls = allVars[i - 1][j].nH;
				col.nrb_ls = allVars[i - 1][j].nrb;
			}
			else
			{
				col.nfc_ls = 0;
				col.nhc_ls = 0;
				col.nac42_ls = 0;
				col.nac32_ls = 0;
				col.nfs_ls = 0;
				col.nhs_ls = 0;
				/*col.nrb_ls = col.nrb;*/
				col.nrb_ls = 0;
			}


			vector<vector<GRBVar>>& RemainBits = col.RemainBits;
			vector<vector<GRBVar>>& AC42 = col.AC42;
			vector<vector<GRBVar>>& AC32 = col.AC32;
			vector<vector<GRBVar>>& f = col.f;
			vector<vector<GRBVar>>& h = col.h;
			vector<GRBVar>& Probs = col.Probs;
			for (int k = 0; k < RemainBits.size(); k++)
			{
				stringstream ss; string s;
				ss << "p" << i << "_" << j << "_" << k;
				ss >> s;
				ss.clear(); ss.str("");
				Probs[k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
				for (int m = 0; m < RemainBits[k].size(); m++)
				{
					//stringstream ss; string s;
					ss << "rb" << i << "_" << j << "_" << k << "_" << m;
					ss >> s;
					ss.clear(); ss.str("");
					RemainBits[k][m] = model.addVar(0, 1, 0, GRB_BINARY, s);
				}
				for (int m = 0; m < AC42[k].size(); m++)
				{
					//stringstream ss; string s;
					int iAC4 = m / 4 + 1;
					int iInput = m % 4;
					// i-th stage, j-th column, k-th bit, iAC4-th ac42 compressor, iInput-th input
					ss << "a4c" << i << "_" << j << "_" << k << "_" << iAC4 << "_" << iInput;
					ss >> s;
					ss.clear(); ss.str("");
					AC42[k][m] = model.addVar(0, 1, 0, GRB_BINARY, s);
				}
				for (int m = 0; m < AC32[k].size(); m++)
				{
					//stringstream ss; string s;
					int iAC3 = m / 3 + 1;
					int iInput = m % 3;
					// i-th stage, j-th column, k-th bit, iAC3-th ac32 compressor, iInput-th input
					ss << "a3c" << i << "_" << j << "_" << k << "_" << iAC3 << "_" << iInput;
					ss >> s;
					ss.clear(); ss.str("");
					AC32[k][m] = model.addVar(0, 1, 0, GRB_BINARY, s);
				}
				for (int m = 0; m < f[k].size(); m++)
				{
					//stringstream ss; string s;
					int iF = m / 3 + 1;
					int iInput = m % 3;
					// i-th stage, j-th column, k-th bit, iF-th 3:2 compressor, iInput-th input
					ss << "f" << i << "_" << j << "_" << k << "_" << iF << "_" << iInput;
					ss >> s;
					ss.clear(); ss.str("");
					f[k][m] = model.addVar(0, 1, 0, GRB_BINARY, s);
				}
				for (int m = 0; m < h[k].size(); m++)
				{
					//stringstream ss; string s;
					int iH = m / 2 + 1;
					int iInput = m % 2;
					// i-th stage, j-th column, k-th bit, iH-th 2:2 compressor, iInput-th input
					ss << "h" << i << "_" << j << "_" << k << "_" << iH << "_" << iInput;
					ss >> s;
					ss.clear(); ss.str("");
					h[k][m] = model.addVar(0, 1, 0, GRB_BINARY, s);
				}

			}
			allVars_eachstage[j] = col;
		}
		allVars[i] = allVars_eachstage;
	}
}

void Initial_probs_of_V(vector<Vars_Column> Vars_stage0, GRBModel& model)
{
	for (int i = 0; i < Vars_stage0.size(); i++)
	{
		Vars_Column col = Vars_stage0[i];
		vector<GRBVar> probs = col.Probs;
		for (auto prob : probs)
			model.addConstr(prob == 0.25);
	}
}

vector<vector<GRBVar>> generate_constraints_of_ConnectOrderOpt(vector<vector<Vars_Column>> allVars, int isApprox, GRBModel& model)
{
	vector<vector<GRBVar>> allMEDs;
	allMEDs.resize(allVars.size() - 1);
	for (int i = 0; i < allMEDs.size(); i++)
		allMEDs[i].resize(allVars[i].size());

	for (int i = 0; i < allVars.size(); i++)
	{
		for (int j = 0; j < allVars[i].size(); j++)
		{
			Vars_Column col = allVars[i][j];
			vector<vector<GRBVar>> RemainBits = col.RemainBits;
			vector<vector<GRBVar>> AC42 = col.AC42;
			vector<vector<GRBVar>> AC32 = col.AC32;
			vector<vector<GRBVar>> f = col.f;
			vector<vector<GRBVar>> h = col.h;
			vector<GRBVar> Probs = col.Probs;

			// constraints: all bits' probabilities of the first stage are same 
			// such that all selection variables can be determined directly;
			if (i == 0)
			{
				//order: remaining bits, AC42 bits, AC32 bits, fa bits, ha bits
				int startRbs = 0, endRbs = col.nrb - 1;
				for (int k = 0; k < RemainBits.size(); k++)
				{
					for (int l = 0; l < RemainBits[k].size(); l++)
					{
						if (l == k)
							model.addConstr(RemainBits[k][l] == 1);
						else
							model.addConstr(RemainBits[k][l] == 0);
					}
				}
				int startAC42 = endRbs + 1, endAC42 = startAC42 + col.nac42 * 4 - 1;
				for (int k = 0; k < AC42.size(); k++)
				{
					for (int l = 0; l < AC42[k].size(); l++)
					{
						if (l == k - startAC42)
							model.addConstr(AC42[k][l] == 1);
						else
							model.addConstr(AC42[k][l] == 0);
					}
				}
				int startAC32 = endAC42 + 1, endAC32 = startAC32 + col.nac32 * 3 - 1;
				for (int k = 0; k < AC32.size(); k++)
				{
					for (int l = 0; l < AC32[k].size(); l++)
					{
						if (l == k - startAC32)
							model.addConstr(AC32[k][l] == 1);
						else
							model.addConstr(AC32[k][l] == 0);
					}
				}
				int startF = endAC32 + 1, endF = startF + col.nF * 3 - 1;
				for (int k = 0; k < f.size(); k++)
				{
					for (int l = 0; l < f[k].size(); l++)
					{
						if (l == k - startF)
							model.addConstr(f[k][l] == 1);
						else
							model.addConstr(f[k][l] == 0);
					}
				}
				int startH = endF + 1, endH = startH + col.nH - 1;
				for (int k = 0; k < h.size(); k++)
				{
					for (int l = 0; l < h[k].size(); l++)
					{
						if (l == k - startH)
							model.addConstr(h[k][l] == 1);
						else
							model.addConstr(h[k][l] == 0);
					}
				}
			}
			else
			{
				// constraints about all selection variables
				vector<GRBLinExpr> sums_rb(col.nrb, 0);
				for (int k = 0; k < RemainBits.size(); k++)
				{
					for (int l = 0; l < RemainBits[k].size(); l++)
						sums_rb[l] += RemainBits[k][l];
				}
				for (int k = 0; k < sums_rb.size(); k++)
					model.addConstr(sums_rb[k] == 1);
				vector<GRBLinExpr> sums_a4(col.nac42 * 4, 0);
				for (int k = 0; k < AC42.size(); k++)
				{
					for (int l = 0; l < AC42[k].size(); l++)
						sums_a4[l] += AC42[k][l];
				}
				for (int k = 0; k < sums_a4.size(); k++)
					model.addConstr(sums_a4[k] == 1);
				vector<GRBLinExpr> sums_a3(col.nac32 * 3, 0);
				for (int k = 0; k < AC32.size(); k++)
				{
					for (int l = 0; l < AC32[k].size(); l++)
						sums_a3[l] += AC32[k][l];
				}
				for (int k = 0; k < sums_a3.size(); k++)
					model.addConstr(sums_a3[k] == 1);

				vector<GRBLinExpr> sums_f(col.nF * 3, 0);
				for (int k = 0; k < f.size(); k++)
				{
					for (int l = 0; l < f[k].size(); l++)
						sums_f[l] += f[k][l];
				}
				for (int k = 0; k < sums_f.size(); k++)
					model.addConstr(sums_f[k] == 1);

				vector<GRBLinExpr> sums_h(col.nH * 2, 0);
				for (int k = 0; k < h.size(); k++)
				{
					for (int l = 0; l < h[k].size(); l++)
						sums_h[l] += h[k][l];
				}
				for (int k = 0; k < sums_h.size(); k++)
					model.addConstr(sums_h[k] == 1);
				// constraints: each bit must be allocated in one position
				// i-th stage, j-th column, k-th bit
				vector<GRBLinExpr> sums_each_bit(RemainBits.size(), 0);
				for (int k = 0; k < RemainBits.size(); k++)
				{
					for (int l = 0; l < RemainBits[k].size(); l++)
						sums_each_bit[k] += RemainBits[k][l];
					for (int l = 0; l < AC42[k].size(); l++)
						sums_each_bit[k] += AC42[k][l];
					for (int l = 0; l < AC32[k].size(); l++)
						sums_each_bit[k] += AC32[k][l];
					for (int l = 0; l < f[k].size(); l++)
						sums_each_bit[k] += f[k][l];
					for (int l = 0; l < h[k].size(); l++)
						sums_each_bit[k] += h[k][l];
					model.addConstr(sums_each_bit[k] == 1);
				}

				// constraints: probability propagation
				// remaining bits
				Vars_Column last_col = allVars[i - 1][j];
				for (int k = 0; k < last_col.nrb; k++)
				{
					GRBQuadExpr prob_sum = 0;
					for (int l = 0; l < last_col.RemainBits.size(); l++)
						prob_sum += last_col.Probs[l] * last_col.RemainBits[l][k];
					model.addQConstr(prob_sum == Probs[k]);
				}

				// AC42 output bits
				GRBLinExpr medExpr = 0;
				vector<GRBVar> probs4;
				probs4.resize(last_col.nac42 * 4);
				for (int k = 0; k < probs4.size(); k++)
				{
					GRBQuadExpr prob_sum = 0;
					for (int l = 0; l < last_col.RemainBits.size(); l++)
						prob_sum += last_col.Probs[l] * last_col.AC42[l][k];
					stringstream ss; string s;
					ss << "p4a" << i << "_" << j << "_" << k / 4 + 1 << "_" << k % 4;
					ss >> s;
					ss.clear(); ss.str("");
					probs4[k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(probs4[k] == prob_sum);
				}
				int startAC42 = col.nrb_ls;
				for (int k = 0; k < last_col.nac42; k++)
				{
					stringstream ss; string s;
					// ac42's quadratic term: i-th stage, j-th column, k-th ac42, 0-th term
					ss << "q4t" << i << "_" << j << "_" << k << "_" << 0;
					ss >> s; ss.clear(); ss.str("");
					GRBVar Qterm0 = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(Qterm0 == probs4[2] * probs4[3]);
					ss << "q4t" << i << "_" << j << "_" << k << "_" << 1;
					ss >> s; ss.clear(); ss.str("");
					GRBVar Qterm1 = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(Qterm1 == probs4[0] * probs4[3]);
					ss << "q4t" << i << "_" << j << "_" << k << "_" << 2;
					ss >> s; ss.clear(); ss.str("");
					GRBVar Qterm2 = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(Qterm2 == probs4[0] * probs4[1]);
					ss << "q4t" << i << "_" << j << "_" << k << "_" << 3;
					ss >> s; ss.clear(); ss.str("");
					GRBVar Qterm3 = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(Qterm3 == probs4[1] * probs4[3]);
					ss << "q4t" << i << "_" << j << "_" << k << "_" << 4;
					ss >> s; ss.clear(); ss.str("");
					GRBVar Qterm4 = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(Qterm4 == probs4[1] * probs4[2]);

					vector<GRBQuadExpr> output_probs;
					if (isApprox == 1)
						output_probs = Approx_AC42_Output_Probs(probs4, { Qterm0,Qterm1,Qterm2,Qterm3,Qterm4 });
					else if (isApprox == 2)
						output_probs = Approx_AC42_Output_Probs_2(probs4, { Qterm0,Qterm1,Qterm2,Qterm3,Qterm4 });
					else
						output_probs = Exact_AC42_Output_Probs(probs4, { Qterm0,Qterm1,Qterm2,Qterm3,Qterm4 });
					model.addQConstr(Probs[startAC42 + 2 * k] == output_probs[0]);
					model.addQConstr(Probs[startAC42 + 2 * k + 1] == output_probs[1]);

					// approximate AC42 Output probs
					//model.addQConstr(Probs[startAC42 + 2 * k] == (probs4[1] + probs4[0] - 1) * (1 - Qterm0) + 2 * Qterm1 + Qterm2 * (1 - 2 * probs4[3]) + 1);
					//model.addQConstr(Probs[startAC42 + 2 * k + 1] == (probs4[2] + probs4[3] - 1) * (1 - Qterm2) - Qterm0 + 1);

					ss << "med4a" << i << "_" << j << "_" << k;
					ss >> s; ss.clear(); ss.str("");
					GRBVar med = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					// approximate AC42 MED
					//model.addQConstr(med == Qterm2 * (probs4[2] + probs4[3]) + Qterm3 * (probs4[0] + probs4[2]));
					model.addQConstr(med == output_probs[2]);
					medExpr += med;
				}

				// AC32 output bits
				vector<GRBVar> probs3;
				probs3.resize(last_col.nac32 * 3);
				for (int k = 0; k < probs3.size(); k++)
				{
					GRBQuadExpr prob_sum = 0;
					for (int l = 0; l < last_col.RemainBits.size(); l++)
						prob_sum += last_col.Probs[l] * last_col.AC32[l][k];
					stringstream ss; string s;
					ss << "p3a" << i << "_" << j << "_" << k / 3 + 1 << "_" << k % 3;
					ss >> s;
					ss.clear(); ss.str("");
					probs3[k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(probs3[k] == prob_sum);
				}
				int StartAC32 = startAC42 + col.nac42_ls;
				for (int k = 0; k < last_col.nac32; k++)
				{
					stringstream ss; string s;
					// ac32's quadratic term: i-th stage, j-th column, k-th ac32, 0-th term
					ss << "q3t" << i << "_" << j << "_" << k << "_" << 0;
					ss >> s; ss.clear(); ss.str("");
					GRBVar Qterm0 = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(Qterm0 == probs3[0] * probs3[1]);

					vector<GRBQuadExpr> output_probs;
					if (isApprox == 1)
						output_probs = Approx_AC32_Output_Probs(probs3, { Qterm0 });
					else if (isApprox == 2)
						output_probs = Approx_AC32_Output_Probs_2(probs3, { Qterm0 });
					else
						output_probs = Exact_AC32_Output_Probs(probs3, { Qterm0 });
					model.addQConstr(Probs[StartAC32 + 2 * k] == output_probs[0]);
					model.addQConstr(Probs[StartAC32 + 2 * k + 1] == output_probs[1]);
					// approximate AC32 output probs
					//model.addConstr(Probs[StartAC32 + 2 * k] == probs3[0] + probs3[1] - Qterm0);
					//model.addConstr(Probs[StartAC32 + 2 * k + 1] == probs3[2] + Qterm0);	

					ss << "med3a" << i << "_" << j << "_" << k;
					ss >> s; ss.clear(); ss.str("");
					GRBVar med = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(med == output_probs[2]);
					// approximate AC32 MED
					//model.addQConstr(med == Qterm0 * probs3[2]);
					medExpr += med;
				}
				stringstream ss; string s;
				ss << "medcol" << i << "_" << j;
				ss >> s; ss.clear(); ss.str("");
				GRBVar medcol = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
				model.addConstr(medExpr == medcol);
				allMEDs[i - 1][j] = medcol;

				//FA output bits
				vector<GRBVar> probsf;
				probsf.resize(last_col.nF * 3);
				for (int k = 0; k < probsf.size(); k++)
				{
					GRBQuadExpr prob_sum = 0;
					for (int l = 0; l < last_col.RemainBits.size(); l++)
						prob_sum += last_col.Probs[l] * last_col.f[l][k];
					stringstream ss; string s;
					ss << "pf" << i << "_" << j << "_" << k / 3 + 1 << "_" << k % 3;
					ss >> s;
					ss.clear(); ss.str("");
					probsf[k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(probsf[k] == prob_sum);
				}
				int Startf = StartAC32 + col.nac32_ls;
				int Startfc;
				if (j > 0)
				{
					Startfc = allVars[i][j - 1].nrb_ls + allVars[i][j - 1].nac42_ls + allVars[i][j - 1].nac32_ls
						+ allVars[i][j - 1].nfs_ls + allVars[i][j - 1].nhs_ls;
				}
				for (int k = 0; k < last_col.nF; k++)
				{
					stringstream ss; string s;
					// f's quadratic term: i-th stage, j-th column, k-th f, 0-th term
					ss << "qft" << i << "_" << j << "_" << k << "_" << 0;
					ss >> s; ss.clear(); ss.str("");
					GRBVar Qterm0 = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(Qterm0 == probsf[0] * probsf[1]);

					ss << "qft" << i << "_" << j << "_" << k << "_" << 1;
					ss >> s; ss.clear(); ss.str("");
					GRBVar Qterm1 = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(Qterm1 == probsf[0] * probsf[2]);

					ss << "qft" << i << "_" << j << "_" << k << "_" << 2;
					ss >> s; ss.clear(); ss.str("");
					GRBVar Qterm2 = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(Qterm2 == probsf[1] * probsf[2]);

					vector<GRBQuadExpr> output_probs;
					if (isApprox)
						output_probs = Approx_F_Output_Probs(probsf, { Qterm0,Qterm1,Qterm2 });
					else
						output_probs = Exact_F_Output_Probs(probsf, { Qterm0,Qterm1,Qterm2 });
					model.addQConstr(Probs[Startf + k] == output_probs[0]);
					model.addQConstr(allVars[i][j - 1].Probs[Startfc + k] == output_probs[1]);
					//approximate FA output probs
					//model.addConstr(Probs[Startf + k] == (probsf[0] + probsf[1] + probsf[2]) - (Qterm0 + Qterm1 + Qterm2));
					//model.addConstr(allVars[i][j - 1].Probs[Startfc + k] == Qterm0 + Qterm1 + Qterm2);
				}

				//HA output bits
				vector<GRBVar> probsh;
				probsh.resize(last_col.nH * 2);
				for (int k = 0; k < probsh.size(); k++)
				{
					GRBQuadExpr prob_sum = 0;
					for (int l = 0; l < last_col.RemainBits.size(); l++)
						prob_sum += last_col.Probs[l] * last_col.h[l][k];
					stringstream ss; string s;
					ss << "ph" << i << "_" << j << "_" << k / 2 + 1 << "_" << k % 2;
					ss >> s;
					ss.clear(); ss.str("");
					probsh[k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, s);
					model.addQConstr(probsh[k] == prob_sum);
				}
				int Starth = Startf + col.nfs_ls;
				int Starthc;
				if (j > 0)
				{
					Starthc = Startfc + last_col.nF;
				}
				for (int k = 0; k < last_col.nH; k++)
				{
					model.addQConstr(Probs[Starth + k] == probsh[0] + probsh[1] - 2 * probsh[0] * probsh[1]);
					model.addQConstr(allVars[i][j - 1].Probs[Starthc + k] == probsh[0] * probsh[1]);
				}
			}

		}
	}
	return allMEDs;
}