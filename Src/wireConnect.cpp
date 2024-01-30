#include "wireConnect.h"
#include <algorithm>

using namespace std;

void column_sol_t::read_GRB_Solu(Vars_Column col, Vars_Column col_left)
{
	vector<vector<GRBVar>> RemainBits = col.RemainBits;
	vector<vector<GRBVar>> AC42 = col.AC42;
	vector<vector<GRBVar>> AC32 = col.AC32;
	vector<vector<GRBVar>> f = col.f;
	vector<vector<GRBVar>> h = col.h;
	

	for (int i = 0; i < col.nrb; i++)
	{
		for (int j = 0; j < RemainBits.size(); j++)
		{
			if (RemainBits[j][i].get(GRB_DoubleAttr_X))
			{
				col_iRB_sel[i] = j;
				col_oRB_pos[i] = i;
				break;
			}
		}
	}

	for (int i = 0; i < col.nac42 * 4; i++)
	{
		for (int j = 0; j < AC42.size(); j++)
		{
			if (AC42[j][i].get(GRB_DoubleAttr_X))
			{
				int iAC42Var = i % 4;
				int iAC42 = i / 4;
				col_iAC42_sel[iAC42][iAC42Var] = j;
			}
		}
	}
	int oStartAC42 = col.nrb, oEndAC42 = oStartAC42 + col.nac42 * 2 - 1;
	for (int i = 0; i < col.nac42; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			int k = oStartAC42 + 2 * i + j;
			col_oAC42_pos[i][j] = k;
		}
	}

	for (int i = 0; i < col.nac32 * 3; i++)
	{
		for (int j = 0; j < AC32.size(); j++)
		{
			if (AC32[j][i].get(GRB_DoubleAttr_X))
			{
				int iAC32Var = i % 3;
				int iAC32 = i / 3;
				col_iAC32_sel[iAC32][iAC32Var] = j;
			}
		}
	}
	int oStartAC32 = oEndAC42 + 1, oEndAC32 = oStartAC32 + col.nac32 * 2 - 1;
	for (int i = 0; i < col.nac32; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			int k = oStartAC32 + 2 * i + j;
			col_oAC32_pos[i][j] = k;
		}
	}

	for (int i = 0; i < col.nF * 3; i++)
	{
		for (int j = 0; j < f.size(); j++)
		{
			if (f[j][i].get(GRB_DoubleAttr_X))
			{
				int iFVar = i % 3;
				int iF = i / 3;
				col_iF_sel[iF][iFVar] = j;
			}
		}
	}
	int oStartF = oEndAC32 + 1, oEndF = oStartF + col.nF - 1;
	int oStartF_left = col_left.nrb + col_left.nac42 * 2 + col_left.nac32 * 2 + col_left.nF + col_left.nH;
	int oEndF_left = oStartF_left + col.nF - 1;
	for (int i = 0; i < col.nF; i++)
	{
		col_oF_pos[i][0] = oStartF + i;
		col_oF_pos[i][1] = oStartF_left + i;
	}

	for (int i = 0; i < col.nH * 2; i++)
	{
		for (int j = 0; j < h.size(); j++)
		{
			if (h[j][i].get(GRB_DoubleAttr_X))
			{
				int iHVar = i % 2;
				int iH = i / 2;
				col_iH_sel[iH][iHVar] = j;
			}
		}
	}
	int oStartH = oEndF + 1, oEndH = oStartH + col.nH - 1;
	int oStartH_left = oEndF_left + 1;
	int oEndH_left = oStartH_left + col.nH - 1;
	for (int i = 0; i < col.nH; i++)
	{
		// sum, cout
		col_oH_pos[i][0] = oStartH + i;
		col_oH_pos[i][1] = oStartH_left + i;
	}
}

vector<int> AC32_Op(vector<int> inputs)
{
	vector<int> outputs;
	outputs.resize(2, 0);
	//outputs[0] = inputs[0] * inputs[1] || inputs[2];
	//outputs[1] = inputs[0] || inputs[1];
	outputs[0] = inputs[0] || inputs[1];
	outputs[1] = inputs[0] * inputs[1] || inputs[2];
	return outputs;
}

vector<int> AC42_Op(vector<int> inputs)
{
	vector<int> outputs;
	outputs.resize(2, 0);
	//outputs[0] = inputs[0] * inputs[1] || inputs[2] || inputs[3];
	//outputs[1] = inputs[2] * inputs[3] || inputs[0] || inputs[1];
	outputs[0] = inputs[2] * inputs[3] || inputs[0] || inputs[1];
	outputs[1] = inputs[0] * inputs[1] || inputs[2] || inputs[3];
	return outputs;
}

vector<int> F_Op(vector<int> inputs)
{
	vector<int> outputs;
	outputs.resize(2, 0);
	// sum cout
	outputs[0] = inputs[0] ^ inputs[1] ^ inputs[2];
	outputs[1] = inputs[0] * inputs[1] || inputs[0] * inputs[2] || inputs[1] * inputs[2];
	return outputs;
}

vector<int> H_Op(vector<int> inputs)
{
	vector<int> outputs;
	outputs.resize(2, 0);
	// sum cout
	outputs[0] = inputs[0] ^ inputs[1];
	outputs[1] = inputs[0] * inputs[1];
	return outputs;
}

vector<vector<column_sol_t>> Init_Sols(vector<vector<Vars_Column>> allVars)
{
	vector<vector<column_sol_t>> allSols;
	// initialize allSols from allVars
	allSols.resize(stages_num + 1, vector<column_sol_t>(allVars[0].size(), column_sol_t()));
	for (int i = 0; i < allVars.size(); i++)
	{
		for (int j = 0; j < allVars[i].size(); j++)
		{
			allSols[i][j] = column_sol_t(allVars[i][j]);
			if (j > 0)
				allSols[i][j].read_GRB_Solu(allVars[i][j], allVars[i][j - 1]);
		}

	}
	return allSols;
}

vector<vector<int>> OptCompressTree(vector<vector<int>> stageBits, vector<vector<column_sol_t>> allSols)
{
	// do compressing
	vector<vector<int>> nextStageBits;
	for (int i = 0; i < allSols.size() - 1; i++)
	{
		nextStageBits.clear();
		for (int j = 0; j < allSols[i + 1].size(); j++)
			nextStageBits.push_back(allSols[i + 1][j].col_bits);
		for (int j = 0; j < allSols[i].size(); j++)
		{
			column_sol_t sol = allSols[i][j];
			// propagate remaining bits from current stage to the next stage
			for (int k = 0; k < sol.col_iRB_sel.size(); k++)
				nextStageBits[j][sol.col_oRB_pos[k]] = stageBits[j][sol.col_iRB_sel[k]];
			// propagate AC42 bits from current stage to the next stage
			for (int k = 0; k < sol.col_iAC42_sel.size(); k++)
			{
				vector<int> inputs(4, 0), outputs;
				inputs[0] = stageBits[j][sol.col_iAC42_sel[k][0]];
				inputs[1] = stageBits[j][sol.col_iAC42_sel[k][1]];
				inputs[2] = stageBits[j][sol.col_iAC42_sel[k][2]];
				inputs[3] = stageBits[j][sol.col_iAC42_sel[k][3]];
				outputs = AC42_Op(inputs);
				nextStageBits[j][sol.col_oAC42_pos[k][0]] = outputs[0];
				nextStageBits[j][sol.col_oAC42_pos[k][1]] = outputs[1];
				//cout << "AC42 out1: " << i + 1 << "," << j << "," << sol.col_oAC42_pos[k][0] << " out2: " << i + 1 << "," << j << "," << sol.col_oAC42_pos[k][1] << endl;
			}
			// propagate AC32 bits from current stage to the next stage
			for (int k = 0; k < sol.col_iAC32_sel.size(); k++)
			{
				vector<int> inputs(3, 0), outputs;
				inputs[0] = stageBits[j][sol.col_iAC32_sel[k][0]];
				inputs[1] = stageBits[j][sol.col_iAC32_sel[k][1]];
				inputs[2] = stageBits[j][sol.col_iAC32_sel[k][2]];
				outputs = AC32_Op(inputs);
				nextStageBits[j][sol.col_oAC32_pos[k][0]] = outputs[0];
				nextStageBits[j][sol.col_oAC32_pos[k][1]] = outputs[1];
				//cout << "AC32 out1: " << i + 1 << "," << j << "," << sol.col_oAC32_pos[k][0] << " out2: " << i + 1 << "," << j << "," << sol.col_oAC32_pos[k][1] << endl;
			}
			// propagate FA bits from current stage to the next stage
			for (int k = 0; k < sol.col_iF_sel.size(); k++)
			{
				vector<int> inputs(3, 0), outputs;
				inputs[0] = stageBits[j][sol.col_iF_sel[k][0]];
				inputs[1] = stageBits[j][sol.col_iF_sel[k][1]];
				inputs[2] = stageBits[j][sol.col_iF_sel[k][2]];
				outputs = F_Op(inputs);
				nextStageBits[j][sol.col_oF_pos[k][0]] = outputs[0];
				nextStageBits[j - 1][sol.col_oF_pos[k][1]] = outputs[1];
				//cout << "F out1: " << i + 1 << "," << j << "," << sol.col_oF_pos[k][0] << " out2: " << i + 1 << "," << j - 1 << "," << sol.col_oF_pos[k][1] << endl;
			}
			// propagate HA bits from current stage to the next stage
			for (int k = 0; k < sol.col_iH_sel.size(); k++)
			{
				vector<int> inputs(2, 0), outputs;
				inputs[0] = stageBits[j][sol.col_iH_sel[k][0]];
				inputs[1] = stageBits[j][sol.col_iH_sel[k][1]];
				outputs = H_Op(inputs);
				nextStageBits[j][sol.col_oH_pos[k][0]] = outputs[0];
				nextStageBits[j - 1][sol.col_oH_pos[k][1]] = outputs[1];
				//cout << "H out1: " << i + 1 << "," << j << "," << sol.col_oH_pos[k][0] << " out2: " << i + 1 << "," << j - 1 << "," << sol.col_oH_pos[k][1] << endl;
			}

		}
		stageBits = nextStageBits;
	}
	return stageBits;
	
}

int final_addition(vector<vector<int>> final_outputs)
{
	int sum = 0;
	for (int i = 0; i < final_outputs.size(); i++)
	{
		int order = final_outputs.size() - i - 1;
		int sum_col = 0;
		if (final_outputs[i].size() == 1)
			sum_col += final_outputs[i][0];
		else if (final_outputs[i].size() == 2)
			sum_col += final_outputs[i][0] + final_outputs[i][1];
		sum += pow(2, order) * (sum_col);
	}
	return sum;
}

int final_mult(int mult1, int mult2, vector<vector<column_sol_t>> allSols)
{
	vector<int> mult_bin1, mult_bin2;
	for (int i = 0; i < MULT_SIZE; i++)
	{
		mult_bin1.push_back(mult1 % 2);
		mult1 /= 2;
		mult_bin2.push_back(mult2 % 2);
		mult2 /= 2;
	}
	//initialize bit matrix
	vector<vector<int>> bv_init(2 * MULT_SIZE + stages_num - 1, vector<int>{});
	for (int i = 2 * MULT_SIZE + stages_num - 2; i >= 0; i--)
	{
		int index = 2 * MULT_SIZE + stages_num - 2 - i;
		for (int j = 0; j <= min(index, MULT_SIZE - 1); j++)
		{
			int k = index - j;
			if (k > MULT_SIZE - 1)
				continue;
			bv_init[i].push_back(mult_bin1[j] * mult_bin2[k]);
			//cout << j << "*" << k << " ";
		}
		//cout << endl;
	}
	return final_addition(OptCompressTree(bv_init, allSols));
}

