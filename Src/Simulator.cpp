#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include "Multiplier_Tree_Optimization.h"
#include <fstream>
using namespace std;

void sci2db(string& strSci, string& value_int)
{
	int value_s;
	stringstream ss;
	int value;
	int    iPower = 0;  //ÃÝ
	double dMntsa = 0;  //Î²Êý
	double dCoefficient = 1;  //ÏµÊý

	std::string strPower, strMntsa;

	if (std::string::npos == strSci.find("e"))
	{
		ss << round(atof(strSci.c_str()));
		ss >> value_int;
	}
	else
	{
		strMntsa = strSci.substr(0, strSci.find("e"));
		strPower = strSci.substr(strSci.find("e") + 1);

		dMntsa = atof(strMntsa.c_str());
		iPower = atoi(strPower.c_str());

		while (iPower != 0)
		{
			if (iPower > 0)
			{
				dCoefficient *= 10;
				iPower--;
			}
			else
			{
				dCoefficient *= 0.1;
				iPower++;
			}
		}
		value_s = round(dMntsa * dCoefficient);
		ss << value_s;
		ss >> value_int;
		//return dMntsa * dCoefficient;
	}

}

//input: p0 p1
//output: sum cout
void sim_h()
{
	int nInputs = 2;
	double sum_p = 0.0, co_p = 0.0;
	cout << "p0 p1 sum cout" << endl;
	for (int i = 0; i < pow(2, nInputs); i++)
	{
		int input0 = i % 2;
		int input1 = (i / 2) % 2;
		int sum, co;
		sum = input0 ^ input1;
		co = input0 * input1;
		cout << input0 << "  " << input1 << "   " << sum << "   " << co << endl;
		if (sum == 1)
			sum_p += (0.25 * input0 + 0.75 * (1 - input0)) * (0.25 * input1 + 0.75 * (1 - input1));
		if (co == 1)
			co_p += (0.25 * input0 + 0.75 * (1 - input0)) * (0.25 * input1 + 0.75 * (1 - input1));
	}
	cout << "The probability of sum to be 1 is " << sum_p << endl;
	cout << "The probability of co to be 1 is " << co_p << endl;
}

//input: p0 p1
//output: sum cout
void sim_f()
{
	int nInputs = 3;
	double sum_p = 0.0, co_p = 0.0;
	cout << "p0 p1 p2 sum cout" << endl;
	for (int i = 0; i < pow(2, nInputs); i++)
	{
		int input0 = i % 2;
		int input1 = (i / 2) % 2;
		int input2 = (i / 2 / 2) % 2;
		int sum, co;
		sum = input0 ^ input1 ^ input2;
		co = input0 * input1 || input0 * input2 || input1 * input2;
		cout << input0 << "  " << input1 << "  " << input2 << "   " << sum << "   " << co << endl;
		if (sum == 1)
			sum_p += (0.25 * input0 + 0.75 * (1 - input0)) * (0.25 * input1 + 0.75 * (1 - input1)) * (0.25 * input2 + 0.75 * (1 - input2));
		if (co == 1)
			co_p += (0.25 * input0 + 0.75 * (1 - input0)) * (0.25 * input1 + 0.75 * (1 - input1)) * (0.25 * input2 + 0.75 * (1 - input2));
	}
	cout << "The probability of sum to be 1 is " << sum_p << endl;
	cout << "The probability of co to be 1 is " << co_p << endl;
}

vector<double> max_error_AC42_propagation(vector<double> probs)
{
	int nInputs = 4;
	int inputs[4];
	double MED = 0.0;
	double ER1 = 0.0, ER2 = 0.0;
	//outputs: output0,output1,output2,output3
	int output0, output1, output2, output3;
	cout << "p0 p1 p2 p3 w1 w2" << endl;
	for (int i = 0; i < pow(2, nInputs); i++)
	{
		inputs[0] = i % 2;
		inputs[1] = i / 2 % 2;
		inputs[2] = i / 2 / 2 % 2;
		inputs[3] = i / 2 / 2 / 2 % 2;
		output0 = (inputs[0] * inputs[1]) * (inputs[2] || inputs[3]);
		output1 = (inputs[2] * inputs[3]) * (inputs[0] || inputs[1]);
		output2 = (inputs[0] * inputs[1]) || (inputs[2] || inputs[3]);
		output3 = (inputs[2] * inputs[3]) || (inputs[0] || inputs[1]);
		cout << inputs[0] << "  " << inputs[1] << "  " << inputs[2] << "   " << inputs[3] << "   " << output2 << "   " << output3 << endl;
		double prob_term = 1.0;
		for (int j = 0; j < probs.size(); j++)
		{
			prob_term *= probs[j] * inputs[j] + (1 - probs[j]) * (1 - inputs[j]);
		}
		MED += (output0 + output1) * prob_term;
		ER1 += (output2)*prob_term;
		ER2 += (output3)*prob_term;
	}
	cout << "The mean error distance is " << MED << endl;
	cout << "The probability of w1 to be 1 is " << ER1 << endl;
	cout << "The probability of w2 to be 1 is " << ER2 << endl;
	return { MED,ER1,ER2 };
}

vector<double> max_error_AC32_propagation(vector<double> probs)
{
	int nInputs = 3;
	int inputs[3];
	double MED = 0.0;
	double ER1 = 0.0, ER2 = 0.0;
	//outputs: output0,output1,output2
	int output0, output1, output2;
	cout << "p0 p1 p2 w1 w2" << endl;
	for (int i = 0; i < pow(2, nInputs); i++)
	{
		inputs[0] = i % 2;
		inputs[1] = i / 2 % 2;
		inputs[2] = i / 2 / 2 % 2;
		output0 = inputs[0] * inputs[1] * inputs[2];
		output1 = (inputs[0] * inputs[1]) || inputs[2];
		output2 = inputs[0] || inputs[1];
		cout << inputs[0] << "  " << inputs[1] << "  " << inputs[2] << "   " << output1 << "   " << output2 << endl;
		double prob_term = 1.0;
		for (int j = 0; j < probs.size(); j++)
		{
			prob_term *= probs[j] * inputs[j] + (1 - probs[j]) * (1 - inputs[j]);
		}
		MED += (output0) * prob_term;
		ER1 += (output1)*prob_term;
		ER2 += (output2)*prob_term;
	}
	cout << "The mean error distance is " << MED << endl;
	cout << "The probability of w1 to be 1 is " << ER1 << endl;
	cout << "The probability of w2 to be 1 is " << ER2 << endl;
	return { MED,ER1,ER2 };
}

void get_name_value(string tmp, int& f, int& h, int& ac32, int& ac42, double& error, int& V, ifstream& f0)
{
	stringstream ss;
	//fi_j
	string tmp_s, tmp_v;
	tmp_s =	tmp.substr(tmp.find(' ') + 1, tmp.size() - tmp.find(' '));
	sci2db(tmp_s, tmp_v);
	ss << tmp_v;
	ss >> f;
	ss.clear(); ss.str("");
	//hi_j
	getline(f0, tmp);
	tmp_s = tmp.substr(tmp.find(' ') + 1, tmp.size() - tmp.find(' '));
	sci2db(tmp_s, tmp_v);
	ss << tmp_v;
	ss >> h;
	ss.clear(); ss.str("");
	//32aci_j
	getline(f0, tmp);
	tmp_s = tmp.substr(tmp.find(' ') + 1, tmp.size() - tmp.find(' '));
	sci2db(tmp_s, tmp_v);
	ss << tmp_v;
	ss >> ac32;
	ss.clear(); ss.str("");
	//42aci_j
	getline(f0, tmp);
	tmp_s = tmp.substr(tmp.find(' ') + 1, tmp.size() - tmp.find(' '));
	sci2db(tmp_s, tmp_v);
	ss << tmp_v;
	ss >> ac42;
	ss.clear(); ss.str("");
	//ei_j
	getline(f0, tmp);
	tmp_s = tmp.substr(tmp.find(' ') + 1, tmp.size() - tmp.find(' '));
	sci2db(tmp_s, tmp_v);
	ss << tmp_v;
	ss >> error;
	ss.clear(); ss.str("");
	//Vi_j
	getline(f0, tmp);
	tmp_s = tmp.substr(tmp.find(' ') + 1, tmp.size() - tmp.find(' '));
	sci2db(tmp_s, tmp_v);
	ss << tmp_v;
	ss >> V;
	ss.clear(); ss.str("");
}

void get_index(string tmp, int& index)
{
	stringstream ss;
	string s_index = tmp.substr(tmp.find('_') + 1, tmp.find(' ') - tmp.find('_'));
	ss << s_index;
	ss >> index;
	ss.clear();
}


void ConstructCT(vector<vector<int>>& F, vector<vector<int>>& H, vector<vector<int>>& AC32, 
	vector<vector<int>>& AC42, vector<vector<double>>& E, vector<vector<int>>& V)
{
	ifstream f0;
	string tmp;

	f0.open((filename + ".sol").c_str());
	int i = 0, j = 0;
	while (getline(f0, tmp))
	{
		if (tmp[0] == 'f' && i < stages_num)
		{
			int f, h, ac32, ac42, v;
			double error;
			int index;
			get_name_value(tmp, f, h, ac32, ac42, error, v, f0);
			F[i][j] = f; H[i][j] = h; AC32[i][j] = ac32; AC42[i][j] = ac42; E[i][j] = error;
			V[i][j] = v;
			get_index(tmp, index);
			if (index == 2 * MULT_SIZE + stages_num - 2)
			{
				i++; j = 0;
			}
			else
				j++;
		}
	}
	f0.close();


}

vector<pair<int, int>> CompressorTree(vector<int> PI1, vector<int> PI2)
{
	vector<vector<int>> F(stages_num, vector<int>(2 * MULT_SIZE + stages_num - 1, 0));
	vector<vector<int>> H(stages_num, vector<int>(2 * MULT_SIZE + stages_num - 1, 0));
	vector<vector<int>> AC32(stages_num, vector<int>(2 * MULT_SIZE + stages_num - 1, 0));
	vector<vector<int>> AC42(stages_num, vector<int>(2 * MULT_SIZE + stages_num - 1, 0));
	vector<vector<double>> E(stages_num, vector<double>(2 * MULT_SIZE + stages_num - 1, 0.0));
	vector<vector<int>> V(stages_num, vector<int>(2 * MULT_SIZE + stages_num - 1, 0));
	ConstructCT(F, H, AC32, AC42, E, V);
	/*cout << "F H AC32 AC42 V" << endl;
	for (int i = 0; i < stages_num; i++)
	{
		for (int j = 0; j < V[0].size(); j++)
		{
			cout << F[i][j] << " " <<  H[i][j] << " " << AC32[i][j] << " " << AC42[i][j] << " " << V[i][j] << "  ";
		}
		cout << endl;
	}*/
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
			bv_init[i].push_back(PI1[j] * PI2[k]);
		}
	}

	vector<vector<int>> bv_stage, bv_nextStage = bv_init;
	for (int i = 0; i < stages_num; i++)
	{
		bv_stage = bv_nextStage;
		for (int t = 0; t < bv_nextStage.size(); t++)
			bv_nextStage[t].clear();
		// operate order: ac42 ac32 f h
		// dot order: remaining carrys ac42 ac32 f h 
		for (int j = bv_stage.size() - 1; j >= 0; j--)
		{
			int ac42 = AC42[i][j]; int ac32 = AC32[i][j]; int f = F[i][j]; int h = H[i][j];
			int remaining = bv_stage[j].size() - 4 * ac42 - 3 * ac32 - 3 * f - 2 * h;
			// operation for remaining dots
			//bv_nextStage[j].insert(bv_nextStage[j].begin(), bv_stage[j].end() - remaining, bv_stage[j].end());
			
			// operation of compressor ac42
			int startInd_a4 = 0, endInd_a4 = startInd_a4 + 4 * ac42 - 1;
			for (int k = 0; k < ac42; k++)
			{
				bv_nextStage[j].push_back(bv_stage[j][4 * k] * bv_stage[j][4 * k + 1] || bv_stage[j][4 * k + 2] || bv_stage[j][4 * k + 3]);
				bv_nextStage[j].push_back(bv_stage[j][4 * k + 2] * bv_stage[j][4 * k + 3] || bv_stage[j][4 * k] || bv_stage[j][4 * k + 1]);
			}
			int startInd_a3 = endInd_a4 + 1, endInd_a3 = startInd_a3 + 3 * ac32 - 1;
			// operation of compressor f
			int startInd_f = endInd_a3 + 1, endInd_f = startInd_f + 3 * f - 1;
			for (int k = 0; k < f; k++)
			{
				bv_nextStage[j - 1].push_back(bv_stage[j][3 * k + startInd_f] * bv_stage[j][3 * k + 1 + startInd_f]
					|| bv_stage[j][3 * k + startInd_f] * bv_stage[j][3 * k + 2 + startInd_f]
					|| bv_stage[j][3 * k + 1 + startInd_f] * bv_stage[j][3 * k + 2 + startInd_f]);
				bv_nextStage[j].push_back(bv_stage[j][3 * k + startInd_f] ^ bv_stage[j][3 * k + 1 + startInd_f] ^ bv_stage[j][3 * k + 2 + startInd_f]);
			}

			// operation of compressor ac32
			
			for (int k = 0; k < ac32; k++)
			{
				bv_nextStage[j].push_back(bv_stage[j][3 * k + startInd_a3] || bv_stage[j][3 * k + 1 + startInd_a3]);
				//bv_nextStage[j].push_back(bv_stage[j][3 * k + startInd_a3] * bv_stage[j][3 * k + 1 + startInd_a3] || bv_stage[j][3 * k + 2 + startInd_a3]);
			}

			// operation of compressor h
			int startInd_h = endInd_f + 1, endInd_h = startInd_h + 2 * h - 1;
			for (int k = 0; k < h; k++)
			{
				bv_nextStage[j - 1].push_back(bv_stage[j][2 * k + startInd_h] * bv_stage[j][2 * k + 1 + startInd_h]);
				bv_nextStage[j].push_back(bv_stage[j][2 * k + startInd_h] ^ bv_stage[j][2 * k + 1 + startInd_h]);
			}

			for (int k = 0; k < ac32; k++)
			{
				//bv_nextStage[j].push_back(bv_stage[j][3 * k + startInd_a3] || bv_stage[j][3 * k + 1 + startInd_a3]);
				bv_nextStage[j].push_back(bv_stage[j][3 * k + startInd_a3] * bv_stage[j][3 * k + 1 + startInd_a3] || bv_stage[j][3 * k + 2 + startInd_a3]);
			}
			// replace the first carry signals to the end 
			/*if (j != bv_stage.size() - 1)
			{
				int nCarry = F[i][j + 1] + H[i][j + 1];
				for (int ca = 0; ca < nCarry; ca++)
					bv_nextStage[j].push_back(bv_nextStage[j][ca]);
				for (int ca = 0; ca < nCarry; ca++)
					bv_nextStage[j].erase(bv_nextStage[j].begin());
			}*/
			

			// operation of compressor ac42
			/*int startInd_a4 = 0, endInd_a4 = 4 * ac42 - 1;
			for (int k = 0; k < ac42; k++)
			{
				bv_nextStage[j].push_back(bv_stage[j][4 * k] * bv_stage[j][4 * k + 1] || bv_stage[j][4 * k + 2] || bv_stage[j][4 * k + 3]);
				bv_nextStage[j].push_back(bv_stage[j][4 * k + 2] * bv_stage[j][4 * k + 3] || bv_stage[j][4 * k] || bv_stage[j][4 * k + 1]);
			}*/
			//// operation of compressor ac32
			//int startInd_a3 = endInd_a4 + 1, endInd_a3 = startInd_a3 + 3 * ac32 - 1;
			//for (int k = 0; k < ac32; k++)
			//{
			//	bv_nextStage[j].push_back(bv_stage[j][3 * k + startInd_a3] * bv_stage[j][3 * k + 1 + startInd_a3] || bv_stage[j][3 * k + 2 + startInd_a3]);
			//	bv_nextStage[j].push_back(bv_stage[j][3 * k + startInd_a3] || bv_stage[j][3 * k + 1 + startInd_a3]);
			//}
			
			

			// operation for remaining dots
			bv_nextStage[j].insert(bv_nextStage[j].end(), bv_stage[j].end() - remaining, bv_stage[j].end());
		}
	}
	vector<pair<int, int>> ret;
	for (auto vec : bv_nextStage)
	{
		if (vec.size() == 1)
			ret.push_back(make_pair(vec[0], 0));
		else if (vec.size() > 1)
			ret.push_back(make_pair(vec[0], vec[1]));
	}
	return ret;
}

int final_addition(vector<pair<int, int>> final_outputs)
{
	int sum = 0;
	for (int i = 0; i < final_outputs.size(); i++)
	{
		int order = final_outputs.size() - i - 1;
		sum += pow(2, order) * (final_outputs[i].first + final_outputs[i].second);
	}
	return sum;
}

int final_mult(int mult1, int mult2)
{
	vector<int> mult_bin1, mult_bin2;
	for (int i = 0; i < MULT_SIZE; i++)
	{
		mult_bin1.push_back(mult1 % 2);
		mult1 /= 2;
		mult_bin2.push_back(mult2 % 2);
		mult2 /= 2;
	}
	return final_addition(CompressorTree(mult_bin1, mult_bin2));
}

void generate_MEDs(vector<double>& MED_ac42, vector<double>& MED_ac32)
{
	vector<double> probs_ac42 = { 0.25,0.25,0.25,0.25 };
	vector<double> probs_ac32 = { 0.25,0.25,0.25 };
	/*vector<double> MED_ac42(stages_num, 0.0);
	vector<double> MED_ac32(stages_num, 0.0);*/
	for (int st = 0; st < stages_num; st++)
	{
		vector<double> outProbs_ac32 = max_error_AC32_propagation(probs_ac32);
		vector<double> outProbs_ac42 = max_error_AC42_propagation(probs_ac42);
		cout << endl;
		cout << "The probability of " << st + 1 << "-th stage is " << probs_ac42[0] << endl;
		cout << "The MED of " << st + 1 << "-th stage's ac42 is " << outProbs_ac42[0] << endl;
		cout << "The MED of " << st + 1 << "-th stage's ac32 is " << outProbs_ac32[0] << endl;
		MED_ac42[st] = outProbs_ac42[0];
		MED_ac32[st] = outProbs_ac32[0];
		probs_ac42 = vector<double>(4, outProbs_ac42[1]);
		probs_ac32 = vector<double>(3, outProbs_ac42[1]);
	}
}