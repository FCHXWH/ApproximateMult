#pragma once
#ifndef __CONNECTIONORDEROPT_H__
#define __CONNECTIONORDEROPT_H__
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "gurobi_c++.h"
using namespace std;
typedef struct sols_t
{
	vector<vector<int>> Vs;
	vector<vector<int>> Fs, Hs;
	vector<vector<int>> AC42, AC32;
	sols_t()
	{

	}

	sols_t(vector<vector<int>> vs, vector<vector<int>> fs, vector<vector<int>> hs, vector<vector<int>> ac42, vector<vector<int>> ac32)
	{
		Vs = vs;
		Fs = fs;
		Hs = hs;
		AC42 = ac42;
		AC32 = ac32;
	}
}Sols;

typedef struct vars_t
{
	vector<vector<GRBVar>> RemainBits;
	vector<vector<GRBVar>> AC42;
	vector<vector<GRBVar>> AC32;
	vector<vector<GRBVar>> f;
	vector<vector<GRBVar>> h;
	//vector<vector<GRBVar>> fc;
	//vector<vector<GRBVar>> hc;
	vector<GRBVar> Probs;
	int nrb_ls, nac42_ls, nac32_ls, nfs_ls, nhs_ls, nfc_ls, nhc_ls;
	int nrb, nac42, nac32, nF, nH;
	vars_t()
	{

	}

	vars_t(int nBits, int nAC42, int nAC32, int nf, int nh)
	{
		int nRemainbits = nBits - 4 * nAC42 - 3 * nAC32 - 3 * nf - 2 * nh;
		nrb = nRemainbits; nac42 = nAC42; nac32 = nAC32; nF = nf; nH = nh;
		// reserve the memory space of probabilities
		Probs.resize(nBits);
		// reserve the memory space of remaining bits
		RemainBits.resize(nBits);
		for (auto& RB : RemainBits)
			RB.resize(nRemainbits);
		// reserve the memory space of AC42
		AC42.resize(nBits);
		for (auto& ac42 : AC42)
			ac42.resize(4 * nAC42);
		// reserve the memory space of AC32
		AC32.resize(nBits);
		for (auto& ac32 : AC32)
			ac32.resize(3 * nAC32);
		// reserve the memory space of f
		f.resize(nBits);
		for (auto& fa : f)
			fa.resize(3 * nf);
		// reserve the memory space of h
		h.resize(nBits);
		for (auto& ha : h)
			ha.resize(2 * nh);
	}
}Vars_Column;
extern vector<vector<double>> allcoefs;
void generate_variables_of_ConnectOrderOpt(vector<vector<Vars_Column>>& allVars, Sols sol, GRBModel& model);
void Initial_probs_of_V(vector<Vars_Column> Vars_stage0, GRBModel& model);
vector<vector<GRBVar>> generate_constraints_of_ConnectOrderOpt(vector<vector<Vars_Column>> allVars, int isApprox, GRBModel& model);
void Init_EP_Coefs(vector<double>& coefs, string solname);
#endif // !__CONNECTIONORDEROPT_H__

