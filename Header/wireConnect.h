#pragma once
#ifndef __WIRECONNECT_H__
#define __WIRECONNECT_H__
#include <iostream>
#include <fstream>
#include "ConnectionOrderOpt.h"
#include "Multiplier_Tree_Optimization.h"
#include <vector>
using namespace std;
typedef struct column_sol_t_ column_sol_t;
struct column_sol_t_
{
	vector<int> col_bits;
	vector<vector<int>> col_iAC42_sel;
	vector<vector<int>> col_iAC32_sel;
	vector<vector<int>> col_iF_sel;
	vector<vector<int>> col_iH_sel;
	vector<int> col_iRB_sel;
	
	vector<vector<int>> col_oAC42_pos;
	vector<vector<int>> col_oAC32_pos;
	vector<vector<int>> col_oF_pos;
	vector<vector<int>> col_oH_pos;
	vector<int> col_oRB_pos;
	int nrb_ls, nac42_ls, nac32_ls, nfs_ls, nhs_ls, nfc_ls, nhc_ls;
	struct column_sol_t_()
	{

	}

	struct column_sol_t_(Vars_Column col)
	{
		int size = col.Probs.size();
		// for this stage
		col_bits.resize(size, 0);
		col_iAC42_sel.resize(col.nac42, vector<int>(4, 0));
		col_iAC32_sel.resize(col.nac32, vector<int>(3, 0));
		col_iF_sel.resize(col.nF, vector<int>(3, 0));
		col_iH_sel.resize(col.nH, vector<int>(2, 0));
		col_iRB_sel.resize(col.nrb, 0);
		// for next stage
		col_oAC42_pos.resize(col.nac42, vector<int>(2, 0));
		col_oAC32_pos.resize(col.nac32, vector<int>(2, 0));
		col_oF_pos.resize(col.nF, vector<int>(2, 0));
		col_oH_pos.resize(col.nH, vector<int>(2, 0));
		col_oRB_pos.resize(col.nrb, 0);
		
		nrb_ls = col.nrb_ls;
		nac42_ls = col.nac42_ls;
		nac32_ls = col.nac32_ls;
		nfs_ls = col.nfs_ls;
		nhs_ls = col.nhs_ls;
		nfc_ls = col.nfc_ls;
		nhc_ls = col.nhc_ls;
	};
	void read_GRB_Solu(Vars_Column col, Vars_Column col_left);
};

vector<vector<column_sol_t>> Init_Sols(vector<vector<Vars_Column>> allVars);
vector<vector<int>> OptCompressTree(vector<vector<int>> stageBits, vector<vector<column_sol_t>> allSols);
int final_mult(int mult1, int mult2, vector<vector<column_sol_t>> allSols);
#endif // !__WIRECONNECT_H__

