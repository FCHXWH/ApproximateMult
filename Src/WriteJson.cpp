#include "WriteJson.h"
#pragma warning(disable : 4996)
void WriteToJson(vector<vector<column_sol_t>> allSols, string file)
{
	StyledWriter writer;
	Value root;
	Value CompressorTree;
	for (int i = 0; i < allSols.size(); i++)
	{
		// each stage
		Value Stage;
		for (int j = 0; j < allSols[i].size(); j++)
		{
			Value Col, RBs, AC42, AC32, F, H;
			Col["col_len"] = allSols[i][j].col_bits.size();
			// Remaining bits input indexs
			vector<int> iRBs = allSols[i][j].col_iRB_sel;
			for (auto index : iRBs)
				RBs.append(index);
			Col["remainingBits"] = RBs;
			// AC42 input bits
			vector<vector<int>> iAC42s = allSols[i][j].col_iAC42_sel;
			for (auto iAC42 : iAC42s)
			{
				Value ac42;
				for (auto input : iAC42)
					ac42.append(input);
				AC42.append(ac42);
			}
			Col["AC42"] = AC42;
			// AC32 input bits
			vector<vector<int>> iAC32s = allSols[i][j].col_iAC32_sel;
			for (auto iAC32 : iAC32s)
			{
				Value ac32;
				for (auto input : iAC32)
					ac32.append(input);
				AC32.append(ac32);
			}
			Col["AC32"] = AC32;
			// F input bits
			vector<vector<int>> iFs = allSols[i][j].col_iF_sel;
			for (auto iF : iFs)
			{
				Value f;
				for (auto input : iF)
					f.append(input);
				F.append(f);
			}
			Col["F"] = F;
			//H input bits
			vector<vector<int>> iHs = allSols[i][j].col_iH_sel;
			for (auto iH : iHs)
			{
				Value h;
				for (auto input : iH)
					h.append(input);
				H.append(h);
			}
			Col["H"] = H;

			Col["outRemainingBitsNum"] = allSols[i][j].nrb_ls;
			Col["outAC42BitsNum"] = allSols[i][j].nac42_ls;
			Col["outAC32BitsNum"] = allSols[i][j].nac32_ls;
			Col["outFBitsNum"] = allSols[i][j].nfs_ls;
			Col["outHBitsNum"] = allSols[i][j].nhs_ls;
			Col["outFcBitsNum"] = allSols[i][j].nfc_ls;
			Col["outHcBitsNum"] = allSols[i][j].nhc_ls;
			Stage.append(Col);
		}
		string s = "stage" + to_string(i);
		CompressorTree[s] = Stage;
	}
	ofstream os;
	os.open(file.c_str());
	os << writer.write(CompressorTree);
	os.close();
}