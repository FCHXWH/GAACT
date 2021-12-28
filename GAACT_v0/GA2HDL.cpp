#include "GA2HDL.h"
#include <assert.h>
#include <json/json.h>
#include <algorithm>

using namespace std;
// we need to obtain the data structure column_sol_t from GA solution
vector<vector<column_sol_t>> GA2ColSol(Gene g, ACT_Par_t Pars)
{
	// each column_sol_t contains the information of yi_j () and 
	vector<vector<column_sol_t>> allSols;
	int nStages = Pars.nCompressors.size(), nColumns = Pars.nCompressors[0].size();
	allSols.resize(nStages + 1, vector<column_sol_t>(nColumns, column_sol_t()));
	// add the first stage for gene g
	for (int j = 0; j < nColumns; j++)
	{
		int size;
		if (j < Pars.BitWidth)
			size = j + 1;
		else
			size = ((2 * Pars.BitWidth - j - 1) > 0) ? (2 * Pars.BitWidth - j - 1) : 0;
		Perm permu;
		permu.perm.resize(size);
		for (int k = 0; k < size; k++)
			permu.perm[k] = k;
		g.addGene(j, permu);
	}
	// initialize allSols
	for (int i = 0; i < nStages; i++)
	{
		for (int j = 0; j < nColumns; j++)
		{
			int iPerm = i * nColumns + j;
			int startIndex = 0, endIndex = 0;
			if (i)
				assert(allSols[i][j].TotalBits() == g.g[iPerm].perm.size());
			allSols[i][j].col_bits.resize(g.g[iPerm].perm.size());
			// exact 2:2 compressors
			allSols[i + 1][j].nhs_ls = Pars.nCompressors[i][j][0]; allSols[i + 1][j + 1].nhc_ls = Pars.nCompressors[i][j][0];
			allSols[i][j].col_iH_sel.resize(Pars.nCompressors[i][j][0]);
			for (int k = 0; k < Pars.nCompressors[i][j][0]; k++)
			{
				endIndex = startIndex + 2;
				allSols[i][j].col_iH_sel[k].resize(2);
				allSols[i][j].col_iH_sel[k][0] = g.g[iPerm].perm[startIndex]; allSols[i][j].col_iH_sel[k][1] = g.g[iPerm].perm[startIndex + 1];
				startIndex = endIndex;
			}
			// exact 3:2 compressors
			allSols[i + 1][j].nfs_ls = Pars.nCompressors[i][j][1]; allSols[i + 1][j + 1].nfc_ls = Pars.nCompressors[i][j][1];
			allSols[i][j].col_iF_sel.resize(Pars.nCompressors[i][j][1]);
			for (int k = 0; k < Pars.nCompressors[i][j][1]; k++)
			{
				endIndex = startIndex + 3;
				allSols[i][j].col_iF_sel[k].resize(3);
				allSols[i][j].col_iF_sel[k][0] = g.g[iPerm].perm[startIndex]; allSols[i][j].col_iF_sel[k][1] = g.g[iPerm].perm[startIndex + 1];
				allSols[i][j].col_iF_sel[k][2] = g.g[iPerm].perm[startIndex + 2];
				startIndex = endIndex;
			}
			// approximate 3:2 compressors
			allSols[i + 1][j].nac32_ls = 2 * Pars.nCompressors[i][j][2];
			allSols[i][j].col_iAC32_sel.resize(Pars.nCompressors[i][j][2]);
			for (int k = 0; k < Pars.nCompressors[i][j][2]; k++)
			{
				endIndex = startIndex + 3;
				allSols[i][j].col_iAC32_sel[k].resize(3);
				allSols[i][j].col_iAC32_sel[k][0] = g.g[iPerm].perm[startIndex]; allSols[i][j].col_iAC32_sel[k][1] = g.g[iPerm].perm[startIndex + 1];
				allSols[i][j].col_iAC32_sel[k][2] = g.g[iPerm].perm[startIndex + 2];
				startIndex = endIndex;
			}
			// approximate 4:2 compressors
			allSols[i + 1][j].nac42_ls = 2 * Pars.nCompressors[i][j][3];
			allSols[i][j].col_iAC42_sel.resize(Pars.nCompressors[i][j][3]);
			for (int k = 0; k < Pars.nCompressors[i][j][3]; k++)
			{
				endIndex = startIndex + 4;
				allSols[i][j].col_iAC42_sel[k].resize(4);
				allSols[i][j].col_iAC42_sel[k][0] = g.g[iPerm].perm[startIndex]; allSols[i][j].col_iAC42_sel[k][1] = g.g[iPerm].perm[startIndex + 1];
				allSols[i][j].col_iAC42_sel[k][2] = g.g[iPerm].perm[startIndex + 2]; allSols[i][j].col_iAC42_sel[k][3] = g.g[iPerm].perm[startIndex + 3];
				startIndex = endIndex;
			}
			// dummy compressors
			allSols[i + 1][j].nrb_ls = g.g[iPerm].perm.size() - startIndex;
			for (int k = startIndex; k < g.g[iPerm].perm.size(); k++)
				allSols[i][j].col_iRB_sel.push_back(g.g[iPerm].perm[k]);
		}
	}
	for (int j = 0; j < nColumns; j++)
	{
		int total_bits = allSols[nStages][j].nac32_ls + allSols[nStages][j].nac42_ls + allSols[nStages][j].nfc_ls
			+ allSols[nStages][j].nfs_ls + allSols[nStages][j].nhc_ls + allSols[nStages][j].nhs_ls
			+ allSols[nStages][j].nrb_ls;
		allSols[nStages][j].col_bits.resize(total_bits);
		allSols[nStages][j].col_iRB_sel.resize(total_bits);
		for (int k = 0; k < total_bits; k++)
			allSols[nStages][j].col_iRB_sel[k] = k;
	}
	return allSols;
}

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