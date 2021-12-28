#pragma once
#ifndef __GA2HDL_H__
#define __GA2HDL_H__
#include <json/json.h>
#include <iostream>
#include <fstream>
#include "GAACT.h"
#include "Multiplier_Tree_Optimization.h"
#include <vector>
using namespace std;
using namespace Json;
typedef struct column_sol_t_ column_sol_t;
struct column_sol_t_
{
	vector<int> col_bits;
	// this stage: yi_j
	int nrb_ls, nac42_ls, nac32_ls, nfs_ls, nhs_ls, nfc_ls, nhc_ls;
	// this stage: zi_j
	vector<vector<int>> col_iAC42_sel;
	vector<vector<int>> col_iAC32_sel;
	vector<vector<int>> col_iF_sel;
	vector<vector<int>> col_iH_sel;
	vector<int> col_iRB_sel;
	// next stage: yi+1_j
	/*vector<vector<int>> col_oAC42_pos;
	vector<vector<int>> col_oAC32_pos;
	vector<vector<int>> col_oF_pos;
	vector<vector<int>> col_oH_pos;
	vector<int> col_oRB_pos;*/
	column_sol_t_()
	{
		nrb_ls = 0; nac42_ls = 0; nac32_ls = 0; nfs_ls = 0; nhs_ls = 0; nfc_ls = 0; nhc_ls = 0;
	}

	//struct column_sol_t_(Vars_Column col)
	//{
	//	int size = col.Probs.size();
	//	
	//	col_bits.resize(size, 0);
	//	col_iAC42_sel.resize(col.nac42, vector<int>(4, 0));
	//	col_iAC32_sel.resize(col.nac32, vector<int>(3, 0));
	//	col_iF_sel.resize(col.nF, vector<int>(3, 0));
	//	col_iH_sel.resize(col.nH, vector<int>(2, 0));
	//	col_iRB_sel.resize(col.nrb, 0);
	//	//// for next stage yi+1_j
	//	//col_oAC42_pos.resize(col.nac42, vector<int>(2, 0));
	//	//col_oAC32_pos.resize(col.nac32, vector<int>(2, 0));
	//	//col_oF_pos.resize(col.nF, vector<int>(2, 0));
	//	//col_oH_pos.resize(col.nH, vector<int>(2, 0));
	//	//col_oRB_pos.resize(col.nrb, 0);
	//	
	//	nrb_ls = col.nrb_ls;
	//	nac42_ls = col.nac42_ls;
	//	nac32_ls = col.nac32_ls;
	//	nfs_ls = col.nfs_ls;
	//	nhs_ls = col.nhs_ls;
	//	nfc_ls = col.nfc_ls;
	//	nhc_ls = col.nhc_ls;
	//};
	int TotalBits()
	{
		return nrb_ls + nac42_ls + nac32_ls + nfs_ls + nhs_ls + nfc_ls + nhc_ls;
	}
};

vector<vector<column_sol_t>> GA2ColSol(Gene g, ACT_Par_t Pars);
void WriteToJson(vector<vector<column_sol_t>> allSols, string file);
#endif // !__GA2HDL_H__

