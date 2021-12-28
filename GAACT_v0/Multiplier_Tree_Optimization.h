#pragma once
#ifndef __MULTIPLIER_TREE_OPTIMIZATION_H__
#define __MULTIPLIER_TREE_OPTIMIZATION_H__
#include <iostream>
#include <vector>
#include <sstream>
#include "gurobi_c++.h"
using namespace std;
typedef struct MUOPT_Par_t_ MUOPT_Par_t;
struct MUOPT_Par_t_
{
	int MULT_SIZE, nStages;
	int Opt_Mode; //0->Constrain error & Optimize area
				  //1->Constrain area & Optimize error
				  //2->CoOptimize area & error
	int MED_bound, Area_bound, Weight;
	string filename;
	vector<int> input_patterns;
	int simulator_approach;
	int Time_Bound_s;
	vector<vector<double>> allcoefs;
	int EP_Approx;
	double MED_ac42, MED_ac32;
	string json_path, json_filename, log_path, log_filename;
	string cmp_path;
};

void init_Pars(int argc, char* argv[], MUOPT_Par_t& Par);
void generate_variables_of_multiplier(MUOPT_Par_t Par, vector<vector<GRBVar>>& variables_fa, vector<vector<GRBVar>>& variables_ha,
	vector<vector<GRBVar>>& variables_V, vector<vector<GRBVar>>& variables_ac32,
	vector<vector<GRBVar>>& variables_ac42, vector<vector<GRBVar>>& variables_error, GRBModel& model);
void Initial_constraints_of_V(MUOPT_Par_t Pars, vector<vector<GRBVar>> variables_V, GRBModel& model);
void generate_constraints_of_h_f_V(MUOPT_Par_t Pars, vector<vector<GRBVar>> variables_fa, vector<vector<GRBVar>> variables_ha, vector<vector<GRBVar>> variables_V,
	vector<vector<GRBVar>> variables_ac32, vector<vector<GRBVar>> variables_ac42, vector<vector<GRBVar>> variables_error, GRBModel& model);
void generate_constraints_of_last_stage(MUOPT_Par_t Pars, vector<vector<GRBVar>> variables_V, GRBModel& model);
void generate_cost_of_adders(vector<vector<GRBVar>> variables_fa, vector<vector<GRBVar>> variables_ha, GRBVar& fa_num, GRBVar& ha_num,
	vector<vector<GRBVar>> variables_ac32, vector<vector<GRBVar>> variables_ac42, GRBVar& ac32_num, GRBVar& ac42_num, GRBModel& model);
void generate_total_error(vector<vector<GRBVar>> variables_error, GRBVar& E, GRBModel& model);

#endif // !__MULTIPLIER_TREE_OPTIMIZATION_H__

