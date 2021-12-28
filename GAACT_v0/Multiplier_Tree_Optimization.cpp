#include "Multiplier_Tree_Optimization.h"
#include <cmath>
using namespace std;

void init_Pars(int argc, char* argv[], MUOPT_Par_t& Par)
{
	Par.filename = "gurobi";
	Par.Time_Bound_s = 1200;
	Par.MULT_SIZE = 16;
	Par.nStages = 4;
	Par.simulator_approach = 1; // 0->MyHDL, 1->Simulator.out
	Par.Opt_Mode = 0;
	Par.Area_bound = 190;
	//Par.MED_bound = 115;
	//Par.MED_bound = 5400;
	Par.MED_bound = 200000;
	// Construct input pattern
	for (int i = 1; i <= 2 * Par.MULT_SIZE - 1; i++)
		Par.input_patterns.push_back(Par.MULT_SIZE - abs(i - Par.MULT_SIZE));
	Par.MED_ac42 = double(14) / 256;
	Par.MED_ac32 = double(1) / 64;
	// set the path of .json file
	Par.json_path = "/mnt/d/Microsoft\\ VS\\ Studio\\ 2019/Project/ApproxMULT_MyHDL_v4/ApproxMULT_MyHDL/";
	Par.json_filename = (!Par.Opt_Mode) ? "CompressorTree_CEOA" + to_string(int(Par.MED_bound)) + ".json"
		: ((Par.Opt_Mode == 1) ? ("CompressorTree_CAOE" + to_string(int(Par.Area_bound)) + ".json") :
			("CompressorTree_CoOpt" + to_string(Par.Weight) + ".json"));
	// set the path of .txt file
	if (!Par.Opt_Mode)
	{
		Par.log_path = "../../../CEOA/" + to_string(Par.MULT_SIZE) + "/" + to_string(int(Par.MED_bound)) + "/";
		Par.log_filename = "simulator_opt_error_" + to_string(Par.MED_bound) + ".txt";
	}
	else if (Par.Opt_Mode == 1)
	{
		Par.log_path = "../../../CAOE/" + to_string(Par.MULT_SIZE) + "/" + to_string(int(Par.Area_bound)) + "/";
		Par.log_filename = "simulator_opt_error_" + to_string(Par.Area_bound) + ".txt";
	}
	else
	{
		Par.log_path = "../../../CoOpt/" + to_string(Par.MULT_SIZE) + "/" + to_string(Par.Weight) + "/";
		Par.log_filename = "simulator_opt_error_" + to_string(Par.Weight) + ".txt";
	}
	Par.cmp_path = "../../../ExactMult/";
}

void generate_variables_of_multiplier(MUOPT_Par_t Par, vector<vector<GRBVar>>& variables_fa, vector<vector<GRBVar>>& variables_ha,
	vector<vector<GRBVar>>& variables_V, vector<vector<GRBVar>>& variables_ac32,
	vector<vector<GRBVar>>& variables_ac42, vector<vector<GRBVar>>& variables_error, GRBModel& model)
{
	for (int i = 0; i < Par.nStages + 1; i++)
	{
		vector<GRBVar> variables_fa_of_each_stage;
		vector<GRBVar> variables_ha_of_each_stage;
		vector<GRBVar> variables_ac32_of_each_stage;
		vector<GRBVar> variables_ac42_of_each_stage;
		vector<GRBVar> variables_error_of_each_stage;
		vector<GRBVar> variables_V_of_each_stage;
		stringstream ss;
		string tmp_s;
		for (int j = 0; j < Par.input_patterns.size() + Par.nStages; j++)
		{
			if (i != 0)
			{
				ss << "f" << i << "_" << j;
				ss >> tmp_s;
				ss.clear(); ss.str("");
				GRBVar variable_of_fa = model.addVar(0, INFINITY, 0, GRB_INTEGER, tmp_s);
				variables_fa_of_each_stage.push_back(variable_of_fa);

				ss << "h" << i << "_" << j;
				ss >> tmp_s;
				ss.clear(); ss.str("");
				GRBVar variable_of_ha = model.addVar(0, INFINITY, 0, GRB_INTEGER, tmp_s);
				variables_ha_of_each_stage.push_back(variable_of_ha);

				ss << "32ac" << i << "_" << j;
				ss >> tmp_s;
				ss.clear(); ss.str("");
				GRBVar variable_of_32ac = model.addVar(0, INFINITY, 0, GRB_INTEGER, tmp_s);
				variables_ac32_of_each_stage.push_back(variable_of_32ac);

				ss << "42ac" << i << "_" << j;
				ss >> tmp_s;
				ss.clear(); ss.str("");
				GRBVar variable_of_42ac = model.addVar(0, INFINITY, 0, GRB_INTEGER, tmp_s);
				variables_ac42_of_each_stage.push_back(variable_of_42ac);

				ss << "e" << i << "_" << j;
				ss >> tmp_s;
				ss.clear(); ss.str("");
				GRBVar variable_of_e = model.addVar(0, INFINITY, 0, GRB_CONTINUOUS, tmp_s);
				variables_error_of_each_stage.push_back(variable_of_e);
			}

			ss << "V" << i << "_" << j;
			ss >> tmp_s;
			ss.clear();
			GRBVar variable_of_V = model.addVar(0, INFINITY, 0, GRB_INTEGER, tmp_s);
			variables_V_of_each_stage.push_back(variable_of_V);
		}
		if (i != 0)
		{
			variables_fa.push_back(variables_fa_of_each_stage);
			variables_ha.push_back(variables_ha_of_each_stage);
			variables_ac32.push_back(variables_ac32_of_each_stage);
			variables_ac42.push_back(variables_ac42_of_each_stage);
			variables_error.push_back(variables_error_of_each_stage);
		}
		variables_V.push_back(variables_V_of_each_stage);

	}

}

void Initial_constraints_of_V(MUOPT_Par_t Pars, vector<vector<GRBVar>> variables_V, GRBModel& model)
{
	for (int i = 0; i < Pars.input_patterns.size() + Pars.nStages; i++)
	{
		if (i < Pars.input_patterns.size())
			model.addConstr(variables_V[0][i] == Pars.input_patterns[i]);
		else
			model.addConstr(variables_V[0][i] == 0);
	}
}

void generate_constraints_of_h_f_V(MUOPT_Par_t Pars, vector<vector<GRBVar>> variables_fa, vector<vector<GRBVar>> variables_ha, vector<vector<GRBVar>> variables_V,
	vector<vector<GRBVar>> variables_ac32, vector<vector<GRBVar>> variables_ac42, vector<vector<GRBVar>> variables_error, GRBModel& model)
{
	for (int i = 1; i < variables_V.size(); i++)
	{
		for (int j = 0; j < variables_V[i].size(); j++)
		{
			if (j == 0)
			{
				model.addConstr(variables_V[i][j] == 1);
				model.addConstr(variables_fa[i - 1][j] == 0);
				model.addConstr(variables_ha[i - 1][j] == 0);
				model.addConstr(variables_ac32[i - 1][j] == 0);
				model.addConstr(variables_ac42[i - 1][j] == 0);
				model.addConstr(variables_error[i - 1][j] == 0);
			}
			else
			{
				model.addConstr(variables_error[i - 1][j] == Pars.MED_ac32 * variables_ac32[i - 1][j] + Pars.MED_ac42 * variables_ac42[i - 1][j]);
				model.addConstr(3 * variables_fa[i - 1][j] + 2 * variables_ha[i - 1][j] + 3 * variables_ac32[i - 1][j] + 4 * variables_ac42[i - 1][j] <= variables_V[i - 1][j]);
				model.addConstr(variables_V[i][j] == variables_V[i - 1][j] - variables_ha[i - 1][j] - 2 * variables_fa[i - 1][j] + (variables_ha[i - 1][j - 1] + variables_fa[i - 1][j - 1]) - variables_ac32[i - 1][j] - 2 * variables_ac42[i - 1][j]);
			}


		}
	}
}

void generate_constraints_of_last_stage(MUOPT_Par_t Pars, vector<vector<GRBVar>> variables_V, GRBModel& model)
{
	for (int i = 0; i < variables_V[0].size() - 1; i++)
	{
		model.addConstr(variables_V[Pars.nStages][i] <= 2);
		if (i < Pars.input_patterns.size())
			model.addConstr(variables_V[Pars.nStages][i] >= 1);
	}
}

void generate_cost_of_adders(vector<vector<GRBVar>> variables_fa, vector<vector<GRBVar>> variables_ha, GRBVar& fa_num, GRBVar& ha_num,
	vector<vector<GRBVar>> variables_ac32, vector<vector<GRBVar>> variables_ac42, GRBVar& ac32_num, GRBVar& ac42_num, GRBModel& model)
{
	GRBLinExpr fs = 0;
	GRBLinExpr hs = 0;
	GRBLinExpr ac32s = 0;
	GRBLinExpr ac42s = 0;
	for (int i = 0; i < variables_fa.size(); i++)
	{
		for (int j = 0; j < variables_fa[i].size(); j++)
		{
			fs += variables_fa[i][j];
			hs += variables_ha[i][j];
			ac32s += variables_ac32[i][j];
			ac42s += variables_ac42[i][j];
			//obj += (3 * variables_fa[i][j] + 2 * variables_ha[i][j]);
		}
	}
	ha_num = model.addVar(0, INFINITY, 0, GRB_INTEGER, "hs");
	fa_num = model.addVar(0, INFINITY, 0, GRB_INTEGER, "fs");
	ac32_num = model.addVar(0, INFINITY, 0, GRB_INTEGER, "ac32s");
	ac42_num = model.addVar(0, INFINITY, 0, GRB_INTEGER, "ac42s");
	model.addConstr(fs == fa_num);
	model.addConstr(hs == ha_num);
	model.addConstr(ac32s == ac32_num);
	model.addConstr(ac42s == ac42_num);
}

void generate_total_error(vector<vector<GRBVar>> variables_error, GRBVar& E, GRBModel& model)
{
	vector<GRBVar> Es;
	for (int j = 0; j < variables_error[0].size(); j++)
		Es.push_back(model.addVar(0, INFINITY, 0, GRB_CONTINUOUS, "E_" + to_string(j)));
	for (int j = 0; j < variables_error[0].size(); j++)
	{
		GRBLinExpr columnerror = 0;
		for (int i = 0; i < variables_error.size(); i++)
			columnerror += variables_error[i][j];
		model.addConstr(Es[j] == columnerror);
	}
	E = model.addVar(0, INFINITY, 0, GRB_CONTINUOUS, "E");
	GRBLinExpr sumerror = 0;
	for (int j = 0; j < Es.size(); j++)
		sumerror += pow(2, j) * Es[j];
	model.addConstr(E == sumerror);
}
