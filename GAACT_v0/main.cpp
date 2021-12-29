#include <iostream>
#include <vector>
#include <stdlib.h>
#include "Multiplier_Tree_Optimization.h"
#include "GAACT.h"
#include "GA2HDL.h"
#include <time.h>
using namespace std;

int main(int argc, char* argv[])
{
	MUOPT_Par_t Par;
	/*int nFA, nHA, nAC42, nAC32;*/
	float Opt1_Analyzed_Error, Opt2_Analyzed_Error;
	ACT_Par_t ActPar;
	setenv("GUROBI_HOME", "/mnt/d/gurobi950_linux/linux64", 0);
	setenv("GRB_LICENSE_FILE", "/mnt/d/gurobi950_linux/linux64/gurobi.lic", 0);
	init_Pars(argc, argv, Par);
	cout << getenv("GUROBI_HOME") << endl;
	system(("mkdir -p " + Par.log_path).c_str());
	//optimization process
	try {
		/*******Set up environment******/
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);
		model.set(GRB_StringParam_LogFile, Par.log_path + Par.filename + ".log");
		//model.set(GRB_INT_PAR_MIPFOCUS, "1");
		model.set(GRB_INT_PAR_PRESOLVE, "2");
		/*model.set(GRB_DBL_PAR_TIMELIMIT, to_string(Time_Bound_s));*/
		//model.set(GRB_DBL_PAR_HEURISTICS, "0.5");
		//model.set(GRB_DBL_PAR_MIPGAP, "0.005");

		/*******Create variables******/
		//Compressor tree
		vector<vector<GRBVar>> variables_fa, variables_ha;
		vector<vector<GRBVar>> variables_V;
		vector<vector<GRBVar>> variables_ac32, variables_ac42;
		vector<vector<GRBVar>> variables_error;
		GRBVar fa_num, ha_num, ac32_num, ac42_num, E;
		generate_variables_of_multiplier(Par, variables_fa, variables_ha, variables_V, variables_ac32, variables_ac42, variables_error, model);

		/*Add general constraints*/
		//Compressor Tree
		Initial_constraints_of_V(Par, variables_V, model);
		generate_constraints_of_h_f_V(Par, variables_fa, variables_ha, variables_V, variables_ac32, variables_ac42, variables_error, model);
		generate_constraints_of_last_stage(Par, variables_V, model);

		// Set objective
		GRBLinExpr obj = 0;
		generate_cost_of_adders(variables_fa, variables_ha, fa_num, ha_num, variables_ac32, variables_ac42, ac32_num, ac42_num, model);
		generate_total_error(variables_error, E, model);
		// mode 0: constrain error & minimize area
		if (!Par.Opt_Mode)
		{
			model.addConstr(E <= Par.MED_bound);
			obj += (5.05 * fa_num + 2.66 * ha_num + 2.66 * ac32_num + 4.78 * ac42_num);
			//obj += (5 * fa_num + 2 * ha_num + 3 * ac32_num + 5 * ac42_num);
		}
		//mode 1: constrain area & minimize error
		else if (Par.Opt_Mode == 1)
		{
			model.addConstr((5.05 * fa_num + 2.66 * ha_num + 2.66 * ac32_num + 4.78 * ac42_num) <= Par.Area_bound);
			obj += E;
		}
		else
		{
			obj += (10 * fa_num + 5 * ha_num + 5 * ac32_num + 9 * ac42_num) + Par.Weight * E;
		}
		model.setObjective(obj, GRB_MINIMIZE);

		// Save problem
		model.write(Par.log_path + Par.filename + ".lp");

		// Optimize model
		model.update();

		model.optimize();

		//Save solution

		model.write(Par.log_path + Par.filename + ".sol");

		init_GA(ActPar, Par, variables_V, variables_ha, variables_fa, variables_ac32, variables_ac42);
		Opt1_Analyzed_Error = E.get(GRB_DoubleAttr_X);
		cout << "MED is:" << " "
			<< Opt1_Analyzed_Error << endl;
		cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}

	// start GA optimization
	clock_t start, finish;
	double totaltime;
	start = clock();
	ActPar.counter++;
	ACT_Man_t ActMan = initialization(ActPar);
	ActMan.Gene_RecordTable.clear();
	for (int i = 0; i < ActPar.nGenerations; i++)
	{
		ActMan.Pars.counter++;
		crossover(ActMan);
		mutation(ActMan);
		cout << "After " << i << " iterations: " << endl;
		cout << "the best gene's fitness is: " << ActMan.Populations[0].fitness << endl;
		cout << "the best gene's UCB fitness is: " << ActMan.Populations[0].fitness_UCB << endl;
		cout << "the best gene's recorded time is: " << ActMan.Gene_RecordTable_Populations[ActMan.Populations[0].usigns] << endl;
		ActMan.Gene_RecordTable.clear();
	}
	finish = clock();
	totaltime = double(finish - start) / CLOCKS_PER_SEC;
	Gene best_gene = ActMan.Populations[0];
	Opt2_Analyzed_Error = best_gene.fitness;
	/*Gene best_gene = ActMan.Populations[0];
	double fitness = best_gene.fitness;
	cout << "After " << ActMan.Pars.nGenerations << " iterations: " << endl;
	cout << "the best gene's fitness is: " << fitness << endl;*/
	// start to translate Gene to HDL
	vector<vector<column_sol_t>> ColSols = GA2ColSol(best_gene, ActMan.Pars);
	WriteToJson(ColSols, Par.log_path + Par.json_filename);
	ofstream f;
	f.open((Par.log_path + Par.log_filename).c_str());
	f << "Runtime of GA : " << totaltime << "s" << endl;
	if (!Par.Opt_Mode)
		f << "Error bound : " << Par.MED_bound << endl;
	else if (Par.Opt_Mode == 1)
		f << "Area bound : " << Par.Area_bound << endl;
	else
		f << "Weight : " << Par.Weight << endl;
	f << "Analytical error of opt1 : " << Opt1_Analyzed_Error << "\t" << "Analytical error of opt2 : " << Opt2_Analyzed_Error << endl;
	f << "============================================================" << endl;
	system("python --version");
	system(("python " + Par.json_path + "ApproxMULT_MyHDL.py " + Par.log_path + Par.json_filename + " ILP_ApproxMult " + Par.log_path + Par.log_filename + " " + to_string(Par.MULT_SIZE) + " " + to_string(Par.simulator_approach)).c_str());
	string command_yosys = "yosys -p \"read_verilog Multiplier_KoggeStone.v; synth; write_blif " + Par.log_path + "Multiplier_KoggeStone.blif;\"";
	system(command_yosys.c_str());
	string command_simulator = "simulator.out -a " + Par.cmp_path + to_string(Par.MULT_SIZE) + ".blif " + "-b " + Par.log_path + "Multiplier_KoggeStone.blif -f " + Par.log_path + Par.log_filename;
	system(command_simulator.c_str());
	f.close();
	
	return 0;
}