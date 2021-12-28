#pragma once
#ifndef __GAACT_H__
#define __GAACT_H__
#include <iostream>
#include <vector>
#include <algorithm>
#include "Multiplier_Tree_Optimization.h"
using namespace std;
typedef struct ACT_Par_t_ ACT_Par_t;
typedef struct ACT_Man_t_ ACT_Man_t;
typedef struct Gene_t Gene;
typedef struct Perm_t Perm;
//typedef vector<vector<vector<int>>> Gene;

struct ACT_Par_t_
{
	double crossoverate;
	double mutationrate;
	int nGenerations;
	int nPopulations;
	vector<vector<int>> geneShape; // size: nStages-1
	vector<vector<vector<int>>> nCompressors; // size: nStages; 
											  // {{{n22,n32,nA32,nA42}}}
	int BitWidth;
};

struct Perm_t
{
	vector<int> perm;
	Perm_t()
	{
		perm = {};
	}
	Perm_t(int num)
	{
		perm.resize(num);
		for (int i = 0; i < num; i++)
			perm[i] = i;
		if (perm.size() > 1)
			random_shuffle(perm.begin(), perm.end());
	}
	void PermShow()
	{
		for (auto it : perm)
			cout << it;
	}
};

struct Gene_t
{
	vector<Perm> g;
	int nPerms;
	//vector<vector<vector<int>>> g;
	double fitness;
	Gene_t()
	{
		g = {};
		fitness = 0.0;
	}

	Gene_t(ACT_Par_t Pars)
	{
		nPerms = Pars.geneShape.size() * Pars.geneShape[0].size();
		g.resize(nPerms);
		for (int i = 0; i < g.size(); i++)
		{
			int iStage = i / Pars.geneShape[0].size();
			int iColumn = i % Pars.geneShape[0].size();
			int PermSize = Pars.geneShape[iStage][iColumn];
			g[i] = Perm(PermSize);
		}
	}

	void show(ACT_Par_t Pars)
	{
		int nColumn = Pars.geneShape[0].size();
		for (int i = 0; i < Pars.geneShape.size(); i++)
		{
			for (int j = 0; j < Pars.geneShape[0].size(); j++)
			{
				g[nColumn * i + j].PermShow();
				cout << "|";
			}
			cout << endl;
		}
	}

	void addGene(int Pos, Perm permu)
	{
		g.insert(g.begin() + Pos, permu);
		nPerms++;
	}
};

struct ACT_Man_t_
{
	ACT_Par_t Pars;
	vector<Gene> Populations;
	int nGenes;
	ACT_Man_t_()
	{
		nGenes = 0;
	}
	ACT_Man_t_(ACT_Par_t p)
	{
		nGenes = 0;
		Pars = p;
		Populations.resize(p.nPopulations*10);
	}
	/*vector<double> Fitness;*/
};

static inline vector<double> Ex22(double a, double b) { return { a + b - 2 * a * b,a * b }; }
static inline vector<double> Ex32(double a, double b, double c) { return { a + b + c - 2 * (a * b + a * c + b * c) + 4 * a * b * c,a * b + a * c + b * c - 2 * a * b * c }; }
static inline vector<double> Ax32(double a, double b, double c) { return { a + b - a * b,c + a * b - a * b * c }; }
static inline vector<double> Ax42(double a, double b, double c, double d) {
	return { a + b + c * d + 2 * a * d + a * b - 2 * a * b * d - b * c * d - a * c * d + a * b * c * d,
			 c + d + a * b - c * d - a * b * c - a * b * d + a * b * c * d };
		   }
static inline double mAx32(double a, double b, double c) { return a * b * c; }
static inline double mAx42(double a, double b, double c, double d) { return a * b * c + b * c * d + 2 * a * b * d * (1 - c); }
double ComputeAnalyticMED(Gene& g, ACT_Par_t Pars);
void init_GA(ACT_Par_t& Pars, MUOPT_Par_t MuPars, vector<vector<GRBVar>> variables_V,
	vector<vector<GRBVar>> variables_ha, vector<vector<GRBVar>> variables_fa,
	vector<vector<GRBVar>> variables_ac32, vector<vector<GRBVar>> variables_ac42);
ACT_Man_t initialization(ACT_Par_t Pars);
void crossover(ACT_Man_t& ActMan);
void mutation(ACT_Man_t& ActMan);
#endif // !__GAACT_H__

