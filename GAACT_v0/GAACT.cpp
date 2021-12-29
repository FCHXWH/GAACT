#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <assert.h>
#include <algorithm>
#include "GAACT.h"
using namespace std;
//void GA_GeneSort(ACT_Man_t& ActMan, Gene gnew)
//{
//    ActMan.Populations[ActMan.nGenes] = gnew;
//    for (int i = ActMan.nGenes - 1; i >= 0; i--)
//    {
//        if (gnew.fitness >= ActMan.Populations[i].fitness)
//            break;
//        ActMan.Populations[i + 1] = ActMan.Populations[i];
//        ActMan.Populations[i] = gnew;
//    }
//    if (ActMan.nGenes < ActMan.Pars.nPopulations)
//        ActMan.nGenes++;
//    /*for (int i = 0; i < ActMan.nGenes - 1; i++)
//        assert(ActMan.Populations[i].fitness <= ActMan.Populations[i + 1].fitness);*/
//}

void GA_GeneSort(ACT_Man_t& ActMan, Gene gnew)
{
    ActMan.Populations[ActMan.nGenes] = gnew;
    for (int i = ActMan.nGenes - 1; i >= 0; i--)
    {
        if (gnew.fitness_UCB >= ActMan.Populations[i].fitness_UCB)
            break;
        ActMan.Populations[i + 1] = ActMan.Populations[i];
        ActMan.Populations[i] = gnew;
    }
    if (ActMan.nGenes < ActMan.Pars.nPopulations)
        ActMan.nGenes++;
    /*for (int i = 0; i < ActMan.nGenes - 1; i++)
        assert(ActMan.Populations[i].fitness <= ActMan.Populations[i + 1].fitness);*/
}

void init_GA(ACT_Par_t& Pars,MUOPT_Par_t MuPars, vector<vector<GRBVar>> variables_V, 
    vector<vector<GRBVar>> variables_ha, vector<vector<GRBVar>> variables_fa, 
    vector<vector<GRBVar>> variables_ac32, vector<vector<GRBVar>> variables_ac42)
{
    Pars.mutationrate = 0.05;
    Pars.crossoverate = 0.8;
    Pars.nPopulations = 100;
    Pars.nGenerations = 300;
    Pars.BitWidth = MuPars.MULT_SIZE;
    Pars.geneShape.resize(variables_ha.size() - 1);
    Pars.nCompressors.resize(variables_ha.size());
    for (int i = 0; i < variables_ha.size(); i++)
    {
        for (int j = 0; j < variables_ha[i].size(); j++)
        {
            if (i != variables_ha.size() - 1)
                Pars.geneShape[i].push_back(int(variables_V[i + 1][j].get(GRB_DoubleAttr_X)));
            Pars.nCompressors[i].push_back({ int(variables_ha[i][j].get(GRB_DoubleAttr_X)),int(variables_fa[i][j].get(GRB_DoubleAttr_X)),
                int(variables_ac32[i][j].get(GRB_DoubleAttr_X)),int(variables_ac42[i][j].get(GRB_DoubleAttr_X)) });
        }
    }
}

double ComputeAnalyticMED(Gene& g, ACT_Par_t Pars)
{
    vector<vector<double>> z_cur,z_next;
    vector<vector<double>> y;
    g.fitness = 0.0;
    int nStages = Pars.geneShape.size() + 1;
    int nColumns = Pars.geneShape[0].size();
    // initialize probs as 0.25
    z_cur.resize(nColumns); z_next.resize(nColumns); y.resize(nColumns);
    for (int j = 0; j < nColumns; j++)
        if (j < Pars.BitWidth)
            z_cur[j] = vector<double>(j + 1, 0.25);
        else
            z_cur[j] = vector<double>(((2 * Pars.BitWidth - j - 1) > 0) ? (2 * Pars.BitWidth - j - 1) : 0, 0.25);
    // computation process
    for (int i = 0; i < nStages; i++)
    {
        // generate this stage's y
        for (int j = 0; j < nColumns; j++)
        {
            int startIndex = 0, endIndex;
            // derive all ys
            // exact 2:2 compressors
            for (int k = 0; k < Pars.nCompressors[i][j][0]; k++)
            {
                endIndex = startIndex + 2;
                vector<double> outs = Ex22(z_cur[j][startIndex], z_cur[j][startIndex + 1]);
                y[j].push_back(outs[0]);
                y[j + 1].push_back(outs[1]);
                startIndex = endIndex;
            }
            // exact 3:2 compressors
            for (int k = 0; k < Pars.nCompressors[i][j][1]; k++)
            {
                endIndex = startIndex + 3;
                vector<double> outs = Ex32(z_cur[j][startIndex], z_cur[j][startIndex + 1], z_cur[j][startIndex + 2]);
                y[j].push_back(outs[0]);
                y[j + 1].push_back(outs[1]);
                startIndex = endIndex;
            }
            // aprroximate 3:2 compressors
            for (int k = 0; k < Pars.nCompressors[i][j][2]; k++)
            {
                endIndex = startIndex + 3;
                vector<double> outs = Ax32(z_cur[j][startIndex], z_cur[j][startIndex + 1], z_cur[j][startIndex + 2]);
                g.fitness += pow(2, j) * mAx32(z_cur[j][startIndex], z_cur[j][startIndex + 1], z_cur[j][startIndex + 2]);
                y[j].push_back(outs[0]);
                y[j].push_back(outs[1]);
                startIndex = endIndex;
            }
            // aprroximate 4:2 compressors
            for (int k = 0; k < Pars.nCompressors[i][j][3]; k++)
            {
                endIndex = startIndex + 4;
                vector<double> outs = Ax42(z_cur[j][startIndex], z_cur[j][startIndex + 1], z_cur[j][startIndex + 2], z_cur[j][startIndex + 3]);
                g.fitness += pow(2, j) * mAx42(z_cur[j][startIndex], z_cur[j][startIndex + 1], z_cur[j][startIndex + 2], z_cur[j][startIndex + 3]);
                y[j].push_back(outs[0]);
                y[j].push_back(outs[1]);
                startIndex = endIndex;
            }
            // dummy compressors
            for (int k = startIndex; k < z_cur[j].size(); k++)
                y[j].push_back(z_cur[j][k]);
        }
        // generate this stage's z according to gene g
        if (i == nStages - 1)
            continue;
        z_next = y;
        for (int j = 0; j < nColumns; j++)
        {
            int iPerm = i * nColumns + j;
            Perm p = g.g[iPerm];
            for (int k = 0; k < p.perm.size(); k++)
                z_next[j][p.perm[k]] = y[j][k];
        }
        z_cur = z_next;
    }
    g.fitness_UCB = g.fitness - sqrt(double(2) * double(log(Pars.counter)) / g.nRecord);
}

ACT_Man_t initialization(ACT_Par_t Pars)
{
    ACT_Man_t ActMan(Pars);
    ActMan.Pars.counter++;
    // initialize the first population
    for (int i = 0; i < ActMan.Pars.nPopulations; i++)
    {
        Gene g_new = Gene(ActMan.Pars);
        if (!ActMan.Gene_RecordTable[g_new.usigns])
        {
            ActMan.Gene_RecordTable[g_new.usigns] = 1;
            ActMan.Gene_RecordTable_Populations[g_new.usigns] += 1;
            g_new.nRecord = ActMan.Gene_RecordTable_Populations[g_new.usigns];
            ComputeAnalyticMED(g_new, Pars);
            GA_GeneSort(ActMan, g_new);
        }
        /*ActMan.Populations[i].show(Pars);
        cout << "Fitness is : " << ActMan.Populations[i].fitness << endl;*/
    }
    return ActMan;
}

void crossover(ACT_Man_t& ActMan)
{
    // do crossover between 2*i and 2*i+1 genes
    for (int i = 0; i < ActMan.Pars.nPopulations / 2; i++)
    {
        int prob = rand() % 10;
        if (prob > int(10 * ActMan.Pars.crossoverate))
            continue;
        Gene father = ActMan.Populations[2 * i], mother = ActMan.Populations[2 * i + 1];
        // sub-gene for one column
        for (int j = 0; j < father.g.size(); j++)
        {
            vector<int> permu1 = father.g[j].perm, permu2 = mother.g[j].perm;
            if (permu1.size() < 2)
                continue;
            int n1 = rand() % permu1.size();
            int tmp = rand() % (permu1.size() - 1);
            int n2 = (tmp < n1) ? tmp : tmp + 1;
            int startIndex = min(n1, n2), endIndex = max(n1, n2);
            for (int k = 0; k < permu1.size(); k++)
            {
                auto it = find(mother.g[j].perm.begin() + startIndex, mother.g[j].perm.begin() + endIndex + 1, permu1[k]);
                if (it != mother.g[j].perm.begin() + endIndex + 1)
                {
                    permu1.erase(permu1.begin() + k);
                    k--;
                }
            }
            for (int k = 0; k < permu2.size(); k++)
            {
                auto it = find(father.g[j].perm.begin() + startIndex, father.g[j].perm.begin() + endIndex + 1, permu2[k]);
                if (it != father.g[j].perm.begin() + endIndex + 1)
                {
                    permu2.erase(permu2.begin() + k);
                    k--;
                }
            }    
            for (int k = startIndex; k <= endIndex; k++)
            {
                permu1.insert(permu1.begin() + k, mother.g[j].perm[k]);
                permu2.insert(permu2.begin() + k, father.g[j].perm[k]);
            }

            father.g[j].perm = permu1;
            mother.g[j].perm = permu2;
        
        }
        father.Gene_Sign(); mother.Gene_Sign();
        // compute the new two genes' fitness
        if (!ActMan.Gene_RecordTable[father.usigns])
        {
            ActMan.Gene_RecordTable[father.usigns] = 1;
            ActMan.Gene_RecordTable_Populations[father.usigns] += 1;
            father.nRecord = ActMan.Gene_RecordTable_Populations[father.usigns];
            ComputeAnalyticMED(father, ActMan.Pars);
            GA_GeneSort(ActMan, father);
        }
        if (!ActMan.Gene_RecordTable[mother.usigns])
        {
            ActMan.Gene_RecordTable[mother.usigns] = 1;
            ActMan.Gene_RecordTable_Populations[mother.usigns] += 1;
            mother.nRecord = ActMan.Gene_RecordTable_Populations[mother.usigns];
            ComputeAnalyticMED(mother, ActMan.Pars);
            GA_GeneSort(ActMan, mother);
        }
    }
}

void mutation(ACT_Man_t& ActMan)
{
    for (int i = 0; i < ActMan.Pars.nPopulations; i++)
    {
        Gene g = ActMan.Populations[i];
        int n3 = rand() % 20;
        if (n3 == 2)
        {
            for (int j = 0; j < g.g.size(); j++)
            {
                vector<int>& permu = g.g[j].perm;
                if (permu.size() < 2)
                    continue;
                int b1 = rand() % permu.size(), b2 = rand() % permu.size();
                int tmp = permu[b1];
                permu[b1] = permu[b2]; permu[b2] = tmp;
            }
        }
        g.Gene_Sign();
        if (!ActMan.Gene_RecordTable[g.usigns])
        {
            ActMan.Gene_RecordTable[g.usigns] = 1;
            ActMan.Gene_RecordTable_Populations[g.usigns] += 1;
            g.nRecord = ActMan.Gene_RecordTable_Populations[g.usigns];
            ComputeAnalyticMED(g, ActMan.Pars);
            GA_GeneSort(ActMan, g);
        }
    }
}