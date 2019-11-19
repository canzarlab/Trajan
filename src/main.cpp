/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include <iostream>
#include <string>
#include "Timer.h"
#include "Greedy.h"
#include "LP.h"
#include "BnB.h"
#include "BnG.h"
#include "LPInt.h"
#include "LPCP.h"
#include "LPFInt.h"
#include "Similarity.h"
#include "Parallel.h"
#include "read_csv.h"

// global varialbe for cost matrix name
std::string costMatrixFileName;
// global variable for score name
std::string outScoreFile;

int Solver::cf;
bool Solver::tt;

void EditDist(Graph& g1, Graph& g2, string & fileName, string & outSolution, std::map<string, string> & t1Label2Node, std::map<string, string> & t2Label2Node, string & outScore);

Solver* MakeSolver(Graph& t1, Graph& t2, int argc, char** argv, std::map<string, string> & t1Label2Node, std::map<string, string> & t2Label2Node)
{
    int s = stoi(argv[argc - 1]);
    Solver::cf = stoi(argv[4 + (argc == 9) + 2 * (argc == 12)]);
    Solver::tt = argc == 12;
    string d = argv[5 + (argc == 9) + 2 * (argc == 12)];
    double k = stod(argv[6 + (argc == 9) + 2 * (argc == 12)]);
    double c = (s != 2) ? 0 : stod(argv[8 + 2 * (argc == 12)]);
    var_eps = (argc == 9) ? 0 : stod(argv[7 + 2 * (argc == 12)]);
    bool dag = argc == 9;

    assert(LP::cf >= 0 && LP::cf <= 2);
    assert(d == "j" || d == "s" || d == "e");
    assert(s >= 0 && s <= 9);

    string outFileName = argv[5];

    if (s == 0)
        return new Greedy(t1, t2, d, k, dag);
    else if (s == 1)
        return new LP(t1, t2, d, k, dag);
    else if (s == 2)
        EditDist(t1, t2, costMatrixFileName, outFileName, t1Label2Node, t2Label2Node, outScoreFile);
    else if (s == 3)
        return new LPCP(t1, t2, d, k, dag);
    else if (s == 4)
        return new DFBnBSolver(t1, t2, d, k, dag); // was test, to delete
    else if (s == 5)
        return new LPInt(t1, t2, d, k, dag);
    else if (s == 6)
        return new LPFInt(t1, t2, d, k, dag);
    else if (s == 7)
        return new BFBnBSolver(t1, t2, d, k, dag);
    else if (s == 8)
        return new DFBnBSolver(t1, t2, d, k, dag);
    else if (s == 9)
        return new HybridBnBSolver(t1, t2, d, k, dag);
    return nullptr;
}

Graph* MakeDAG(const char* f1, const char* f2, int s)
{
    return ((s == 0) ? new DAG(f1, f2) : new LDAG(f1, f2))->Init();
}

pair<Graph*, Graph*> MakeGraphs(int argc, char** argv)
{
    if (argc == 10)
        return {new Tree(argv[1]), new Tree(argv[2])};
    else if (argc == 12)
        return {new Tree(argv[1], argv[2]), new Tree(argv[3], argv[4])};

    int s = stoi(argv[argc - 1]);
    return {MakeDAG(argv[1], argv[2], s), MakeDAG(argv[3], nullptr, s)};
}

int main(int argc, char** argv)
{
    if (argc < 9 || argc > 14)
    {
        cout << "tree usage: " << argv[0] << " <filename.newick> <filename.newick> <align> <constraints> <weightfunc> <k> <vareps> <coneps> <solver>" << endl;
        cout << "dag usage: " << argv[0] << " <yeastnet> <mapping> <go> <align> <constraints> <weightfunc> <k> <solver>" << endl;
        cout << "tree usage (2): " << argv[0] << " <tree> <map> <tree> <map> <align> <constraints> <weightfunc> <k> <vareps> <coneps> <solver>" << endl;
        return EXIT_FAILURE;
    }
    // change the parser for cell hali inputs.
    if (argc == 14){
        costMatrixFileName = argv[5];
        outScoreFile = argv[6];
        argc = argc-2;
        for (int i = 5; i < 15; ++i){
            argv[i] = argv[i+2];
        }
    }
    if (argc == 13){
        costMatrixFileName = argv[5];
        argc = argc-1;
        for (int i = 5; i < 14; ++i){
            argv[i] = argv[i+1];
        }
    }

    Timer T;
    T.start();
    Graph *t1, *t2;
    if (stoi(argv[argc - 1]) == 2){
        // create two new tree.
        string t1fName   = argv[1];
        string t1MapfName = argv[2];
        string t2fName   = argv[3];
        string t2MapfName = argv[4];
        CSVReader t1_treeReader(t1fName, " ");
        vector<vector<string>> t1_tree = t1_treeReader.getStringData();
        CSVReader t1_mapReader(t1MapfName, " ");
        vector<vector<string>> t1_map = t1_mapReader.getStringData();
        
        CSVReader t2_treeReader(t2fName, " ");
        vector<vector<string>> t2_tree = t2_treeReader.getStringData();
        CSVReader t2_mapReader(t2MapfName, " ");
        vector<vector<string>> t2_map = t2_mapReader.getStringData();
        
        std::map<string, string> t1Node2Label;
        std::map<string, string> t1Label2Node;
        std::map<string, string> t2Node2Label;
        std::map<string, string> t2Label2Node;
        for (const auto &u: t1_map){
            t1Label2Node[u[0]] = u[1];
            t1Node2Label[u[1]] = u[0];
        }
        for (const auto &u: t2_map){
            t2Label2Node[u[0]] = u[1];
            t2Node2Label[u[1]] = u[0];
        }
        // change the nodes for tree.
        for (auto & u: t1_tree){
            u[0] = t1Node2Label[u[0]];
            u[1] = t1Node2Label[u[1]];
        }
        for (auto & u: t2_tree){
            u[0] = t2Node2Label[u[0]];
            u[1] = t2Node2Label[u[1]];
        }
        
        string redundantT1 = "redundant"+outScoreFile+argv[5]+"doNotDeleteT1.ooo";
        string redundantT2 = "redundant"+outScoreFile+argv[5]+"doNotDeleteT2.ooo";
        ofstream t1NewTree(redundantT1);
        for (auto & u: t1_tree){
            for (auto & v: u){
                t1NewTree << v << " ";
            }
            t1NewTree << endl;
        }
        t1NewTree.close();
        
        ofstream t2NewTree(redundantT2);
        for (auto & u: t2_tree){
            for (auto & v: u){
                t2NewTree << v << " ";
            }
            t2NewTree << endl;
        }
        t2NewTree.close();
        
        // argv[1] = const_cast<char*>(redundantT1);
        // argv[3] = const_cast<char*>(redundantT2);
        cout << redundantT1 << ", " << redundantT2 << endl;
        argv[1] = const_cast<char*>(redundantT1.c_str());
        argv[3] = const_cast<char*>(redundantT2.c_str());
        tie(t1, t2) = MakeGraphs(argc, argv);
        // remove files 
        // remove(redundantT1);
        // remove(redundantT2);
        Solver* solver = MakeSolver(*t1, *t2, argc, argv, t1Label2Node, t2Label2Node);
        if (solver) solver->Solve(argv[3 + (argc == 9) + 2 * (argc == 12)]);
        delete solver;
        delete t1;
        delete t2;
    }
    else
    {
        tie(t1, t2) = MakeGraphs(argc, argv);
        if (stoi(argv[argc - 1]) > 9)
        {
            int s = stoi(argv[argc - 1]);
            Solver::cf = stoi(argv[4 + (argc == 9) + 2 * (argc == 12)]);
            Solver::tt = argc == 12;
            string d = argv[5 + (argc == 9) + 2 * (argc == 12)];
            double k = stod(argv[6 + (argc == 9) + 2 * (argc == 12)]);
            double c = (s != 2) ? 0 : stod(argv[8 + 2 * (argc == 12)]);
            var_eps = (argc == 9) ? 0 : stod(argv[7 + 2 * (argc == 12)]);

            assert(LP::cf >= 0 && LP::cf <= 2);
            assert(d == "j" || d == "s" || d == "e");
            ParallelSolver(*t1, *t2, d, k, argc == 9, s - 9).Solve(argv[3 + (argc == 9) + 2 * (argc == 12)], outScoreFile);
        }
        else
        {
            std::map<string, string> redundantMap;
            Solver* solver = MakeSolver(*t1, *t2, argc, argv, redundantMap, redundantMap);
            if (solver) solver->Solve(argv[3 + (argc == 9) + 2 * (argc == 12)], outScoreFile);
            delete solver;
            delete t1;
            delete t2;
        }
    }
    
    T.stop();
    clog << "TIME: " << T.secs() << " secs" << endl;
}
