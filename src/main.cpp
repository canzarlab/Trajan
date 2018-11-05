/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
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

#include <iostream>

int Solver::cf;
bool Solver::tt;

void EditDist(Graph& g1, Graph& g2);

Solver* MakeSolver(Graph& t1, Graph& t2, int argc, char** argv)
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
    assert(d == "j" || d == "s");
    assert(s >= 0 && s <= 9);

    if (s == 0)
        return new Greedy(t1, t2, d, k, dag);
    else if (s == 1)
        return new LP(t1, t2, d, k, dag);
    else if (s == 2)
        EditDist(t1, t2);
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
    if (argc < 9 || argc > 12)
    {
        cout << "tree usage: " << argv[0] << " <filename.newick> <filename.newick> <align> <constraints> <weightfunc> <k> <vareps> <coneps> <solver>" << endl;
        cout << "dag usage: " << argv[0] << " <yeastnet> <mapping> <go> <align> <constraints> <weightfunc> <k> <solver>" << endl;
        cout << "tree usage (2): " << argv[0] << " <tree> <map> <tree> <map> <align> <constraints> <weightfunc> <k> <vareps> <coneps> <solver>" << endl;
        return EXIT_FAILURE;
    }

    Timer T;
    T.start();
    Graph *t1, *t2;
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
        assert(d == "j" || d == "s");

        ParallelSolver(*t1, *t2, d, k, argc == 9, s - 9).Solve(argv[3 + (argc == 9) + 2 * (argc == 12)]);
    }
    else
    {
        Solver* solver = MakeSolver(*t1, *t2, argc, argv);
        if (solver) solver->Solve(argv[3 + (argc == 9) + 2 * (argc == 12)]);
        delete solver;
        delete t1;
        delete t2;
    }
    T.stop();
    clog << "TIME: " << T.secs() << " secs" << endl;
}
