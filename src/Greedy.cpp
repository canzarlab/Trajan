/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include "Greedy.h"
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <numeric>
using namespace std::placeholders;

Greedy::Greedy(Graph& t1, Graph& t2, string d, double k, bool dag) : Solver(t1, t2, d, k, dag), A(t1.GetNumNodes(), vd(t2.GetNumNodes()))
{
    vb P(t1.GetNumNodes());
    DFSLeft(t1.GetRoot(), P, [&](newick_node* nodel, newick_node* noder, double w)
    {
        if (w != 0.0 && (dag || nodel->parent) && (dag || nodel->child) && (dag || noder->parent) && (dag || noder->child))
            E.emplace_back(nodel->taxoni, noder->taxoni, w);
    });
    sort(E.begin(), E.end(), [](const iid& a, const iid& b){return get<2>(a) > get<2>(b);});
}

Greedy::Greedy(Graph& t1, Graph& t2, string d, double k, vvi& K, map<size_t, bool>& N) : Solver(t1, t2, d, k, false), A(t1.GetNumNodes(), vd(t2.GetNumNodes()))
{
    vb P(t1.GetNumNodes());
    DFSLeft(t1.GetRoot(), P, [&](newick_node* nodel, newick_node* noder, double w)
    {
        if (w != 0.0 && (dag || nodel->parent) && (dag || nodel->child) && (dag || noder->parent) && (dag || noder->child))
        {
            auto e = N.find(K[nodel->taxoni][noder->taxoni]);
            if (e == N.end())
                E.emplace_back(nodel->taxoni, noder->taxoni, w);
            else if (e->second)
            {
                this->M.emplace_back(nodel->taxoni, noder->taxoni, w);
                A[nodel->taxoni][noder->taxoni] = w;
            }
        }
    });
    sort(E.begin(), E.end(), [](const iid& a, const iid& b){return get<2>(a) > get<2>(b);});
}

void Greedy::Solve(string filename, string outScoreFile)
{
    for (const iid& e : E)
        if (all_of(M.begin(), M.end(), bind(&Greedy::CC, this, _1, cref(e))))
            M.push_back(e), A[get<0>(e)][get<1>(e)] = get<2>(e);
    if (filename != "") WriteSolution(filename);
}

void Greedy::WriteSolution(string fileName)
{
    ofstream sol_file(fileName);
    double weight = 0.0;
    vn n1 = t1.GetNodes(), n2 = t2.GetNodes();
    for (const iid& e : M)
        sol_file << n1[get<0>(e)]->taxon << " " << n2[get<1>(e)]->taxon << " 1\n", weight += get<2>(e);
    PrintScore(weight);
}

double Greedy::GetSolution()
{
    return accumulate(M.begin(), M.end(), 0., [](double a, iid& e) { return a + get<2>(e); });
}

bool Greedy::CC(const iid& a, const iid& b) const
{
    int i = get<0>(a), j = get<0>(b);
    int x = get<1>(a), y = get<1>(b);
    return IsNotInConflict(i, j, x, y);
}

void Greedy::GetSolution(vvi& K, Vector& v, double& d)
{
    double weight = 0.0;
    v = Vector::Zero(v.size());
    for (const iid& e : M)
        if (K[get<0>(e)][get<1>(e)] != -1)
            weight += get<2>(e), v(K[get<0>(e)][get<1>(e)]) = 1;
    d = -weight;
}
