/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include "Graph.h"
#include "AntichainConstraint.h"
#include <thread>

AntichainConstraint::AntichainConstraint(vector<ET>& Triplets, Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp) : Constraint(Triplets, t1, t2, K, x, swp), G(((LDAG*)&t2)->G), ncr(0), B(t1.GetNumNodes()), pi(0)
{
    Z = t2.GetNumNodes();
    SZ = Z * 2 + 2;
    S = SZ - 2;
    T = SZ - 1;
}

int AntichainConstraint::AddTriplets(int nr_rows)
{
    this->nr_rows = nr_rows;
    RunParallel();
    return ncr;
}

void AntichainConstraint::RunParallel()
{
    vector<thread> vt;
    for (int i = 0; i < NR_THREADS; ++i)
        vt.emplace_back(&AntichainConstraint::AntichainJob, this, i);

    for (int i = 0; i < NR_THREADS; ++i)
        vt[i].join();
}

void AntichainConstraint::AntichainJob(int id)
{
    LDAG &g1 = (LDAG&)t1, &g2 = (LDAG&)t2;
    while (true)
    {
        int i;
        {
            lock_guard<mutex> g(qmutex);
            if (pi == g1.P.size())
                break;
            i = pi++;
            if (all_of(g1.P[i].begin(), g1.P[i].end(), [&](newick_node* node){return B[node->taxoni];}))
                continue;
        }
        Antichain(g1.P[i], g2.R[id]);
    }
}

double AntichainConstraint::Push(int x, double flow, vvd& R, vi& D)
{
    if (x == T)
        return flow;
    for (int y : G[x])
        if (R[x][y] > 0 && D[y] == D[x] + 1)
            if (double f = Push(y, min(flow, R[x][y]), R, D))
                return R[x][y] -= f, R[y][x] += f, f;
    return D[x] = -1, 0;
}

bool AntichainConstraint::LevelGraph(vi& D, vvd& R)
{
    queue<int> W;
    W.push(S);
    D[S] = 0;
    while (!W.empty())
    {
        int k = W.front(); W.pop();
        if (k == T) return true;
        for (int x : G[k])
            if (D[x] == -1 && R[k][x] > 0)
                D[x] = D[k] + 1, W.push(x);
    }
    return false;
}

double AntichainConstraint::MaxFlow(vi& D, vvd& R)
{
    double flow = 0;
    while (LevelGraph(D, R))
    {
        while (double f = Push(S, INF, R, D))
            flow += f;
        fill(D.begin(), D.end(), -1);
    }
    return flow;
}

double AntichainConstraint::Reset(vn& P, vvd& R)
{
    double sum = 0;
    for (int i = 0; i < Z; ++i)
    {
        R[S][i] = R[i][S] = R[T][i + Z] = 0;
        for (int j : G[i + Z])
            R[i + Z][j] = 0;

        for (newick_node* nodel : P)
        {
            R[S][i] += GetWeight(nodel->taxoni, i);
            R[i + Z][T] += GetWeight(nodel->taxoni, i);
        }
        sum += R[S][i];
    }
    return sum;
}

void AntichainConstraint::Antichain(vn& P, vvd& R)
{
    vi D(SZ, -1);
    double max = Reset(P, R);
    double flow = MaxFlow(D, R);
    if (max - flow <= 1 + EPS)
        return;

    vii PN;
    for (int i = 0; i < Z; ++i)
        if (D[i] != -1 && D[i + Z] == -1)
            for (newick_node* nodel : P)
                PN.emplace_back(nodel->taxoni, i);

    lock_guard<mutex> g(qmutex);
    if (!all_of(P.begin(), P.end(), [&](newick_node* node){return B[node->taxoni];}))
    {
        AddConstraint(nr_rows + ncr++, PN);
        for (newick_node* node : P)
            B[node->taxoni] = true;
    }
}
