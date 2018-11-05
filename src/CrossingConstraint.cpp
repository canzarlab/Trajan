/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include "CrossingConstraint.h"
#include <thread>

CrossingConstraint::CrossingConstraint(vector<ET>& Triplets, Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp) : Constraint(Triplets, t1, t2, K, x, swp)
{
    DP.resize(t1.GetNumNodes());
    PA.resize(t1.GetNumNodes());
    for (auto& v : DP)
        v.resize(t2.GetNumNodes());

    vb C(t1.GetNumNodes());
    DFSLeft(t1.GetRoot(), C);
    Q.push(t1.GetRoot());
    RunParallel();
}

int CrossingConstraint::AddTriplets(int nr_rows)
{
    int ncr = 0;
    for (auto node : t1.L)
    {
        vii P;
        Reconstruct(P, node, t2.GetRoot());

        double sum = 0;
        for (auto k : P)
            sum += GetWeight(k.first, k.second);

        if (sum - EPS > 1)
            AddConstraint(nr_rows + ncr++, P);
    }
    return ncr;
}

void CrossingConstraint::RunParallel()
{
    const size_t NR_THREADS = min(::NR_THREADS, size_t(64));
    ax = (NR_THREADS < 64 ? (uint64_t(1) << NR_THREADS) - 1 : ~0);
    vector<thread> vt;
    for (int i = 0; i < NR_THREADS; ++i)
        vt.emplace_back(&CrossingConstraint::CrossingJob, this, i);

    for (int i = 0; i < NR_THREADS; ++i)
        vt[i].join();
}

void CrossingConstraint::CrossingJob(int i)
{
    while (ax)
    {
        newick_node* node = nullptr;
        {
            lock_guard<mutex> g(qmutex);
            if (Q.empty())
            {
                ax &= ~(1 << i);
                continue;
            }
            else
                ax |= 1 << i;
            node = Q.front();
            Q.pop();
        }
        DFSRight(t2.GetRoot(), node);
        lock_guard<mutex> g(qmutex);
        for (newick_child* child = node->child; child; child = child->next)
            if (--PA[child->node->taxoni] == 0)
                Q.push(child->node);
    }
}

pair<newick_node*, double> CrossingConstraint::GetMaxPC(newick_node* nodel, newick_child* noder, bool s)
{
    double mx = 0;
    newick_node* mc = nullptr;
    for (newick_child* pc = noder; pc; pc = pc->next)
    {
        double cw = GetDP(nodel, pc->node, s);
        if (cw >= mx)
        {
            mx = cw;
            mc = pc->node;
        }
    }
    return make_pair(mc, mx);
}

inline double& CrossingConstraint::GetDP(newick_node* nodel, newick_node* noder, bool s = false)
{
    return s ? DP[noder->taxoni][nodel->taxoni] : DP[nodel->taxoni][noder->taxoni];
}

inline pair<newick_node*, double> CrossingConstraint::GetMaxChild(newick_node* nodel, newick_node* noder)
{
    return GetMaxPC(nodel, noder->child, false);
}

inline pair<newick_node*, double> CrossingConstraint::GetMaxParent(newick_node* nodel, newick_node* noder)
{
    return GetMaxPC(nodel, noder->parent, true);
}

void CrossingConstraint::DFSLeft(newick_node* node, vb& C)
{
    C[node->taxoni] = true;
    for (newick_child* child = node->child; child; child = child->next)
    {
        PA[child->node->taxoni]++;
        if (!C[child->node->taxoni])
            DFSLeft(child->node, C);
    }
}

double CrossingConstraint::DFSRight(newick_node* node, newick_node* nodel)
{
    double mx = 0;
    for (newick_child* child = node->child; child; child = child->next)
        mx = max(mx, DFSRight(child->node, nodel));
    mx = max(mx, GetMaxParent(node, nodel).second);
    return GetDP(nodel, node) = mx + GetWeight(nodel, node);
}

void CrossingConstraint::Reconstruct(vii& P, newick_node* nodel, newick_node* noder)
{
    double pw, cw;
    newick_node *child, *parent;
    tie(child, cw) = GetMaxChild(nodel, noder);
    tie(parent, pw) = GetMaxParent(noder, nodel);
    P.emplace_back(nodel->taxoni, noder->taxoni);
    if (nodel->parent && (!child || pw > cw))
        Reconstruct(P, parent, noder);
    else if (child && (!parent || cw >= pw))
        Reconstruct(P, nodel, child);
    else
        assert(!parent && !child);
}
