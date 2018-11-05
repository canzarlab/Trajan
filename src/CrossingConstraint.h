/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#ifndef CROSSING_CONSTRAINT_H
#define CROSSING_CONSTRAINT_H

#include "Constraint.h"
#include <mutex>
#include <queue>
#include <atomic>
#include <cstdint>

class CrossingConstraint : Constraint
{
public:
    CrossingConstraint(vector<ET>& Triplets, Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp);
    
    int AddTriplets(int nr_rows);
private:
    pair<newick_node*, double> GetMaxPC(newick_node* nodel, newick_child* noder, bool s);
    inline pair<newick_node*, double> GetMaxChild(newick_node* nodel, newick_node* noder);
    inline pair<newick_node*, double> GetMaxParent(newick_node* nodel, newick_node* noder);

    void RunParallel();
    void CrossingJob(int i);
    void DFSLeft(newick_node* node, vb& C);
    double DFSRight(newick_node* node, newick_node* nodel);
    void Reconstruct(vii& P, newick_node* nodel, newick_node* noder);
    inline double& GetDP(newick_node* nodel, newick_node* noder, bool s);

    atomic<uint64_t> ax;
    mutex qmutex;
    queue<newick_node*> Q;
    vi PA;
    vvd DP;
};

#endif
