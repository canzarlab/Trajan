/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#ifndef ANTICHAIN_CONSTRAINT_H
#define ANTICHAIN_CONSTRAINT_H

#include "Constraint.h"
#include <mutex>
#include <queue>

class AntichainConstraint : Constraint
{
public:
    AntichainConstraint(vector<ET>& Triplets, Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp);
    int AddTriplets(int nr_rows);

private:
    void RunParallel();
    void AntichainJob(int id);
    double MaxFlow(vi& D, vvd& R);
    double Push(int x, double flow, vvd& R, vi& D);
    bool LevelGraph(vi& D, vvd& R);
    double Reset(vn& P, vvd& R);
    void Antichain(vn& P, vvd& R);

    vvi& G;
    int ncr, nr_rows, S, T, Z, SZ;
    mutex qmutex;
    vb B;
    int pi;
};

#endif
