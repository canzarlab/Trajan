/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#ifndef GREEDY_H
#define GREEDY_H

#include "Solver.h"
#include "Geno.h"
#include <tuple>

typedef tuple<int, int, double> iid;
typedef vector<iid> viid;

class Greedy : public Solver
{
public:
    Greedy(Graph& t1, Graph& t2, string d, double k, bool dag);    
    Greedy(Graph& t1, Graph& t2, string d, double k, vvi& K, map<size_t, bool>& M);

    virtual void Solve(string filename) override;
    void WriteSolution(string fileName) override;
    double GetSolution();
    void GetSolution(vvi& K, Vector& v, double& d);

private:
    bool CC(const iid& a, const iid& b) const;

    vvd A;
    viid E, M;
};

#endif
