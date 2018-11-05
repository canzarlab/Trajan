/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#ifndef LP_H
#define LP_H

#include<iostream>

#include "Solver.h"
#include "Geno.h"
#include "IndependentSetConstraint.h"
#include "AntichainConstraint.h"
#include "CrossingConstraint.h"

extern int c;

class LP : public Solver
{
public:
    LP(Graph& t1, Graph& t2, string d, double k, bool dag);
    ~LP();

    virtual void Solve(string filename);
    void WriteSolution(string fileName);

protected:
    template<class T> int Add();
    template<int N> int Add();
    template<class T> int Add(vector<ET>& Triplets, Vector& x, size_t& nr_rows);
    template<int N> int Add(vector<ET>& Triplets, Vector& x, size_t& nr_rows);
    virtual bool SolveLP();
    void MatchingConstraints();

    vector<ET> Triplets;
    // backup x->warm_x and y->warm_y for two consecutive iterations
    Vector warm_x, warm_y;
    Vector x, y;
    vvi K;
    Vector c;
    int nr_rows, nr_cols;
    int cnt;
};

template<class T>
int LP::Add()
{
    int row_old = nr_rows;
    T c12(Triplets, t1, t2, K, x, false);
    nr_rows += c12.AddTriplets(nr_rows);
    T c21(Triplets, t2, t1, K, x, true);
    nr_rows += c21.AddTriplets(nr_rows);
    return nr_rows - row_old;
}

template<int N>
int LP::Add()
{
    static_assert(N == 1 || N == 2, "invalid constraint");
    if (N == 1)
        return Add<CrossingConstraint>();
    return (dag ? Add<AntichainConstraint>() : Add<IndependentSetConstraint>());
}

template<class T>
int LP::Add(vector<ET>& Triplets, Vector& x, size_t& nr_rows)
{
    int row_old = nr_rows;
    T c12(Triplets, t1, t2, K, x, false);
    nr_rows += c12.AddTriplets(nr_rows);
    T c21(Triplets, t2, t1, K, x, true);
    nr_rows += c21.AddTriplets(nr_rows);
    return nr_rows - row_old;
}

template<int N>
int LP::Add(vector<ET>& Triplets, Vector& x, size_t& nr_rows)
{
    static_assert(N == 1 || N == 2, "invalid constraint");
    if (N == 1)
        return Add<CrossingConstraint>(Triplets, x, nr_rows);
    return (dag ? Add<AntichainConstraint>() : Add<IndependentSetConstraint>(Triplets, x, nr_rows));
}

#endif
