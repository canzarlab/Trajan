/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "Geno.h"
#include "Graph.h"
#include <utility>
using namespace std;

const double EPS = 1e-2;

typedef vector<pair<int, int> > vii;

class Constraint
{
public:
    Constraint(vector<ET>& Triplets, Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp);

    virtual int AddTriplets(int nr_rows) = 0;
protected:
    inline int GetCol(int i, int j) const;
    inline int GetCol(newick_node* nodel, newick_node* noder) const;
    inline double GetWeight(newick_node* nodel, newick_node* noder) const;
    inline double GetWeight(int i, int j) const;
    void AddConstraint(int row, vii& P);

    vector<ET>& Triplets;
    Graph &t1, &t2;
    vvi& K;
    Vector& x;
    bool swp;
};

inline int Constraint::GetCol(int i, int j) const
{
    return swp ? K[j][i] : K[i][j];
}

inline int Constraint::GetCol(newick_node* nodel, newick_node* noder) const
{
    return GetCol(nodel->taxoni, noder->taxoni);
}

inline double Constraint::GetWeight(newick_node* nodel, newick_node* noder) const
{
    return GetWeight(nodel->taxoni, noder->taxoni);
}

inline double Constraint::GetWeight(int i, int j) const
{
    int in = GetCol(i, j);
    return in == -1 ? 0 : x(in);
}

#endif
