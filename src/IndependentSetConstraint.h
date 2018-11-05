/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#ifndef INDEPENDENT_SET_CONSTRAINT_H
#define INDEPENDENT_SET_CONSTRAINT_H

#include "Constraint.h"

typedef list<newick_node*> LN;
typedef pair<double, LN> dLN;

class IndependentSetConstraint : Constraint
{
public:
    IndependentSetConstraint(vector<ET>& Triplets, Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp);

    int AddTriplets(int nr_rows);
private:
    void DFSRight(newick_node* noder);
    void DFSLeft(newick_node* nodel, newick_node* noder, double w);
    dLN DFSRight(newick_node* nodel, newick_node* noder);

    vvd D;
};

#endif
