/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include "IndependentSetConstraint.h"

IndependentSetConstraint::IndependentSetConstraint(vector<ET>& Triplets, Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp) : Constraint(Triplets, t1, t2, K, x, swp), D(t1.GetNumNodes(), vd(t2.GetNumNodes()))
{
}

int IndependentSetConstraint::AddTriplets(int nr_rows)
{
    DFSRight(t2.GetRoot());
    int ncr = 0;
    for (newick_node* node : t1.L)
    {
        dLN L = DFSRight(node, t2.GetRoot());
        if (L.first - EPS <= 1)
            continue;

        vii P;
        for (newick_node* noder : L.second)
            for (newick_node* nodel = node; nodel; nodel = nodel->parent ? nodel->parent->node : nullptr)
                P.emplace_back(nodel->taxoni, noder->taxoni);

        AddConstraint(nr_rows + ncr++, P);
    }
    return ncr;
}

void IndependentSetConstraint::DFSRight(newick_node* noder)
{
    DFSLeft(t1.GetRoot(), noder, 0);
    for (newick_child* child = noder->child; child; child = child->next)
        DFSRight(child->node);
}

void IndependentSetConstraint::DFSLeft(newick_node* nodel, newick_node* noder, double w)
{
    if (!nodel->child)
        D[nodel->taxoni][noder->taxoni] = w + GetWeight(nodel, noder);
    for (newick_child* child = nodel->child; child; child = child->next)
        DFSLeft(child->node, noder, w + GetWeight(nodel, noder));
}

dLN IndependentSetConstraint::DFSRight(newick_node* nodel, newick_node* noder)
{
    double w = D[nodel->taxoni][noder->taxoni], sum = 0;
    LN V;

    for (newick_child* child = noder->child; child; child = child->next)
    {
        double ww;
        LN T;
        tie(ww, T) = DFSRight(nodel, child->node);
        sum += ww;
        V.splice(V.begin(), T);
    }

    if (sum > w)
        return make_pair(sum, V);
    return make_pair(w, LN(1, noder));
}
