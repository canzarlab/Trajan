/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include "Constraint.h"

Constraint::Constraint(vector<ET>& Triplets, Graph& t1, Graph& t2, vvi& K, Vector& x, bool swp) : Triplets(Triplets), t1(t1), t2(t2), K(K), x(x), swp(swp)
{
}

void Constraint::AddConstraint(int row, vii& P)
{
    for (auto k : P)
    {
        int col = GetCol(k.first, k.second);
        if (col != -1)
            Triplets.emplace_back(row, col, 1.);
    }
}
