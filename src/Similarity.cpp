/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include "Similarity.h"
#include <algorithm>
#include <cmath>
#include <vector>

double var_eps;

double JaccardSim(const ls& L1, const ls& L2, double k)
{
    vector<string> I(min(L1.size(), L2.size()));
    auto iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
    I.resize(iit - I.begin());
    double i = I.size(), u = L1.size() + L2.size() - I.size();
    double j = pow(i / u, k);
    return j < var_eps ? 0 : 2 * j;
}

double SymdifSim(const ls& L1, const ls& L2)
{
    vector<string> I(min(L1.size(), L2.size()));
    auto iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
    I.resize(iit - I.begin());
    double j = I.size();
    return j < var_eps ? 0 : 2 * j;
}
