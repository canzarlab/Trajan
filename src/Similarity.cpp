/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Hali.

    Hali is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include "Similarity.h"
#include <algorithm>
#include <cmath>
#include <iostream>
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

double EditDistance(const string & L1, const string & L2, std::vector<std::vector<double>> & cost_matrix)
{
    int i = stoi(L1);
    int j = stoi(L2);
    int n = (int) cost_matrix.size();
    int m = (int) cost_matrix[0].size();
    // cout << (cost_matrix[i][m-1] + cost_matrix[n-1][j] - cost_matrix[i][j]) << " ";
    assert(i < cost_matrix.size() && j < cost_matrix[0].size());
    return (cost_matrix[i][m-1] + cost_matrix[n-1][j] - cost_matrix[i][j]) < var_eps ? 0 : (cost_matrix[i][m-1] + cost_matrix[n-1][j] - cost_matrix[i][j]);
}
