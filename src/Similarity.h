/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#ifndef SIMILARITY_H
#define SIMILARITY_H

#include <list>
#include <string>
using namespace std;

typedef list<string> ls;

extern double var_eps;

double JaccardSim(const ls& L1, const ls& L2, double k);
double SymdifSim(const ls& L1, const ls& L2);

#endif
