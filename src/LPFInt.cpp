/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include "LPFInt.h"

LPFInt::LPFInt(Graph& t1, Graph& t2, string d, double k, bool dag) : LPInt(t1, t2, d, k, dag)
{
}

void LPFInt::Solve(string filename, string outScoreFile)
{
    LP::Solve(filename);
    LPInt::Solve(filename);
}
