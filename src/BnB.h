/*
    Copyright (C) 2018 Luka Borozan

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#ifndef BNB_H
#define BNB_H

#include "LP.h"
#include "Greedy.h"

class BnB : public LP
{
public:
    BnB(Graph& t1, Graph& t2, string d, double k, bool dag, double c);
    virtual void Solve(string filename, string outScoreFile = "score.csv") override;

private:
    void Cleanup(size_t nr_t, size_t nr_r);
    bool SolveLP() override; // override je bilo   
    bool SolveRec(size_t p, bool b);

    Greedy       G;
    double       sys_lb;
    vector<bool> sys_x;
    Vector       sys_s;
    Vector       sys_lo;
    Vector       sys_hi;
    
    double       con_eps;
};

#endif
