/*
    Copyright (C) 2018 Luka Borozan

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include "BnB.h"
#include "Timer.h"
#include <iostream>
#include <algorithm>

// ./trajan inputs/a915 inputs/a1915 align 1 s 1 0 0.001 2 2>/dev/null

// ./trajan data/T_cslogs_4.tree data/T_cslogs_4.map data/T_cslogs_5.tree data/T_cslogs_5.map output_Sol 2 e 1 0 0.001 2 > t4_t5_log.txt

// ./trajan data/T_cslogs_4.tree data/T_cslogs_4.map data/T_cslogs_5.tree data/T_cslogs_5.map output_Sol 2 e 1 0 0 1 > t4_t5.log.txt

#define DEBUG 1

// struct compare_pair{
//     bool operator()(pair<size_t, size_t>& p, pair<size_t, size_t>& q)
//     { 
//         int i = K[p.first][p.second];
//         int j = K[q.first][q.second];
//         return c(i) > c(j); 
//     }
// };

BnB::BnB(Graph& t1, Graph& t2, string d, double k, bool dag, double c) : LP(t1, t2, d, k, dag), G(Greedy(t1, t2, d, k, dag)), con_eps(c)
{
}

#if DEBUG == 1
int    geno_calls = 0;
double geno_time  = 0;
int    gap_count  = 0;
double gap_amount = 0;
bool   gap_flag   = 0;
#endif

void BnB::Solve(string filename, string outScoreFile)
{
#if DEBUG == 1
    Timer T;
    T.start();
#endif

    G.Solve(""); 
    sys_lb = -G.GetSolution();

    x = Vector::Zero(nr_cols);
    sys_lo.conservativeResizeLike(Vector::Zero(nr_cols));
    sys_hi.conservativeResizeLike(Vector::Ones(nr_cols));
    sys_x.resize(nr_cols);    
    sys_x.assign(nr_cols, false);

    if (SolveLP())    
    {
        x = sys_s;
        for (size_t i = 0; i < x.size(); ++i)
            x(i) = round(x(i));
        WriteSolution(filename);
    }    
    else
        G.WriteSolution(filename);

#if DEBUG == 1
    T.stop();
    cout << endl;
    cout << "total geno calls: " << geno_calls << endl;
    cout << "total geno time: " << geno_time << endl;
    cout << "total time: " << T.secs() << endl;    
    double f = 0;
    for (size_t i = 0; i < x.size(); ++i)
        f += x(i) * c(i);
    cout << "sol: " << f << endl;    
    cout << gap_amount << ' ' << gap_count << endl;
#endif
}

void BnB::Cleanup(size_t nr_t, size_t nr_r)
{
    Triplets.resize(nr_t);    
    nr_rows = nr_r;
}

bool BnB::SolveLP()
{
    int nr_t = Triplets.size();
    int nr_r = nr_rows;
    double f;

#if DEBUG == 1
    int qwe = 0;
    double qwer = 0.0;
    Timer T;
#endif

    while(1)
    {
        SpMat A(nr_rows, nr_cols);
        A.setFromTriplets(Triplets.begin(), Triplets.end());
        Vector b = Vector::Ones(nr_rows);       
        y = Vector::Zero(nr_rows);
        Vector d = -c;            

        BranchingJRF simpleJRF(A, b, d, x, y);
        simpleJRF.lo = sys_lo;
        simpleJRF.hi = sys_hi;
        AugmentedLagrangian solver(simpleJRF, 15);
        solver.setParameter("verbose", false);
        solver.setParameter("pgtol", 0.275); 
        solver.setParameter("constraintsTol", 1e-4);

#if DEBUG == 1
        T.start();
#endif
        SolverStatus status = solver.solve();    
#if DEBUG == 1        
        if (status == SUBOPTIMAL)
            cout << "SUBOPTIMAL(" << status << ")" << endl; 
        if (status == UNBOUNDED)
            cout << "UNBOUNDED(" << status << ")" << endl; 
        if (status == INFEASIBLE)
            cout << "INFEASIBLE(" << status << ")" << endl; 
        if (status == NUM_ERROR)
            cout << "NUM_ERROR(" << status << ")" << endl; 
        T.stop();
        qwe++;
        qwer += T.secs();
        geno_calls++;
        geno_time += T.secs();
#endif

        if (status == INFEASIBLE) 
        {                
            Cleanup(nr_t, nr_r);            
            return 0; 
        }

        x = Vector::ConstMapType(solver.x(), nr_cols); 
        
        if (LP::cf == 1 && Add<1>())
            continue;
        else if (LP::cf == 2 && (Add<1>() + Add<2>()))
            continue;

        f = solver.f();    

#if DEBUG == 1
        cout << "rows: " << nr_rows << endl;
        cout << "cols: " << nr_cols << endl;
        cout << "geno calls: " << qwe << endl;
        cout << "geno time: " << qwer << endl;
        cout << "total geno calls: " << geno_calls << endl;
        cout << "total geno time: " << geno_time << endl;
        cout << "upper: " << f << endl;
        cout << "lower: " << sys_lb << endl;
        int qwert = 0;
        for (size_t i = 0; i < x.size(); ++i)
            if (x(i) > 1e-3 && x(i) < 1 - 1e-3 && !sys_x[i])
                qwert++;
        cout << "fracs: " << qwert << endl;
#endif

        if (sys_lb != double(INF) && f >= sys_lb * (1.0 + con_eps)) 
        {    
#if DEBUG == 1            
            cout << "The branch was cut." << endl << endl;
#endif
            Cleanup(nr_t, nr_r);            
            return 0; 
        }            
#if DEBUG == 1
        else cout << endl;        
#endif

        break;
    }
    
    vector<pair<size_t, size_t>> vp;
    size_t p = 0;

    for (size_t i = 0; i < t1.GetNumNodes(); ++i)
        for (size_t j = 0; j < t2.GetNumNodes(); ++j)
            if (K[i][j] != -1 && x(K[i][j]) > 1e-3 && x(K[i][j]) < 1 - 1e-3 && !sys_x[K[i][j]])                
                vp.emplace_back(i, j);                

    sort(vp.begin(), vp.end(), [this](pair<size_t, size_t>& p, pair<size_t, size_t>& q)
    { 
        int i = K[p.first][p.second];
        int j = K[q.first][q.second];
        return c(i) > c(j); 
    });
   
    for (int i = 0; !p && i < vp.size(); ++i)
    {
        bool flag = 1;
        for (size_t k = 0; flag && k < t1.GetNumNodes(); ++k)          
            for (size_t j = 0; flag && j < t2.GetNumNodes(); ++j)
                if (K[k][j] != -1 && sys_x[K[k][j]])                              
                    flag = IsNotInConflict(vp[i].first, k, vp[i].second, j);                    
        if (flag) p = 1 + K[vp[i].first][vp[i].second];
    }

    if (vp.size())
    {
        bool b0 = !p; 
        p = (b0) ? K[vp[0].first][vp[0].second] : p - 1; // bool b1 = b0 || x(--p) >= 0.1;     
        bool s0 = !b0 && SolveRec(p, 1);            
        bool s1 = SolveRec(p, 0);                       
        if (!s0 && !s1)
        {                 
            Cleanup(nr_t, nr_r);            
            return 0;
        }
    }
    else if (sys_lb > f)
    {
#if DEBUG == 1
        if (gap_flag) 
        {
            gap_amount += abs(sys_lb - f);  
            ++gap_count;
        }
        else 
            gap_flag = 1;
#endif
        sys_s = x;
        sys_lb = f; 
    }

    Cleanup(nr_t, nr_r);
    return 1;
}

bool BnB::SolveRec(size_t p, bool b)
{ 
    sys_lo(p) = b;
    sys_hi(p) = b;
    sys_x[p] = 1;
    
    bool f = SolveLP();    
    
    sys_lo(p) = 0;
    sys_hi(p) = 1;
    sys_x[p] = 0;

    return f;
}
