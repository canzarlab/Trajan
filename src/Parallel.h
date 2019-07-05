/*
    Copyright (C) 2018 Luka Borozan

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version. 
*/

#ifndef Parallel_H
#define Parallel_H

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "BnG.h"

#define MAX_THREADS 15
#define GREEDYTOL 0.1

class ParallelSolver
{

public:
  
    ParallelSolver(Graph& t1, Graph& t2, string d, double k, bool dag, int nthreads);
  
    // Solves the model and writes it down to "filename".
    void Solve(string filename, string outScoreFile = "score.csv");

    // Updates the best upper bound and solution with regards to locking.
    bool PushUB(Vector& var, double  val, GenericBnBSolver& solver);
    bool PullUB(Vector& var, double& val, GenericBnBSolver& solver);

    double  GetUBVal(); 
    Vector& GetUBVar(); 

    bool    Finished(); 

private:

    void Callback(string filename, GenericBnBSolver* solver);
    
    // threading locks and return values
    mutex                             thr_block; // Locks the best upper bound.
    mutex                             thr_slock; // Locks the solution.
    condition_variable thr_cond;  // Finishes all threads once one solver finds an optimal solution.
    atomic<thread::id> thr_val;   // Value of the finished thread.
    int                thr_num;   // Number of running threads.

    // Best upper bound and solution.
    Vector sys_sol;
    double sys_ub;

    // Solver related data. Passing this to the constructors of each solver.
    Graph& t1, t2;
    string d;
    double k;
    bool dag;

    bool sol; // 
};

// Different BnB solvers
#define PAR_CLASS(X, P, VS)                                             \
    class X : public P                                                  \
    {                                                                   \
    public:                                                             \
                                                                        \
        X(Graph& t1, Graph& t2, string dist, double k, bool dag, ParallelSolver& par) : \
            P(t1, t2, dist, k, dag), par(par)                           \
        {                                                               \
        }                                                               \
                                                                        \
    protected:                                                          \
                                                                        \
        virtual double VarScore(int i, BnBNode* node) VS                \
                                                                        \
            bool OnNodeLP(BnBNode* node)                                \
        {                                                               \
            par.PullUB(sys_sol, sys_ub, *this);                         \
            return !par.Finished();                                     \
        }                                                               \
                                                                        \
        void OnSolverUpdate()                                           \
        {                                                               \
            par.PushUB(sys_sol, sys_ub, *this);                         \
        }                                                               \
                                                                        \
        void OnSolverInit()                                             \
        {                                                               \
            Greedy G(t1, t2, d, k, dag); G.Solve("");                   \
            G.GetSolution(K, sys_sol, sys_ub); grd_lb = sys_ub;         \
        }                                                               \
                                                                        \
        bool OnSolverFinish()                                           \
        {                                                               \
            par.PullUB(sys_sol, sys_ub, *this);                         \
            return !par.Finished();                                     \
        }                                                               \
                                                                        \
        bool OnNodeFinish(BnBNode* node, bool flag)                     \
        {                                                               \
            if (!flag) return true;                                     \
            Vector& x = node->sol;                                      \
                                                                        \
            map<size_t, bool> M;                                        \
            for (int i = 0; i < x.size(); ++i)                          \
                if (node->IsVarFixed(i)) M[i] = node->var_ub(i); else if (1 - x(i) < 0.01) M[i] = 1; \
                                                                        \
            Greedy T(t1, t2, d, k, K, M); T.Solve(""); double f = -T.GetSolution(); grd_lb = min(grd_lb, f); \
            if (f < sys_ub * 1.001)                                     \
            {                                                           \
                T.GetSolution(K, sys_sol, sys_ub);                      \
                OnSolverUpdate();                                       \
                return true;                                            \
            }                                                           \
            return f * (1 + GREEDYTOL) < grd_lb;                        \
        }                                                               \
                                                                        \
        void OnNodeInit(BnBNode* node, int index, double val)           \
        {                                                               \
            if (val != 1) return;                                       \
            size_t m, n;                                                \
                                                                        \
            for (size_t i = 0; i < t1.GetNumNodes(); ++i)               \
                for (size_t j = 0; j < t2.GetNumNodes(); ++j)           \
                    if (K[i][j] == index) m = i, n = j;                 \
                                                                        \
            for (size_t i = 0; i < t1.GetNumNodes(); ++i)               \
                for (size_t j = 0; j < t2.GetNumNodes(); ++j)           \
                    if (K[i][j] != -1 && !node->IsVarFixed(K[i][j]) && ((m == i && n != j || m != i && n == j) || !IsNotInConflict(m, i, n, j))) \
                        node->FixVar(K[i][j], 0);                       \
        }                                                               \
                                                                        \
        ParallelSolver& par;                                            \
        double grd_lb;                                                  \
    };

PAR_CLASS(BnBBFMF, BFBnBSolver, { return 0.5 - abs(0.5 - node->sol(i)); })
PAR_CLASS(BnBBFLF, BFBnBSolver, { return abs(0.5 - node->sol(i)); })
PAR_CLASS(BnBBFWF, BFBnBSolver, { return c(i) * node->sol(i); })
PAR_CLASS(BnBBFF,  BFBnBSolver, { return node->sol(i); })
PAR_CLASS(BnBBFA,  BFBnBSolver, ;)

PAR_CLASS(BnBDFMF, DFBnBSolver, { return 0.5 - abs(0.5 - node->sol(i)); })
PAR_CLASS(BnBDFLF, DFBnBSolver, { return abs(0.5 - node->sol(i)); })
PAR_CLASS(BnBDFWF, DFBnBSolver, { return c(i) * node->sol(i); })
PAR_CLASS(BnBDFF,  DFBnBSolver, { return node->sol(i); })
PAR_CLASS(BnBDFA,  DFBnBSolver, ;)

PAR_CLASS(BnBHMF, HybridBnBSolver, { return 0.5 - abs(0.5 - node->sol(i)); })
PAR_CLASS(BnBHLF, HybridBnBSolver, { return abs(0.5 - node->sol(i)); })
PAR_CLASS(BnBHWF, HybridBnBSolver, { return c(i) * node->sol(i); })
PAR_CLASS(BnBHF,  HybridBnBSolver, { return node->sol(i); })
PAR_CLASS(BnBHA,  HybridBnBSolver, ;)

#endif
