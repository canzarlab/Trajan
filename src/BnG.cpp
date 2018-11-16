/*
    Copyright (C) 2018 Luka Borozan

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

        Example run:
        ./trajan data/1.tree data/2.tree outs/align 2 j 1 0.25 0 7
    84.7252
*/

#include "BnG.h"
#include "Timer.h"

#include <iomanip>

#define PGTOL 0.01 
#define NUMTOL 0.00125
#define THREADS 4

BnBNode::BnBNode(vector<ET>& Triplets, size_t rows, size_t cols) : Triplets(Triplets), rows(rows), cols(cols), warm(nullptr), obj(1)
{
    var_lb.conservativeResizeLike(Vector::Zero(cols));
    var_ub.conservativeResizeLike(Vector::Ones(cols));
#if DEBUG == 1    
    debug_depth = 0;
#endif
}

BnBNode::BnBNode(BnBNode* node) : Triplets(node->Triplets), rows(node->rows), cols(node->cols), obj(1)
{
    warm = new Vector(node->sol);
    var_lb = node->var_lb;
    var_ub = node->var_ub;
#if DEBUG == 1    
    debug_depth = 1 + node->debug_depth;
#endif
}

BnBNode::~BnBNode()
{
    if (warm) delete warm;
}

bool BnBNode::IsVarFixed(size_t index)
{
    return var_lb(index) == var_ub(index);
}

void BnBNode::FixVar(size_t index, double val)
{
    var_lb(index) = var_ub(index) = val;
}

GenericBnBSolver::GenericBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag) : 
    LP(t1, t2, dist, k, dag), sys_ub(666), finished(0), thr_no(THREADS)
{    
}

void GenericBnBSolver::Solve(string filename)
{
    this->filename = filename;
    
#if DEBUG == 1    
    if (debug_file == "")
        debug_log = ofstream(filename + ".log");
    else
        debug_log = ofstream(debug_file + ".log");
    debug_nodecnt  = 0;
    debug_genocnt  = 0;
    debug_genotime = 0;
    debug_genomin  = 1e20;
    debug_genomax  = 0;
#endif

#if DEBUG == 1    
    Timer T; T.start();
#endif

    sys_sol.conservativeResizeLike(Vector::Zero(nr_cols));
    OnSolverInit();

    PushNode(InitNodeFrom(nullptr));
    while (!OpenEmpty())
    {
        for (auto& it : EvalOpen())
        {
            vector<BnBNode*>* V = EvalBranch(it);
#if DEBUG == 1
            if (V != nullptr)
                debug_log << "fixed: " << V->front()->debug_varid << " (" << V->front()->debug_varval << " <- " << it->sol(V->front()->debug_varid) << ") in node " << it->debug_nodeid << endl << endl;
#endif

            if (V != nullptr)
                PushAll(V);
            else if (sys_ub > it->obj)
            {
                sys_sol = it->sol;
                sys_ub  = it->obj;         
#if DEBUG == 1    
                debug_log << "BOUND " << sys_ub << " (" << it->debug_nodeid << ")" << endl << endl;
#endif            
                OnSolverUpdate();
            }
            delete it;
            
        }         
    }
#if DEBUG == 1    
    T.stop();
#endif

    if (OnSolverFinish())
    {
        x = sys_sol;
        for (size_t i = 0; i < x.size(); ++i)
            x(i) = round(x(i));    
        if (filename != "") 
            WriteSolution(filename);     
    }

#if DEBUG == 1    
    debug_log << "total time: " << T.secs() << endl;
    debug_log << "total nodes: " << debug_nodecnt << endl;
    debug_log << "total geno calls: " << debug_genocnt << endl;
    debug_log << "total geno time: " << debug_genotime << " (" << debug_genomin << ", " << debug_genotime / debug_genocnt << "," << debug_genomax << ")" << endl;
    debug_log << "best: " << -sys_ub << endl << endl;
#endif
}

BnBNode* GenericBnBSolver::InitNodeFrom(BnBNode* node)
{
    BnBNode* newnode;
    if (node == nullptr)
    {
        newnode = new BnBNode(Triplets, nr_rows, nr_cols);
        newnode->warm = new Vector(sys_sol);     
    }    
    else
        newnode = new BnBNode(node);
#if DEBUG == 1
    debug_lock.lock();
    newnode->debug_nodeid = debug_nodecnt++;
    newnode->debug_parent = (node == nullptr) ? 0 : node->debug_nodeid;    
    debug_lock.unlock();
#endif
    return newnode;
}

bool GenericBnBSolver::SolveNode(BnBNode* node, double pgtol, double numtol)
{
    Vector x = (node->warm) ? *(node->warm) : Vector::Zero(nr_cols);

#if DEBUG == 1    
    size_t debug_node_genocnt  = 0;
    double debug_node_genotime = 0;
    debug_log << "id: " << node->debug_nodeid << endl;
    debug_log << "parent: " << node->debug_parent << endl;    
    debug_log << "depth: " << node->debug_depth << endl;
    debug_log << "best: " << -sys_ub << endl;
#endif

    if (!OnNodeStart(node)) return false;

    while(1)
    {
        SpMat A(node->rows, node->cols);
        A.setFromTriplets(node->Triplets.begin(), node->Triplets.end());
        Vector b = Vector::Ones(node->rows);       
        Vector y = Vector::Zero(node->rows);        
        Vector d = -c;                    

        BranchingJRF simpleJRF(A, b, d, x, y);
        simpleJRF.lo = node->var_lb;
        simpleJRF.hi = node->var_ub;
        AugmentedLagrangian solver(simpleJRF, 15);
        solver.setParameter("verbose", false);
        solver.setParameter("pgtol", pgtol); 
        solver.setParameter("constraintsTol", 1e-4);

        if (!OnNodeLP(node))
        { 
#if DEBUG == 1    
            debug_log << "ABORT" << endl << endl;
#endif
            return false;
        }
#if DEBUG == 1    
        Timer debug_T; debug_T.start();        
#endif      
        SolverStatus status = solver.solve();    
#if DEBUG == 1    
        debug_T.stop();
        debug_genocnt++;    
        debug_node_genocnt++;
        debug_node_genotime += debug_T.secs();
        debug_genotime += debug_T.secs();
        debug_genomin = min(debug_genomin, (double)debug_T.secs());
        debug_genomax = max(debug_genomax, (double)debug_T.secs());
        if (status == SUBOPTIMAL)
            debug_log << "SUBOPTIMAL(" << status << ")" << endl; 
        if (status == UNBOUNDED)
            debug_log << "UNBOUNDED(" << status << ")" << endl; 
        if (status == INFEASIBLE)
            debug_log << "INFEASIBLE(" << status << ")" << endl; 
        if (status == NUM_ERROR)
            debug_log << "NUM_ERROR(" << status << ")" << endl; 
#endif
        
        if (status == INFEASIBLE || CheckUB(solver.f(), numtol))
        {
#if DEBUG == 1    
            debug_log << "geno calls: " << debug_node_genocnt << endl;
            debug_log << "geno time: " << debug_node_genotime << endl;
            debug_log << "CUT: " << -solver.f() << endl << endl;
#endif
            return false;
        } 

        x = Vector::ConstMapType(solver.x(), node->cols);     

        if ((LP::cf == 1 && Add<1>(node->Triplets, x, node->rows)) || (LP::cf == 2 && (Add<1>(node->Triplets, x, node->rows) + Add<2>(node->Triplets, x, node->rows))))
            continue;
        
        node->obj = solver.f();        
        node->sol = x;

        if (!OnNodeFinish(node, true))
        {
#if DEBUG == 1    
            debug_log << "geno calls: " << debug_node_genocnt << endl;
            debug_log << "geno time: " << debug_node_genotime << endl;
            debug_log << "GREEDY CUT: " << -solver.f() << endl << endl;
#endif
            return false;
        }

        break;
    }

#if DEBUG == 1    
    debug_log << "geno calls: " << debug_node_genocnt << endl;
    debug_log << "geno time: " << debug_node_genotime << endl;
    debug_log << node->rows << " x " << node->cols << endl;
    debug_log << "lp: " << -node->obj << endl;
    size_t debug_fracs = 0;
    for (size_t i = 0; i < x.size(); ++i)
        debug_fracs += IsVarFrac(node->sol(i)) && !node->IsVarFixed(i);
    debug_log << "fracs: " << debug_fracs << endl << endl;
#endif

    return true;
}

vector<BnBNode*>* GenericBnBSolver::EvalBranch(BnBNode* node)
{
    vector<pair<size_t, size_t>> vp;
    size_t p = 0;

    for (size_t i = 0; i < t1.GetNumNodes(); ++i)
        for (size_t j = 0; j < t2.GetNumNodes(); ++j)
            if (K[i][j] != -1 && IsVarFrac(node->sol(K[i][j])) && !node->IsVarFixed(K[i][j]))                
                vp.emplace_back(i, j);                

    sort(vp.begin(), vp.end(), [this, node](pair<size_t, size_t>& p, pair<size_t, size_t>& q)
    {
        return VarScore(K[p.first][p.second], node) > VarScore(K[q.first][q.second], node);  
    });

    for (int i = 0; !p && i < vp.size(); ++i)
    {
        bool flag = 1;
        for (size_t k = 0; flag && k < t1.GetNumNodes(); ++k)          
            for (size_t j = 0; flag && j < t2.GetNumNodes(); ++j)
                if (K[k][j] != -1 && node->IsVarFixed(K[k][j]) && node->var_ub(K[k][j]) != 0)                              
                    flag = IsNotInConflict(vp[i].first, k, vp[i].second, j);                    
        if (flag) p = 1 + K[vp[i].first][vp[i].second];
    }

    if (vp.size())
    {
        vector<BnBNode*>* V = new vector<BnBNode*>;
        bool b0 = !p; 
        p = (b0) ? K[vp[0].first][vp[0].second] : p - 1; 
        if (CheckFix(node->sol(p)))
        {
            if (!b0) V->push_back(MakeNode(node, p, 1));
            V->push_back(MakeNode(node, p, 0));        
        }
        else
        {
            V->push_back(MakeNode(node, p, 0));        
            if (!b0) V->push_back(MakeNode(node, p, 1));
        }
        return V;
    }

    return nullptr;
}

BnBNode* GenericBnBSolver::MakeNode(BnBNode* parent, size_t index, double val)
{
    BnBNode* node = InitNodeFrom(parent);
    node->FixVar(index, val);
#if DEBUG == 1
    node->debug_varid  = index;
    node->debug_varval = val;  
#endif
    OnNodeInit(node, index, val);
    return node;
}

void GenericBnBSolver::PushAll(vector<BnBNode*>* Nodes)
{ 
    for (auto& it : *Nodes) PushNode(it);
    delete Nodes;
}

// ================= TMP\TEST BNB SOLVER ==================== // 
/*
  TestBnBSolver::TestBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag) : GenericBnBSolver(t1, t2, dist, k, dag)
  {    
  }

  void TestBnBSolver::PushNode(BnBNode* node)
  {
  Open.push(node);
  }

  bool TestBnBSolver::OpenEmpty()
  {
  return Open.empty();
  }

  BnBNode* TestBnBSolver::EvalOpen()
  {
  BnBNode* node = Open.top(); Open.pop();
  return SolveNode(node, PGTOL, NUMTOL) ? node : nullptr;
  }
*/
// ================= BF BNB SOLVER ==================== // 

BFBnBSolver::BFBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag) : GenericBnBSolver(t1, t2, dist, k, dag)
{    
}

void BFBnBSolver::PushNode(BnBNode* node)
{
    Open.push_back(node);
}

bool BFBnBSolver::OpenEmpty()
{
    return !Open.size();
}

void BFBnBSolver::ThreadCallback(vector<BnBNode*>* thr_vec)
{
    for (auto& it : *thr_vec)
        if (!SolveNode(it, PGTOL, NUMTOL))
            it->obj = 2;
}

vector<BnBNode*> BFBnBSolver::EvalOpen()
{
    /*
      1. Solve nodes in open with up to thr_no threads.
      2. When all threads are finished pick the best node and return.

      Old code for reference:

      for (auto& it : Open)
      if (it->obj == 1 && !SolveNode(it, PGTOL, NUMTOL))
      it->obj = 2;

      sort(Open.begin(), Open.end(), [](BnBNode*& l, BnBNode*& r)
      { 
      return l->obj < r->obj; //return l->obj > r->obj; 
      });

      BnBNode* node = Open.back(); Open.pop_back();        
      return (node->obj < 1) ? node : nullptr;
    */
    thread thr_obj[thr_no];    
    vector<BnBNode*> thr_vec[thr_no];
    size_t k = 0;    

    for (auto& it : Open)
        if (it->obj == 1) 
            thr_vec[k % thr_no].push_back(it), ++k;
    k = min(k, thr_no);

    for (size_t i = 0; i < k; ++i)
        thr_obj[i] = thread(&BFBnBSolver::ThreadCallback, this, &thr_vec[i]); 

    for (size_t i = 0; i < k; ++i)
        thr_obj[i].join();

    sort(Open.begin(), Open.end(), [](BnBNode*& l, BnBNode*& r)
    { 
        return l->obj < r->obj; 
    });

    while (Open.size() && Open.back()->obj > 0)
        Open.pop_back();

    vector<BnBNode*> v; 
    for (size_t i = 0; Open.size() && i < thr_no; ++i)
        v.push_back(Open.back()), Open.pop_back();

    return v;
}

// ================= DF BNB SOLVER ==================== //

DFBnBSolver::DFBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag) : GenericBnBSolver(t1, t2, dist, k, dag)
{    
}

void DFBnBSolver::PushNode(BnBNode* node)
{
    Open.push_back(node);
}

bool DFBnBSolver::OpenEmpty()
{
    return !Open.size();
}

vector<BnBNode*> DFBnBSolver::EvalOpen()
{
    BnBNode* lt = Open.back(); Open.pop_back();

    if (Open.size() > 1)
    {
        BnBNode* rt = Open.back(); 
        if (lt->obj == 1 && !SolveNode(lt, PGTOL, NUMTOL)) lt->obj = 2;
        if (lt->warm == rt->warm) // TODO this is horrible.
        {
            Open.pop_back();
            if (rt->obj == 1 && !SolveNode(rt, PGTOL, NUMTOL)) rt->obj = 2;
            if (lt->obj < rt->obj) swap(lt, rt); // if (lt->obj > rt->obj) swap(lt, rt);
            if (rt->obj < 1) Open.push_back(rt);         
        }        
    }

    vector<BnBNode*> v;
    if (!(lt->obj == 2 || (lt->obj == 1 && !SolveNode(lt, PGTOL, NUMTOL))))
        v.push_back(lt);
    return v;
}

// ================= Hybrid BNB SOLVER ==================== //

HybridBnBSolver::HybridBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag) : GenericBnBSolver(t1, t2, dist, k, dag)
{    
}

void HybridBnBSolver::PushNode(BnBNode* node)
{
    Open.push_back(node);
}

bool HybridBnBSolver::OpenEmpty()
{
    return !Open.size();
}

void HybridBnBSolver::ThreadCallback(vector<BnBNode*>* thr_vec)
{
    for (auto& it : *thr_vec)
        if (!SolveNode(it, PGTOL, NUMTOL))
            it->obj = 2;
}

vector<BnBNode*> HybridBnBSolver::EvalOpen()
{
    if (!GetObjective())
    {
        BnBNode* lt = Open.back(); Open.pop_back();

        if (Open.size() > 1)
        {
            BnBNode* rt = Open.back(); 
            if (lt->obj == 1 && !SolveNode(lt, PGTOL, NUMTOL)) lt->obj = 2;
            if (lt->warm == rt->warm) // TODO this is horrible.
            {
                Open.pop_back();
                if (rt->obj == 1 && !SolveNode(rt, PGTOL, NUMTOL)) rt->obj = 2;
                if (lt->obj < rt->obj) swap(lt, rt); // if (lt->obj > rt->obj) swap(lt, rt);
                if (rt->obj < 1) Open.push_back(rt);         
            }        
        }

        vector<BnBNode*> v;
        if (!(lt->obj == 2 || (lt->obj == 1 && !SolveNode(lt, PGTOL, NUMTOL))))
            v.push_back(lt);
        return v;
        //return (lt->obj == 2 || (lt->obj == 1 && !SolveNode(lt, PGTOL, NUMTOL))) ? nullptr : lt;
    }
    else
    {
        thread thr_obj[thr_no];    
        vector<BnBNode*> thr_vec[thr_no];
        size_t k = 0;    

        for (auto& it : Open)
            if (it->obj == 1) 
                thr_vec[k % thr_no].push_back(it), ++k;
        k = min(k, thr_no);

        for (size_t i = 0; i < k; ++i)
            thr_obj[i] = thread(&HybridBnBSolver::ThreadCallback, this, &thr_vec[i]); 

        for (size_t i = 0; i < k; ++i)
            thr_obj[i].join();

        sort(Open.begin(), Open.end(), [](BnBNode*& l, BnBNode*& r)
        { 
            return l->obj < r->obj; 
        });

        while (Open.size() && Open.back()->obj > 0)
            Open.pop_back();

        vector<BnBNode*> v; 
        for (size_t i = 0; Open.size() && i < thr_no; ++i)
            v.push_back(Open.back()), Open.pop_back();

        return v;
        //return (node->obj < 1) ? node : nullptr;
    }
}

