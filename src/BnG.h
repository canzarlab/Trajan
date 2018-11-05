/*
    Copyright (C) 2018 Luka Borozan

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/

#ifndef BNG_H
#define BNG_H

#define DEBUG 1

#include "LP.h"
#include "Greedy.h"

#include <vector>
#include <stack>
#include <iostream>
#include <fstream>
#include <cmath>

#include <thread>
#include <mutex>

/*
  Class BnBNode

  It holds all the information required to evaluate LP at a single BnB node.
*/
class BnBNode
{
public:

    // This constructor creates a node without a parent (eg. for creating BnB root node). 
    BnBNode(vector<ET>& Triplets, size_t rows, size_t cols); 

    // This constructor creates a node with a parent (last argument) and warm starts from it.
    BnBNode(BnBNode* node);

    ~BnBNode();

    // Returns whether has BnB fixed a variable at index.
    bool       IsVarFixed(size_t index);

    // Fixes a variable at index to value val.
    void       FixVar(size_t index, double val);

    vector<ET> Triplets; // Stores the system matrix.
    size_t     rows;     // Stores the number of rows of system matrix.
    size_t     cols;         // Stores the number of columns of system matrix.

    Vector     var_lb;   // Stores lower bounds for the variables.
    Vector     var_ub;     // Stores upper bounds for the variables.
    
    Vector*    warm;     // Stores the warm start vector x.
    Vector     sol;      // Stores the solution vector x after the node has been evaluated.
    double     obj;      // Stores the objective value after the node has been evaluated.

#if DEBUG == 1    
    size_t     debug_depth;
    size_t     debug_nodeid;
    size_t     debug_parent;    
    size_t     debug_varid;
    double     debug_varval;
#endif
};

/*
  Class GenericBnBSolver

  This is the template for the BnB solver which takes care of the basic functionallity.

  The idea.
  1. Create Open (set of BnBNode pointers)
  2. Open.Push(BnB_Root_Node)
  3. WHILE (NOT Open.Empty())
  4.     Let (Node \in Open) be a solved node
  5.  IF there are fractionals in Node, add its children to Open

  Open is an arbitrary data structure for storing nodes waiting to be evaluated.
  In steps 4 and 5, different branching strateges may be applied.     
*/
class GenericBnBSolver : public LP
{
public:

    // Same constructor as for the LP solver.
    GenericBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag);

    // Solves the problem and writes the solution to the file. If no file is required pass an empty string.
    void Solve(string filename) override; 

    Vector GetSolution()  { return sys_sol; } // Returns the current solution (copy).
    double GetObjective() { return sys_ub; }  // Returns the current objective function value.

    // Threading interface.
    int  GetNoThreads()      { return thr_no; }
    void SetNoThreads(int n) { thr_no = n;    }

protected:

    // Creates a BnBNode using current solver data and a designated parent node. 
    // If no parent node is needed, pass a nullptr.
    BnBNode*                  InitNodeFrom   (BnBNode* node);

    // Solves the LP in a node using specified pgtol and numtol.
    //    pgtol  - geno solver parameter
    //    numtol - numerical offset used when cutting branches (cut if obj >= best_obj * (1 + numtol)) 
    virtual bool              SolveNode      (BnBNode* node, double pgtol, double numtol);    

    // Push a node into the Open set.
    virtual void              PushNode       (BnBNode* node)              { }

    // Is the Open set empty? If it is, Solver will terminate.    
    virtual bool              OpenEmpty      ()                           { return true; }    

    // Evaluate the Open set and decide which node to branch.    
    virtual vector<BnBNode*>  EvalOpen       ()                           { return vector<BnBNode*>(); }

    // Branch the node and return a vector of children nodes.
    // Return nullptr when there is nothing to branch into.
    virtual vector<BnBNode*>* EvalBranch     (BnBNode* node);
    
    // Fractionallity of the variable. For the 'most fractional' approach,
    // please replace 'c(i)' with '0.5 - abs(0.5 - node->sol(i))'.
    virtual double            VarScore       (int i, BnBNode* node)       { return c(i); }

    // Initialize a node from a parent node and fix variable at index to val.
    virtual BnBNode*          MakeNode       (BnBNode* parent, size_t index, double val); 

    // Auxilliary function which decides whether a variable is to be considered fractional.
    virtual bool              IsVarFrac      (double val)                 { return val > 0.001 && val < 0.999; }    

    // Checks whether the pruning is in order.
    virtual bool                            CheckUB        (double val, double numtol)  { return val >= sys_ub * (1.0 + numtol); }    

    // Decides whether the variable should be preferred to be fixed into 1 or 0.
    virtual bool              CheckFix       (double val)                 { return val >= 0.5; }

    // Checks whether a feasible solution has been obtained.
    virtual bool              CheckSol       ()                           { return sys_sol.size() == nr_rows; }

    // Event callbacks
    virtual void              OnSolverInit   ()                                     { }
    virtual void              OnSolverUpdate ()                                     { }
    virtual bool                            OnSolverFinish ()                                                            { return true; }
    virtual void              OnNodeInit     (BnBNode* node, int index, double val) { }
    virtual bool              OnNodeStart    (BnBNode* node)                                { return true; }
    virtual bool                            OnNodeLP       (BnBNode* node)                        { return true; }
    virtual bool                            OnNodeFinish   (BnBNode* node, bool solved)           { return true; }           

    double sys_ub;   // Stores the best current upper bound.
    Vector sys_sol;  // Stores the best current solution.
    string filename; // Solution file name.

    size_t thr_no;   // Max number of active threads.

private:

    // Pushes all nodes to the Open set.
    void PushAll(vector<BnBNode*>* Nodes);

    bool finished;

#if DEBUG == 1
public: 
    
    mutex       debug_lock;
    ofstream debug_log;
    string   debug_file;
    size_t   debug_nodecnt;
    size_t   debug_genocnt;
    double   debug_genotime;
    double   debug_genomin;
    double   debug_genomax;
#endif
};

/*
  class TestBnBSolver
    
  Example on how to create BnB solvers using the generic solver. It is equivallent to Trajan option 2. 

  class TestBnBSolver : public GenericBnBSolver
  {
  public:

  TestBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag);

  protected:

  // Override when neccessary.
  void              PushNode   (BnBNode* node) override;
  bool              OpenEmpty  ()              override;
  vector<BnBNode*>  EvalOpen   ()              override;

  stack<BnBNode*> Open; // Open set defined as a stack.
  };
*/

/*
  class BFBnBSolver
    
  Best first BnB solver.
*/
class BFBnBSolver : public GenericBnBSolver
{
public:

    BFBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag);

protected:

    void ThreadCallback(vector<BnBNode*>* thr_vec);

    void              PushNode   (BnBNode* node) override;
    bool              OpenEmpty  ()              override;
    vector<BnBNode*>  EvalOpen   ()              override;

    vector<BnBNode*> Open; 
};

/*
  class DFBnBSolver
    
  Depth first BnB solver.
*/
class DFBnBSolver : public GenericBnBSolver
{
public:

    DFBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag);

protected:

    void              PushNode   (BnBNode* node) override;
    bool              OpenEmpty  ()              override;
    vector<BnBNode*>  EvalOpen   ()              override;

    vector<BnBNode*> Open; 
};

/*
  class HybridBnBSolver
    
  Hybrid BnB solver combines the DF and BF approach. It runs DF until the first
  integral solution is hit. After that, BF takes over.
*/
class HybridBnBSolver : public GenericBnBSolver
{
public:

    HybridBnBSolver(Graph& t1, Graph& t2, string dist, double k, bool dag);

protected:

    void              PushNode   (BnBNode* node) override;
    bool              OpenEmpty  ()              override;
    vector<BnBNode*>  EvalOpen   ()              override;

    vector<BnBNode*> Open; 

private:

    void ThreadCallback(vector<BnBNode*>* thr_vec);
};

#endif
