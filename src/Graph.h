/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <map>
#include <list>
#include "newick.h"
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<double> vd;
typedef vector<vd> vvd;
typedef vector<vb> vvb;

const size_t NR_THREADS = 4;

class Graph
{
public:
    Graph() : _n(0) { }
    virtual ~Graph() { dealloc_dag(root, _n); };

    int GetNumNodes() const { return _n; }
    newick_node* GetRoot() const { return root; }
    virtual Graph* Init();

    vn GetNodes()
    {
        return get_nodes(root, _n);
    }

    vn L;
    mnls clade;
    vvb D;
protected:
    void TransitiveClosure();
    void TransitiveClosure(newick_node* node, newick_node* rnode, vvb& C);
    void Init(newick_node* node);
    virtual void Leaf(newick_node* node) { }
    virtual void Child(newick_node* node, newick_node* child) { }
    virtual void Relation(int d, int a) { }
    newick_node* root;
    long _n;
};

class DAG : public Graph
{
public:
    DAG(const char* f1, const char* f2);
    ~DAG() { }

protected:
    virtual void Relation(int d, int a) { }
};

class Tree : public Graph
{
public:
    Tree() { }
    Tree(const char* f1);
    Tree(const char* f1, const char* f2);
    ~Tree() { }

protected:
    virtual void Leaf(newick_node* node) override;
    virtual void Child(newick_node* node, newick_node* child) override;
    bool t;
};

class LDAG : public DAG
{
public:
    LDAG(const char* f1, const char* f2);
    Graph* Init();

    vvi G;
    vvd R[NR_THREADS];
    vector<vn> P;
protected:
    void GenPaths(newick_node* node, vn& P);
    void TransitiveReduction(newick_node* node, vb& C);
    void TransitiveReduction(newick_node* parent, newick_node* node, vb& C);
    void Reduce(newick_child** childptr, newick_node* node);
    void Wipe(newick_node* node);
    void Renumerate(newick_node* node);
    template<class F>
    void ForeachPair(vn& Q, F f)
    {
        for (int i = 0; i < Q.size(); ++i)
            for (int j = 1; j < Q.size(); ++j)
                f(Q, Q[i]->taxoni, Q[j]->taxoni);
    }
    virtual void Relation(int d, int a) override;
};

#endif
