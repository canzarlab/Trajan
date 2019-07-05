/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Hali.

    Hali is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include "Graph.h"
#include "Geno.h"
#include <algorithm>

DAG::DAG(const char* f1, const char* f2)
{
    msn M;
    Init(root = load_dag(f1, f2, clade, M));
}

LDAG::LDAG(const char* f1, const char* f2) : DAG(f1, f2)
{
    size_t SZ = _n * 2 + 2;
    G.resize(SZ);
    for (int j = 0; j < NR_THREADS; ++j)
        R[j].resize(SZ, vd(SZ));

    int S = SZ - 2, T = SZ - 1;
    for (int i = 0; i < _n; ++i)
    {
        G[S].push_back(i);
        G[i].push_back(S);
        G[T].push_back(i + _n);
        G[i + _n].push_back(T);
    }
}

Graph* Graph::Init()
{
    TransitiveClosure();
    return this;
}

Graph* LDAG::Init()
{
    vn T;
    vb C(_n);
    TransitiveReduction(root, C);
    // Restore dfs ordering of indices
    Wipe(root);
    _n = 0;
    Renumerate(root);
    GenPaths(root, T);
    sort(P.begin(), P.end(), [](const vn& a, const vn& b)
    {
        return a.size() > b.size();
    });
    vvb D(_n, vb(_n));
    vector<vn> Q;
    for (int k = 0; k < P.size(); ++k)
    {
        ForeachPair(P[k], [&](vn& p, int i, int j)
        {
            if (D[i][j]) return;
            Q.push_back(p);
            ForeachPair(P[k], [&](vn& p, int i, int j)
            {
                D[i][j] = D[j][j] = true;
            });
        });
    }
    P = move(Q);
    return Graph::Init();
}

void LDAG::Wipe(newick_node* node)
{
    node->taxoni = -1;
    for (newick_child* child = node->child; child; child = child->next)
        if (child->node->taxoni != -1)
            Wipe(child->node);
}

void LDAG::Renumerate(newick_node* node)
{
    node->taxoni = _n++;
    for (newick_child* child = node->child; child; child = child->next)
        if (child->node->taxoni == -1)
            Renumerate(child->node);
}

void LDAG::TransitiveReduction(newick_node* node, vb& C)
{
    if (C[node->taxoni])
        return;

    C[node->taxoni] = true;
    for (newick_child* child = node->child; child; child = child->next)
    {
        TransitiveReduction(child->node, C);
        for (newick_child* cchild = child->node->child; cchild; cchild = cchild->next)
        {
            vb CC(_n);
            TransitiveReduction(node, cchild->node, CC);
        }
    }
}

void LDAG::TransitiveReduction(newick_node* parent, newick_node* node, vb& C)
{
    C[node->taxoni] = true;
    Reduce(&parent->child, node);
    Reduce(&node->parent, parent);
    for (newick_child* child = node->child; child; child = child->next)
        if (!C[child->node->taxoni])
            TransitiveReduction(parent, child->node, C);
}

void LDAG::Reduce(newick_child** childptr, newick_node* node)
{
    while (*childptr)
    {
        if ((*childptr)->node == node)
            *childptr = (*childptr)->next;
        else childptr = &(*childptr)->next;
    }
}

void LDAG::GenPaths(newick_node* node, vn& T)
{
    T.push_back(node);
    if (!node->child)
        P.push_back(T);
    for (newick_child* child = node->child; child; child = child->next)
        GenPaths(child->node, T);
    T.pop_back();
}

void Graph::TransitiveClosure()
{
    D.resize(_n, vb(_n));
    vvb C(_n, vb(_n));
    for (newick_node* leaf : L)
        TransitiveClosure(leaf, leaf, C);
}

void Graph::TransitiveClosure(newick_node* node, newick_node* rnode, vvb& C)
{
    int l = rnode->taxoni;
    int i = node->taxoni;
    if (l != i)
    {
        D[l][i] = true;
        Relation(l, i);
    }

    C[i][l] = true;
    for (newick_parent* parent = node->parent; parent; parent = parent->next)
    {
        newick_node* pn = parent->node;
        int pnt = pn->taxoni;
        if (!C[pnt][l])
            TransitiveClosure(pn, rnode, C);
        if (!C[pnt][pnt])
            TransitiveClosure(pn, pn, C);
    }
}

/*
#include <algorithm>

for (string& j : L[i])
    if (find(L[r].begin(), L[r].end(), j) == L[r].end())
        L[r].push_back(j);
*/
void Graph::Init(newick_node* node)
{
    if (!node->child)
    {
        L.push_back(node);
        Leaf(node);
    }

    node->taxoni = _n++;
    for (newick_child* child = node->child; child; child = child->next)
    {
        newick_node* cnode = child->node;
        newick_parent** parentptr = &cnode->parent;
        if (!cnode->parent)
            Init(cnode);

        while (*parentptr) parentptr = &(*parentptr)->next;
        *parentptr = new newick_parent(node);
        Child(node, child->node);
    }
    clade[node].sort();
}

Tree::Tree(const char* f1, const char* f2) : t(true)
{
    msn M;
    Init(root = load_dag(f1, f2, clade, M));
    TransitiveClosure();
}

Tree::Tree(const char* f1) : t(false)
{
    Init(root = load_tree(f1));
    TransitiveClosure();
}

void Tree::Leaf(newick_node* node)
{
    if (!t)
        clade[node].push_back(node->taxon);
}

void Tree::Child(newick_node* node, newick_node* child)
{
    ls &cl = clade[node], &cr = clade[child];
//     cl.insert(cl.end(), cr.begin(), cr.end());
}

void LDAG::Relation(int l, int i)
{
    for (int j = 0; j < NR_THREADS; ++j)
        R[j][l][i + _n] = INF;
    G[l].push_back(i + _n);
    G[i + _n].push_back(l);
}
