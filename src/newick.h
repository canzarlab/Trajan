/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#ifndef NEWICK_H
#define NEWICK_H

#include <string>
#include <map>
#include <vector>
#include <list>
using namespace std;

struct newick_node;

struct newick_child
{
    newick_child(newick_node* node);
    ~newick_child();

    newick_node* node;
    newick_child* next;
};

typedef newick_child newick_parent;

struct newick_node
{
    newick_node(const string& taxon = "", float dist = 0, newick_child* child = nullptr);
    ~newick_node();

    newick_child* child;
    string taxon;
    int taxoni;
    float dist;
    newick_parent* parent;
};

typedef map<string, newick_node*> msn;
typedef vector<newick_node*> vn;
typedef list<string> ls;
typedef map<newick_node*, ls> mnls;
typedef vector<bool> vb;

newick_node* load_tree(const char* filename);
newick_node* load_dag(const char* f1, const char* f2, mnls& clade, msn& M);
void dealloc_dag(newick_node* node, int n);
void print_tree(newick_node* root, ostream& file);
vn get_nodes(newick_node* node, int n);

#endif

