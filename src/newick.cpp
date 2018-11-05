/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include "newick.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

newick_child::newick_child(newick_node* node) : node(node), next(nullptr)
{
}

newick_child::~newick_child()
{
    delete next;
}

newick_node::newick_node(const string& taxon, float dist, newick_child* child) : child(child), taxon(taxon), dist(dist), parent(nullptr)
{
    try { taxoni = stoi(taxon); } catch(...) { taxoni = -1; }
}

newick_node::~newick_node()
{
    delete parent;
    delete child;
}

static void get_nodes(newick_node* node, vn& N)
{
    N[node->taxoni] = node;
    for (newick_child* child = node->child; child; child = child->next)
        if (!N[child->node->taxoni])
            get_nodes(child->node, N);
}

vn get_nodes(newick_node* node, int n)
{
    vn N(n);
    get_nodes(node, N);
    return move(N);
}

void dealloc_dag(newick_node* node, int n)
{
    vn N = get_nodes(node, n);
    for (newick_node* node : N)
        delete node;
}

static newick_node* parse_node(string& str, size_t& itr, newick_child* child)
{
    size_t index = str.find_first_of(",);\r\n\0", itr);
    string data = str.substr(itr, index - itr);
    itr = index;

    size_t index2 = data.find_first_of(':');
    string taxon = data.substr(0, index2);
    float dist = 0;
    if (index2 != string::npos && data.back() != ':')
        dist = stof(data.substr(index2 + 1));

    return new newick_node(taxon, dist, child);
}

static newick_node* parse_tree(string& str, size_t& itr)
{
    newick_child* child = nullptr;
    newick_child** childptr = &child;
    for (itr += 1; itr < str.size() && str[itr - 1] != ')'; ++itr)
    {
        newick_node* new_child = nullptr;
        if (str[itr] == '(')
            new_child = parse_tree(str, itr);
        else
            new_child = parse_node(str, itr, nullptr);

        *childptr = new newick_child(new_child);
        childptr = &(*childptr)->next;
    }
    return parse_node(str, itr, child);
}

newick_node* load_tree(const char* filename)
{
    ifstream File(filename);
    if (!File)
    {
        cout << "Cannot read " << filename << endl;
        return nullptr;
    }

    int Size;
    File.seekg(0, ios::end);
    Size = File.tellg();
    File.seekg(0, ios::beg);

    string str(Size, 0);
    File.read(&str[0], Size);

    size_t itr = 0;
    return parse_tree(str, itr);
}

typedef vector<string> vs;
typedef map<string, vs> msvs;
typedef pair<string, vs> svs;
typedef pair<string, string> ss;
typedef vector<ss> vss;

static newick_node* load_dag_internal(const string& r, msn& M, msvs& C)
{
    if (newick_node* node = M[r])
        return node;

    newick_child* child = nullptr;
    newick_child** childptr = &child;
    for (const string& i : C[r])
    {
        *childptr = new newick_child(load_dag_internal(i, M, C));
        childptr = &(*childptr)->next;
    }
    return M[r] = new newick_node(r, 0, child);
}

newick_node* load_dag(const char* f1, const char* f2, mnls& clade, msn& M)
{
    vss S;
    msvs C, P;
    ifstream ef(f1);
    if (!ef)
        return nullptr;

    string n1, n2, s;
    while (ef >> n1 >> n2 >> s)
    {
        if (s == "default")
        {
            if (!f2) swap(n1, n2);
            P[n1].push_back(n2);
            P[n2];
            C[n2].push_back(n1);
        }
        else if (s == "gene" && !f2)
            S.emplace_back(n1, n2);
        else
            cerr << "WARNING: unsupported entry!\n";
    }

    const auto& l = [](const svs& i){return i.second.empty();};
    const string& r = find_if(P.begin(), P.end(), l)->first;
    newick_node* root = load_dag_internal(r, M, C);

    if (f2)
    {
        ifstream lf(f2);
        if (!lf)
            return nullptr;

        while (lf >> n1 >> n2)
            clade[M[n2]].push_back(n1);
    }
    else
        for (const ss& p : S)
            clade[M[p.first]].push_back(p.second);

    return root;
}

void print_tree(newick_node* root, ostream& file)
{
    file << fixed << setprecision(6);
    if (root->child)
    {
        file << "(";
        for (newick_child* child = root->child; child; child = child->next)
        {
            print_tree(child->node, file);
            if (child->next)
                file << ",";
        }
        file << ")" << root->taxon << ":" /*<< root->dist*/;
    }
    else
        file << root->taxon << ":" /*<< root->dist*/;
}
