#include "newick.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <utility>
#include <algorithm>
#include <cassert>
#include <numeric>
using namespace std;

newick_node* yule(int n)
{
    newick_node* root = new newick_node;
    vector<newick_node*> L;
    L.push_back(root);
    while (L.size() < n)
    {
        newick_node* leaf1 = new newick_node;
        newick_node* leaf2 = new newick_node;
        newick_child* child = new newick_child(leaf1);
        child->next = new newick_child(leaf2);
        int i = rand() % L.size();
        L[i]->child = child;
        L[i] = leaf1;
        L.push_back(leaf2);
    }

    unsigned random[n];
    iota(random, random + n, 1);
    random_shuffle(random, random + n);
    for (int i = 0; i < L.size(); ++i)
        L[i]->taxon = to_string(random[i]);

    return root;
}

newick_node* uniform(int n)
{
    unsigned random[n], size = 2;
    iota(random, random + n, 1);
    random_shuffle(random, random + n);

    newick_node* leaf1 = new newick_node(to_string(random[0]));
    newick_node* leaf2 = new newick_node(to_string(random[1]));
    newick_child* child = new newick_child(leaf1);
    child->next = new newick_child(leaf2);

    newick_node* root = new newick_node("", 0, child);
    vector<pair<newick_node*, newick_node*>> E;
    E.emplace_back(root, leaf1);
    E.emplace_back(root, leaf2);

    for (int k = 0; k < n - 2; k++)
    {
        int i = rand() % E.size();
        leaf1 = new newick_node(to_string(random[size++]));
        leaf2 = new newick_node("", 0, new newick_child(E[i].second));
        leaf2->child->next = new newick_child(leaf1);

        newick_child* pchild = E[i].first->child;
        if (pchild->node == E[i].second)
            pchild->node = leaf2;
        else
            pchild->next->node = leaf2;

        E.emplace_back(leaf2, leaf1);
        E.emplace_back(leaf2, E[i].second);
        E[i].second = leaf2;
    }
    return root;
}

newick_node* (*F[])(int) = { uniform, yule };

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        cout << "usage: " << argv[0] << " <n-leaves> <n-trees>\n";
        return 1;
    }
    srand(unsigned(time(0)));
    int n_leaves = stoi(argv[1]);
    int n_trees = stoi(argv[2]);
    for (int i = 0; i < n_trees; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            newick_node* root = F[j](n_leaves);
            ofstream file(string("data") + to_string(j) + string("/a") + to_string(i));
            assert(file);
            print_tree(root, file);
            file << '\n';
            delete root;
        }
    }
}
