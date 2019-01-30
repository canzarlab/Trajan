#include <iostream>
#include <fstream>
#include "newick.h"
using namespace std;

void trans_closure(newick_node* node, newick_node* rnode, vector<vector<int>>& D)
{
    D[rnode->taxoni][node->taxoni] = true;
    for (newick_child* child = node->child; child; child = child->next)
        trans_closure(child->node, rnode, D);
}

void make_trans(newick_node* node, vector<vector<int>>& D)
{
    trans_closure(node, node, D);
    for (newick_child* child = node->child; child; child = child->next)
        make_trans(child->node, D);
}

void make_map(newick_node* node, map<string, newick_node*>& m, int& cnt)
{
    m[node->taxon] = node;
    node->taxoni = cnt++;
    for (newick_child* child = node->child; child; child = child->next)
        make_map(child->node, m, cnt);
}

int cc0, cc1, cc2;

bool check_conflict(tuple<int, int> a, tuple<int, int> b, vector<vector<int>>& D1, vector<vector<int>>& D2)
{
    int i = get<0>(a), j = get<0>(b);
    int x = get<1>(a), y = get<1>(b);
    bool c2 = (D1[j][i] || D1[i][j]) == (D2[x][y] || D2[y][x]);
    bool c1 = (D1[i][j] == D2[x][y]); // assuming c2 is satisfied
    bool c0 = (i != j && x != y);
    if (!c0)
        cc0++;
    else if (!c2)
        cc2++;
    else if (!c1)
        cc1++;
    else
        return true;
    return false;
}


int main(int argc, char** argv)
{
    if (argc != 4)
    {
        cout << "usage: " << argv[0] << " <t1.newick> <t2.newick> <align>" << endl;
        return EXIT_FAILURE;
    }
    newick_node* t1 = load_tree(argv[1]);
    newick_node* t2 = load_tree(argv[2]);
    map<string, newick_node*> m1, m2;
    int cnt1 = 0;
    make_map(t1, m1, cnt1);
    int cnt2 = 0;
    make_map(t2, m2, cnt2);
    vector<vector<int>> D1(cnt1, vector<int>(cnt1)), D2(cnt2, vector<int>(cnt2));
    make_trans(t1, D1);
    make_trans(t2, D2);
    ifstream align(argv[3]);
    string aa, bb;
    vector<tuple<int, int> > E;
    while (align >> aa >> bb)
        E.emplace_back(m1[aa]->taxoni, m2[bb]->taxoni);

    for (int i = 0; i < E.size(); ++i)
        for (int j = i + 1; j < E.size(); ++j)
            check_conflict(E[i], E[j], D1, D2);

    cout << "C0 violations: " << cc0 << endl;
    cout << "C1 violations: " << cc1 << endl;
    cout << "C2 violations: " << cc2 << endl;
}
