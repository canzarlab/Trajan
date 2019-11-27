/*
    Copyright (C) 2019 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>
#include <functional>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <cassert>
#define die() assert(false)
using namespace std;

using mask = unsigned long long;
// TODO: replace with gp_hash_table from __gnu_pbds or some other faster hash table
using dp_table = unordered_map<mask, unordered_map<mask, int>>;

// bitmask operations
inline mask pc(mask x) { return __builtin_popcountll(x); }
inline mask fm(int x) { return (1ull << x) - 1; }
inline mask in(int x) { return 1ull << x; };
inline mask ls(mask x) { return x & -x; };
inline mask tz(mask x) { return __builtin_ctzll(x); }
inline mask em(mask m, mask r) { return m & ~r; }

// maximum number of nodes in a tree
const int MAXN = 500;
// maximum number of leaf/branching nodes in a tree
const int MAXLB = 64;

dp_table fdp;

// should be faster than vector<vector<int>>
struct array2d
{
    int n, m;
    vector<int> v;

    array2d(int n, int m) : n(n), m(m)
    {
        v.resize(n * m, -1);
    }

    int& operator()(int x, int y)
    {
        return v[x * m + y];
    }
};

struct tree
{
    // mapping of node identifiers to nodes
    map<string, int> m;
    // adjacency list of a node
    vector<int> adj[MAXN];
    // number of nodes
    int n;
    // number of leaf/branching nodes
    int nb;
    // root node
    int r;
    // topological ordering of leaf/branching nodes
    int b[MAXLB];
    // unique parent of a node
    int par[MAXN];
    // leaf/branching descendant bitmask of a leaf/branching node
    mask ch[MAXLB];
    // node -> leaf/branching node index (-1 if not leaf/branching)
    int rb[MAXN];
    // branching ancestor mask
    mask ba[MAXLB];

    tree(const char* filename) : n(0), nb(0), r(-1)
    {
        ifstream file(filename);
        if (!file)
            die();

        // map string indices to integers
        auto idx = [&](string& x)
        {
            auto it = m.find(x);
            if (it == m.end())
                return m[x] = n++;
            return it->second;
        };

        // nroot[x] == 1 means x is not the root
        vector<int> nroot(MAXN);

        // read from file
        string sa, sb, type;
        while (file >> sa >> sb >> type)
        {
            if (type != "default")
                die();

            // directed edge sb->sa
            int x = idx(sa);
            int y = idx(sb);
            adj[y].push_back(x);
            nroot[x] = 1;
        }

        // tree too big
        if (m.size() > MAXN)
            die();

        // find root
        for (int i = 0; i < m.size(); ++i)
        {
            // found the root
            if (!nroot[i])
            {
                // multiple roots
                if (r != -1)
                    die();

                r = i;
            }
        }

        // tree has no root
        if (r == -1)
            die();

        // identify branching nodes
        function<void(int)> mark = [&] (int x)
        {
            for (int y : adj[x])
            {
                // fill the parent table
                par[y] = x;
                mark(y);
            }
            // x is a leaf/branching node - record finishing time
            if (adj[x].size() != 1)
                b[nb++] = x;
        };
        // root has no parent
        par[r] = -1;
        mark(r);

        // too many leaf/branching nodes
        if (nb > MAXLB)
            die();

        // obtain topological ordering
        reverse(b, b + nb);

        // create node->leaf/branching node index mapping
        memset(rb, -1, sizeof rb);
        for (int i = 0; i < nb; ++i)
            rb[b[i]] = i;

        // fill the table of immediate leaf/branching descendant masks
        function<mask(int)> child = [&] (int x)
        {
            mask msk = 0;
            for (int y : adj[x])
                msk |= child(y);

            // x is a leaf/branching node
            if (rb[x] != -1)
            {
                // record leaf/branching children mask
                ch[rb[x]] = msk;
                return in(rb[x]);
            }
            return msk;
        };
        child(r);

        // identify branching ancestors
        function<void(int, mask)> branch = [&] (int x, mask m)
        {
            // x is a branching or leaf node
            if (rb[x] != -1)
            {
                ba[rb[x]] = m;
                m |= in(rb[x]);
            }
            for (int y : adj[x])
                branch(y, m);
        };
        branch(r, 0);
    }
};

// TODO: This special case is UNUSED for now
int path_tree(array2d& dp, const tree& t1, const tree& t2, const int lx, int x, int y)
{
    // trivial case: empty tree
    if (x == lx)
        return 0;

    // memoization
    int& sol = dp(x, y);
    if (sol != -1)
        return sol;

    // TODO: select the correct child
    // x's only child
    int w = t1.adj[x].front();

    // delete path node
    sol = max(sol, path_tree(dp, t1, t2, lx, w, y));
    for (int z : t2.adj[y])
    {
        // match x and y
        sol = max(sol, 1 + path_tree(dp, t1, t2, lx, w, z));
        // delete tree node
        sol = max(sol, path_tree(dp, t1, t2, lx, x, z));
    }
    return sol;
}

// lx/ly are parents of roots of paths in t1 and t2 respectively
int path_path(array2d& dp, const tree& t1, const tree& t2, const int lx, const int ly, int x, int y)
{
    // trivial case: empty tree
    if (x == lx || y == ly)
        return 0;

    // memoization
    int& sol = dp(x, y);
    if (sol != -1)
        return sol;

    // delete t1 node
    sol = max(sol, path_path(dp, t1, t2, lx, ly, t1.par[x], y));
    // match x and y
    sol = max(sol, 1 + path_path(dp, t1, t2, lx, ly, t1.par[x], t2.par[y]));
    // delete t2 node
    sol = max(sol, path_path(dp, t1, t2, lx, ly, x, t2.par[y]));
    return sol;
}

// TODO: precompute these
int find_root(const tree& t, const mask m, int x)
{
    while (x != -1 && (t.rb[x] == -1 || (m & in(t.rb[x]))))
        x = t.par[x];
    return x;
}

// generic bipartite matching dp: D is weight matrix, k index on left side, m mask of right side
int bipartite(array2d& D, array2d& dp, int k, mask m)
{
    // nothing left to match
    if (k == D.n || m == 0)
        return 0;

    // memoization
    int& sol = dp(k, m);
    if (sol != -1)
        return sol;

    // try to match kth tree in t1 with any avaliable tree in t1
    for (mask mi = m; mi != 0; mi -= ls(mi))
        sol = max(sol, D(k, tz(mi)) + bipartite(D, dp, k + 1, m ^ ls(mi)));
    return sol;
};

// rx/ry are leaf/branching roots of forests
int forest_forest(const tree& t1, const tree& t2, mask rx, mask ry)
{
    // trivial case: empty forest
    if (rx == 0 || ry == 0)
        return 0;

    // memoization
    auto it = fdp[rx].find(ry);
    if (it != fdp[rx].end())
        return it->second;

    int& sol = fdp[rx][ry];

    // (branching) roots of trees in t1 and t2 respectively
    // for use in bipartite matching
    vector<int> va, vb;

    // remove a path in t1
    for (mask mp = rx; mp != 0; mp -= ls(mp))
    {
        // branching node index
        int u = tz(mp);
        // dont disconnect the current root-to-branching path
        if (rx & t1.ch[u])
            continue;

        // for bipartite matching consider only nodes which are leaf/branching in current subforest
        va.push_back(u);

        // don't remove leaves
        if (!t1.ch[u])
            continue;

        // exclude u, include its branching children, exclude its branching ancestors
        mask m = em(rx ^ ls(mp) ^ t1.ch[u], t1.ba[u]);
        sol = max(sol, forest_forest(t1, t2, m, ry));
    }

    // remove a path in t2
    for (mask mq = ry; mq != 0; mq -= ls(mq))
    {
        // branching node index
        int v = tz(mq);
        // dont disconnect the current root-to-branching path
        if (ry & t2.ch[v])
            continue;

        // for bipartite matching consider only nodes which are leaf/branching in current subforest
        vb.push_back(v);

        // don't remove leaves
        if (!t2.ch[v])
            continue;

        // exclude v, include its branching children, exclude its branching ancestors
        mask m = em(ry ^ ls(mq) ^ t2.ch[v], t2.ba[v]);
        sol = max(sol, forest_forest(t1, t2, rx, m));
    }

    // tree to tree distance table
    array2d D(va.size(), vb.size());

    // compute tree-to-tree distances
    for (int i = 0; i < va.size(); ++i)
    {
        for (int j = 0; j < vb.size(); ++j)
        {
            // branching node index
            int u = va[i], v = vb[j];
            // path  x..z in t1
            int z = t1.b[u];
            int x = find_root(t1, rx, z);
            // path y..w in t2
            int w = t2.b[v];
            int y = find_root(t2, ry, w);

            // corresponding distance table entry
            int& uvs = D(i, j);
            // select a branch at u
            for (mask s = t1.ch[u]; s != 0; s -= ls(s))
                uvs = max(uvs, forest_forest(t1, t2, in(u) | ls(s) | (rx & t1.ba[u]), in(v) | (ry & t2.ba[v])));
            // select a branch at v
            for (mask s = t2.ch[v]; s != 0; s -= ls(s))
                uvs = max(uvs, forest_forest(t1, t2, in(u) | (rx & t1.ba[u]), in(v) | ls(s) | (ry & t2.ba[v])));

            // initialize dp table for path to path matching
            array2d dp(t1.n, t2.n);
            // match the paths
            uvs = max(uvs, forest_forest(t1, t2, t1.ch[u], t2.ch[v]) + path_path(dp, t1, t2, x, y, z, w));
        }
    }

    // dp table for bipartite matching
    array2d bdp(va.size(), 1 << vb.size());

    // compute optimal forest-to-forest bipartite matching
    sol = max(sol, bipartite(D, bdp, 0, fm(vb.size())));
    return sol;
}

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        cout << "usage: " << argv[0] << " <tree1> <tree2>\n";
        return 0;
    }
    tree t1(argv[1]);
    tree t2(argv[2]);

    // first node in topological ordering is the root
    cout << forest_forest(t1, t2, 1, 1);
}
