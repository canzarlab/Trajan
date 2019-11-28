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
#include <cmath>
#include "read_csv.h"
#define die() assert(false)
#define dbg(x) cerr << #x << " = " << x << endl
using namespace std;

using uint = unsigned;
using mask = unsigned long long;

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

// should be faster than vector<vector<int>>
struct array2d
{
    int n, m;
    vector<int> v;

    array2d() : n(0), m(0)
    {
    }

    array2d(int n, int m) : n(n), m(m)
    {
        v.resize(n * m, -1e9);
    }

    int& operator()(int x, int y)
    {
        return v[x * m + y];
    }

    const int& operator()(int x, int y) const
    {
        return v[x * m + y];
    }

    array2d transpose() const
    {
        array2d D(m, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                D(j, i) = (*this)(i, j);
        return D;
    }
};

// TODO: replace with gp_hash_table from __gnu_pbds or some other faster hash table
using dp_table = unordered_map<mask, unordered_map<mask, int>>;
using path_dp_table = map<tuple<int, int, int, int>, array2d>;
using bi_dp_table = unordered_map<mask, unordered_map<mask, vector<int>>>;
dp_table fdp; // forest-to-forest dp table
path_dp_table pdp; // path-to-path dp table
bi_dp_table bdp; // bipartite matching dp table

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

    tree(const char* filename, const char* mapname) : n(0), nb(0), r(-1)
    {
        ifstream file(filename), mapfile(mapname);
        if (!file || !mapfile)
            die();

        // map string indices to integers
        int node_idx; string node_name;
        while (mapfile >> node_idx >> node_name)
            m[node_name] = node_idx;
        n = m.size();

        // nroot[x] == 1 means x is not the root
        vector<int> nroot(MAXN);

        // read from file
        string sa, sb, type;
        while (file >> sa >> sb >> type)
        {
            if (type != "default")
                die();

            // directed edge sb->sa
            int x = m[sa], y = m[sb];
            adj[y].push_back(x);
            nroot[x] = 1;
        }

        // tree too big
        if (m.size() > MAXN)
            die();

        // find root
        for (uint i = 0; i < m.size(); ++i)
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
    // trivial case: empty path
    if (x == lx)
        return 0;

    // memoization
    int& sol = dp(x, y);
    if (sol != -1e9)
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
int path_path(array2d& dp, const array2d& matrix, const tree& t1, const tree& t2, const int lx, const int ly, int x, int y)
{
    // trivial case: empty path
    if (x == lx || y == ly)
        return 0;

    // memoization
    int& sol = dp(x, y);
    if (sol != -1e9)
        return sol;

    // delete t1 node
    sol = max(sol, path_path(dp, matrix, t1, t2, lx, ly, t1.par[x], y));
    // match x and y
    sol = max(sol, matrix(x, y) + path_path(dp, matrix, t1, t2, lx, ly, t1.par[x], t2.par[y]));
    // delete t2 node
    sol = max(sol, path_path(dp, matrix, t1, t2, lx, ly, x, t2.par[y]));
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
int bipartite(const array2d& D, array2d& dp, array2d& pdp, int k, mask m)
{
    // nothing left to match
    if (k == D.n || m == 0)
        return 0;

    // memoization
    int& sol = dp(k, m);
    if (sol != -1e9)
        return sol;

    // don't match kth tree with anyone
    sol = bipartite(D, dp, pdp, k + 1, m), pdp(k, m) = 0;

    // try to match kth tree in t1 with any avaliable tree in t2
    for (mask mi = m; mi != 0; mi -= ls(mi))
    {
        int nsol = D(k, tz(mi)) + bipartite(D, dp, pdp, k + 1, m ^ ls(mi));
        if (nsol > sol)
            sol = nsol, pdp(k, m) = tz(mi);
    }
    return sol;
}

int bipartite_aux(array2d& D, vector<int>& sol)
{
    // complexity of bipartite is O(n2^m) so make sure m < n
    if (D.n < D.m)
        D = D.transpose();

    // dp table for the matching
    array2d bdp(D.n, 1 << D.m), pdp(D.n, 1 << D.m);

    // compute the matching
    int weight = bipartite(D, bdp, pdp, 0, fm(D.m));

    // store the matching
    sol.resize(D.n);
    int m = fm(D.m);
    for (int k = 0; k < D.n; ++k)
    {
        sol[k] = pdp(k, m) - 1;
        m ^= in(pdp(k, m));
    }
    return weight;
}

// rx/ry are leaf/branching roots of forests
int forest_forest(const tree& t1, const tree& t2, const array2d& matrix, mask rx, mask ry)
{
    // trivial case: empty forest
    if (rx == 0 || ry == 0)
        return 0;

    // memoization
    auto it = fdp[rx].find(ry);
    if (it != fdp[rx].end())
        return it->second;

    int& sol = fdp[rx][ry] = -1e9;

    // (branching) roots of trees in t1 and t2 respectively
    vector<int> va, vb;
    // roots of paths in va and vb respectively
    vector<int> ra, rb;

    // don't waste time on reallocation
    va.reserve(MAXLB);
    vb.reserve(MAXLB);
    ra.reserve(MAXLB);
    rb.reserve(MAXLB);

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

        // path  x..z in t1
        int z = t1.b[u];
        int x = find_root(t1, rx, z);
        ra.push_back(x);

        // don't remove leaves
        if (!t1.ch[u])
            continue;

        // exclude u, include its branching children, exclude its branching ancestors
        mask m = em(rx ^ ls(mp) ^ t1.ch[u], t1.ba[u]);
        sol = max(sol, forest_forest(t1, t2, matrix, m, ry));
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

        // path y..w in t2
        int w = t2.b[v];
        int y = find_root(t2, ry, w);
        rb.push_back(y);

        // don't remove leaves
        if (!t2.ch[v])
            continue;

        // exclude v, include its branching children, exclude its branching ancestors
        mask m = em(ry ^ ls(mq) ^ t2.ch[v], t2.ba[v]);
        sol = max(sol, forest_forest(t1, t2, matrix, rx, m));
    }

    // tree to tree distance table
    array2d D(va.size(), vb.size());

    // compute tree-to-tree distances
    for (uint i = 0; i < va.size(); ++i)
    {
        for (uint j = 0; j < vb.size(); ++j)
        {
            // branching node index
            int u = va[i], v = vb[j];
            // path  x..z in t1
            int z = t1.b[u];
            int x = ra[i];
            // path y..w in t2
            int w = t2.b[v];
            int y = rb[j];

            // corresponding distance table entry
            int& uvs = D(i, j);
            // select a branch at u
            for (mask s = t1.ch[u]; s != 0; s -= ls(s))
                uvs = max(uvs, forest_forest(t1, t2, matrix, in(u) | ls(s) | (rx & t1.ba[u]), in(v) | (ry & t2.ba[v])));
            // select a branch at v
            for (mask s = t2.ch[v]; s != 0; s -= ls(s))
                uvs = max(uvs, forest_forest(t1, t2, matrix, in(u) | (rx & t1.ba[u]), in(v) | ls(s) | (ry & t2.ba[v])));

            // initialize or fetch path-to-path matching dp table
            auto it = pdp.find({x, y, z, w});
            if (it == pdp.end())
                pdp[{x, y, z, w}] = array2d(t1.n, t2.n);

            array2d& dp = pdp[{x, y, z, w}];
            // match the paths
            uvs = max(uvs, forest_forest(t1, t2, matrix, t1.ch[u], t2.ch[v]) + path_path(dp, matrix, t1, t2, x, y, z, w));
        }
    }
    // compute optimal forest-to-forest bipartite matching
    return sol = max(sol, bipartite_aux(D, bdp[rx][ry]));
}

int main(int argc, char** argv)
{
    if (argc != 6)
    {
        cout << "usage: " << argv[0] << " <tree1> <map1> <tree2> <map2> <matrix>\n";
        return 0;
    }
    tree t1(argv[1], argv[2]);
    tree t2(argv[3], argv[4]);
    vector<vector<double>> cost_matrix = CSVReader(argv[5]).getDoubleData();

    // last row/column is deletion cost
    array2d minmatrix(t1.n + 1, t2.n + 1);

    // convert scores to integers
    const double scale = 1000.0;
    for (int i = 0; i < t1.n + 1; ++i)
        for (int j = 0; j < t2.n + 1; ++j)
            minmatrix(i, j) = scale * cost_matrix[i][j];

    // convert into maximization problem
    array2d maxmatrix(t1.n, t2.n);
    for (int i = 0; i < t1.n; ++i)
        for (int j = 0; j < t2.n; ++j)
            maxmatrix(i, j) = minmatrix(i, t2.n) + minmatrix(t1.n, j) - minmatrix(i, j);

    // first node in topological ordering is the root
    double weight = forest_forest(t1, t2, maxmatrix, 1, 1);

    // convert to minimization score
    weight = -weight;
    for (int i = 0; i < t1.n; ++i)
        weight += minmatrix(i, t2.n);
    for (int i = 0; i < t2.n; ++i)
        weight += minmatrix(t1.n, i);
    // scale back the score
    cout << weight / scale << endl;
}
