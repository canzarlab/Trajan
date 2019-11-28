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
#include <limits>
#include <fstream>
#include <cassert>
#include <cmath>
#include "read_csv.h"
#define die() assert(false)
#define dbg(x) cerr << #x << " = " << x << endl
using namespace std;

using uint = unsigned;
using mask = unsigned long long;

const int NINF = numeric_limits<int>::min();

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

    array2d(int n, int m, int val = NINF) : n(n), m(m)
    {
        v.resize(n * m, val);
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
using ff_state = tuple<int, int, int, vector<int>, vector<int>, vector<int>, vector<int>, vector<int>, array2d, array2d>;
using dp_table = unordered_map<mask, unordered_map<mask, ff_state>>;
using path_dp_table = map<tuple<int, int>, pair<array2d, array2d>>;
dp_table fdp; // forest-to-forest dp table
path_dp_table pdp; // path-to-path dp table

// dp transitions
enum : int
{
    DEL_1 = 0,
    DEL_2,
    MATCH
};

struct tree
{
    // mapping of node identifiers to nodes
    map<string, int> m;
    // mapping of nodes to node identifiers
    vector<string> rm;
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

        // create reverse mapping
        rm.resize(n);
        for (auto [name, idx] : m)
            rm[idx] = name;

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
    if (sol != NINF)
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
int path_path(array2d& dp, array2d& pdp, const array2d& matrix, const tree& t1, const tree& t2, const int lx, const int ly, int x, int y)
{
    // trivial case: empty path
    if (x == lx || y == ly)
        return 0;

    // memoization
    int& sol = dp(x, y);
    if (sol != NINF)
        return sol;

    // initialize score to negative infinity
    sol = NINF;

    // delete t1 node
    int nsol = path_path(dp, pdp, matrix, t1, t2, lx, ly, t1.par[x], y);
    if (nsol > sol)
        sol = nsol, pdp(x, y) = DEL_1;

    // match x and y
    nsol = matrix(x, y) + path_path(dp, pdp, matrix, t1, t2, lx, ly, t1.par[x], t2.par[y]);
    if (nsol > sol)
        sol = nsol, pdp(x, y) = MATCH;

    // delete t2 node
    nsol = path_path(dp, pdp, matrix, t1, t2, lx, ly, x, t2.par[y]);
    if (nsol > sol)
        sol = nsol, pdp(x, y) = DEL_2;

    return sol;
}

// find the root of the tree containing node x in the forest m
int find_root(const tree& t, const mask m, int x)
{
    // go up if x isn't a branching node or is a branching node included in the forest m
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
    if (sol != NINF)
        return sol;

    // don't match kth tree with anyone
    sol = bipartite(D, dp, pdp, k + 1, m);

    // try to match kth tree in t1 with any avaliable tree in t2
    for (mask mi = m; mi != 0; mi -= ls(mi))
    {
        int nsol = D(k, tz(mi)) + bipartite(D, dp, pdp, k + 1, m ^ ls(mi));
        if (nsol > sol)
            sol = nsol, pdp(k, m) = tz(mi);
    }
    return sol;
}

// forest-to-forest bipartite matching
int bipartite_aux(array2d& D, vector<int>& sol)
{
    // complexity of bipartite is O(n2^m) so make sure m < n
    bool transposed = false;
    if (D.n < D.m)
    {
        D = D.transpose();
        transposed = true;
    }

    // dp table for the matching
    array2d bdp(D.n, 1 << D.m), pdp(D.n, 1 << D.m, -1);

    // compute the matching
    int weight = bipartite(D, bdp, pdp, 0, fm(D.m));

    // store the matching
    // note: -1 is unmatched
    sol.resize(D.n, -1);
    mask m = fm(D.m);
    for (int k = 0; k < D.n && m > 0; ++k)
    {
        sol[k] = pdp(k, m);
        if (pdp(k, m) != -1)
            m ^= in(pdp(k, m));
    }

    // transpose the matching back
    if (transposed)
    {
        vector<int> zol(D.m, -1);
        for (int i = 0; i < D.n; ++i)
            if (sol[i] != -1)
                zol[sol[i]] = i;
        swap(sol, zol);
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
        return get<0>(it->second);

    // score, transition type, transition destination, bipartite matching
    // (branching) roots of trees in t1 and t2 respectively
    // roots of paths in va and vb respectively
    // bipartite matching transition table
    // branch choice table
    auto& [sol, trs, dst, bdp, va, vb, ra, rb, btr, bst] = fdp[rx][ry];

    // initialize score to negative infinity
    sol = NINF;

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
        int nsol = forest_forest(t1, t2, matrix, m, ry);
        if (nsol > sol)
            sol = nsol, trs = DEL_1, dst = u;
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
        int nsol = forest_forest(t1, t2, matrix, rx, m);
        if (nsol > sol)
            sol = nsol, trs = DEL_2, dst = v;
    }

    // tree to tree distance table
    array2d D(va.size(), vb.size());

    // initialize tree to tree transition tables
    btr = array2d(va.size(), vb.size());
    bst = array2d(va.size(), vb.size());

    // compute tree-to-tree weights for use in bipartite matching
    for (uint i = 0; i < va.size(); ++i)
    {
        for (uint j = 0; j < vb.size(); ++j)
        {
            // branching node index
            int u = va[i], v = vb[j];
            // path x..z in t1
            int z = t1.b[u];
            int x = ra[i];
            // path y..w in t2
            int w = t2.b[v];
            int y = rb[j];
            // corresponding distance table entry
            int& uvs = D(i, j);

            // select a branch at u
            for (mask s = t1.ch[u]; s != 0; s -= ls(s))
            {
                // on t1 side: include u, chosen branch s and keep branching ancestors of u
                // on t2 side: include v and keep branching ancestors of v
                int nuvs = forest_forest(t1, t2, matrix, in(u) | ls(s) | (rx & t1.ba[u]), in(v) | (ry & t2.ba[v]));
                if (nuvs > uvs)
                    uvs = nuvs, btr(i, j) = DEL_1, bst(i, j) = tz(s);
            }

            // select a branch at v
            for (mask s = t2.ch[v]; s != 0; s -= ls(s))
            {
                // on t1 side: include u and keep branching ancestors of u
                // on t2 side: include v, chosen branch s and keep branching ancestors of v
                int nuvs = forest_forest(t1, t2, matrix, in(u) | (rx & t1.ba[u]), in(v) | ls(s) | (ry & t2.ba[v]));
                if (nuvs > uvs)
                    uvs = nuvs, btr(i, j) = DEL_2, bst(i, j) = tz(s);
            }

            // initialize or fetch path-to-path matching dp table
            auto it = pdp.find({x, y});
            if (it == pdp.end())
                pdp[{x, y}] = {array2d(t1.n, t2.n), array2d(t1.n, t2.n)};
            auto& [dp, tdp] = pdp[{x, y}];

            // match the paths
            int nuvs = forest_forest(t1, t2, matrix, t1.ch[u], t2.ch[v]) + path_path(dp, tdp, matrix, t1, t2, x, y, z, w);
            if (nuvs > uvs)
                uvs = nuvs, btr(i, j) = MATCH, bst(i, j) = -1;
        }
    }
    // compute optimal forest-to-forest bipartite matching
    int nsol = bipartite_aux(D, bdp);
    if (nsol > sol)
        sol = nsol, trs = MATCH, dst = -1;

    return sol;
}

// recover path-to-path matching
void recover_matching_aux(const tree& t1, const tree& t2, vector<pair<int, int>>& matching, const array2d& tdp, const int lx, const int ly, int x, int y)
{
    if (x == lx || y == ly)
        return;

    switch (tdp(x, y))
    {
        case DEL_1:
            recover_matching_aux(t1, t2, matching, tdp, lx, ly, t1.par[x], y);
            break;
        case DEL_2:
            recover_matching_aux(t1, t2, matching, tdp, lx, ly, x, t2.par[y]);
            break;
        case MATCH:
            matching.emplace_back(x, y);
            recover_matching_aux(t1, t2, matching, tdp, lx, ly, t1.par[x], t2.par[y]);
            break;
    }
}

// rx/ry are leaf/branching roots of forests
void recover_matching(const tree& t1, const tree& t2, vector<pair<int, int>>& matching, mask rx, mask ry)
{
    // trivial case: empty forest
    if (rx == 0 || ry == 0)
        return;

    // score, transition type, transition destination, bipartite matching
    // (branching) roots of trees in t1 and t2 respectively
    // roots of paths in va and vb respectively
    // bipartite matching transition table
    // branch choice table
    const auto& [sol, trs, dst, bdp, va, vb, ra, rb, btr, bst] = fdp[rx][ry];
    switch (trs)
    {
        case DEL_1:
            recover_matching(t1, t2, matching, em(rx ^ in(dst) ^ t1.ch[dst], t1.ba[dst]), ry);
            break;
        case DEL_2:
            recover_matching(t1, t2, matching, rx, em(ry ^ in(dst) ^ t2.ch[dst], t2.ba[dst]));
            break;
        case MATCH:
            for (int i = 0; i < bdp.size(); ++i)
            {
                // ith tree isn't matched with anyone
                if (bdp[i] == -1)
                    continue;

                // ith tree is matched with jth tree
                int j = bdp[i];
                // branching node index
                int u = va[i], v = vb[j];
                // path  x..z in t1
                int z = t1.b[u];
                int x = ra[i];
                // path y..w in t2
                int w = t2.b[v];
                int y = rb[j];
                // branch choice
                int s = bst(i, j);

                switch (btr(i, j))
                {
                    case DEL_1:
                        recover_matching(t1, t2, matching, in(u) | in(s) | (rx & t1.ba[u]), in(v) | (ry & t2.ba[v]));
                        break;
                    case DEL_2:
                        recover_matching(t1, t2, matching, in(u) | (rx & t1.ba[u]), in(v) | in(s) | (ry & t2.ba[v]));
                        break;
                    case MATCH:
                        const auto& [dp, tdp] = pdp[{x, y}];
                        recover_matching_aux(t1, t2, matching, tdp, x, y, z, w);
                        recover_matching(t1, t2, matching, t1.ch[u], t2.ch[v]);
                        break;
                }
            }
            break;
    }
}

// transitive closure auxiliary
void trans_closure(const tree& t, int x, int r, array2d& D)
{
    D(r, x) = 1;
    for (int y : t.adj[x])
        trans_closure(t, y, r, D);
}

// fill transitive closure table
void make_trans(const tree& t, int x, array2d& D)
{
    trans_closure(t, x, x, D);
    for (int y : t.adj[x])
        make_trans(t, y, D);
}

// return true iff edges a and b can simultaneously be in an arboreal matching
bool valid(pair<int, int> a, pair<int, int> b, array2d& D1, array2d& D2)
{
    int i = get<0>(a), j = get<0>(b);
    int x = get<1>(a), y = get<1>(b);
    bool c2 = (D1(j, i) || D1(i, j)) == (D2(x, y) || D2(y, x));
    bool c1 = (D1(i, j) == D2(x, y)); // assuming c2 is satisfied
    bool c0 = (i != j && x != y);
    return c0 && c1 && c2;
}

// return true iff matching is valid for (t1, t2)
bool verify(const tree& t1, const tree& t2, vector<pair<int, int>>& matching)
{
    // ancestry relationship tables
    array2d D1(t1.n, t1.n, 0), D2(t2.n, t2.n, 0);
    make_trans(t1, t1.r, D1);
    make_trans(t2, t2.r, D2);

    // check for pairwise conflicts
    for (int i = 0; i < matching.size(); ++i)
        for (int j = i + 1; j < matching.size(); ++j)
            if (!valid(matching[i], matching[j], D1, D2))
            {
                auto [a, b] = matching[i];
                auto [x, y] = matching[j];
                dbg(t1.rm[a]);
                dbg(t2.rm[b]);
                dbg(t1.rm[x]);
                dbg(t2.rm[y]);
                return false;
            }
    return true;
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
    const double scale = 1000;
    for (int i = 0; i < t1.n + 1; ++i)
        for (int j = 0; j < t2.n + 1; ++j)
            minmatrix(i, j) = scale * cost_matrix[i][j];

    // convert into maximization problem
    array2d maxmatrix(t1.n, t2.n);
    for (int i = 0; i < t1.n; ++i)
        for (int j = 0; j < t2.n; ++j)
            maxmatrix(i, j) = minmatrix(i, t2.n) + minmatrix(t1.n, j) - minmatrix(i, j);

    // compute optimal tree to tree matching
    // note: first node in topological ordering is the root
    double weight = forest_forest(t1, t2, maxmatrix, 1, 1);

    // convert to minimization score
    weight = -weight;
    for (int i = 0; i < t1.n; ++i)
        weight += minmatrix(i, t2.n);
    for (int i = 0; i < t2.n; ++i)
        weight += minmatrix(t1.n, i);

    // scale back the score
    clog << weight / scale << endl;

    // recover the matching
    vector<pair<int, int>> matching;
    recover_matching(t1, t2, matching, 1, 1);

    // verify the matching
    assert(verify(t1, t2, matching));

    // print the matching
    for (auto [x, y] : matching)
        cout << t1.rm[x] << ' ' << t2.rm[y] << '\n';
}
