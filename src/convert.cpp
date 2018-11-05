/*
    Copyright (C) 2018 Mislav Blažević

    This file is part of Trajan.

    Trajan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <string>
#include <sstream>
using namespace boost;
using namespace std;

struct vertex_label_t
{
    typedef vertex_property_tag kind;
};

struct edge_label_t
{
    typedef edge_property_tag kind;
};

typedef property <vertex_name_t, string, property <vertex_label_t, string> > vertex_p;
typedef property <edge_color_t, string, property <edge_label_t, string> > edge_p;
typedef property <graph_name_t, string> graph_p;
typedef adjacency_list <vecS, vecS, directedS, vertex_p, edge_p, graph_p> graph_t;

int main(int argc, char** argv)
{
    if (argc != 4)
    {
        cout << "usage: " << argv[0] << " <input.dot> <out-dag> <out-map>" << endl;
        return EXIT_FAILURE;
    }

    graph_t graph(0);
    dynamic_properties dp(ignore_other_properties);

    property_map<graph_t, vertex_name_t>::type name =
        get(vertex_name, graph);
    dp.property("node_id", name);

    property_map<graph_t, vertex_label_t>::type vertex_label =
        get(vertex_label_t(), graph);
    dp.property("label", vertex_label);

    property_map<graph_t, edge_color_t>::type color =
        get(edge_color, graph);
    dp.property("color", color);

    property_map<graph_t, edge_label_t>::type edge_label =
        get(edge_label_t(), graph);
    dp.property("label", edge_label);

    ref_property_map<graph_t*, string> gname(get_property(graph, graph_name));
    dp.property("name", gname);

    ifstream gvgraph(argv[1]);
    if (!read_graphviz(gvgraph, graph, dp, "node_id"))
    {
        cout << "error reading dot file" << endl;
        return EXIT_FAILURE;
    }

    ofstream fdag(argv[2]);
    graph_traits<graph_t>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei)
        fdag << name[target(*ei, graph)] << " " << name[source(*ei, graph)] << " default\n";

    ofstream fmap(argv[3]);
    graph_traits<graph_t>::vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(graph); vi != vi_end; ++vi)
    {
        string lab = vertex_label[*vi];
        size_t beg = lab.find('=') + 1;
        size_t end = lab.find('|');
        size_t len = end - beg;
        fmap << lab.substr(beg, len) << " " << name[*vi] << '\n';
    }

    return EXIT_SUCCESS;
}
