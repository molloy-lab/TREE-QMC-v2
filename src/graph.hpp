#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "utility.hpp"
#include "dict.hpp"
#include "tree.hpp"
#include "heuristics/maxcut/burer2002.h"
#include "problem/instance.h"
#include "problem/max_cut_instance.h"

class Graph {
    public:
        Graph(std::vector<Tree *> trees, Taxa &subset);
        ~Graph();
        std::string to_string();
        weight_t get_cut(std::vector<index_t> *A, std::vector<index_t> *B);
    private:
        index_t size;
        std::unordered_map<index_t, index_t> index2index;
        std::vector<index_t> indices;
        weight_t ***graph;
        weight_t sdp_cut(weight_t alpha, std::vector<index_t> *A, std::vector<index_t> *B);
};

#endif
