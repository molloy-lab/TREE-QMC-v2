#include "tree.hpp"

Node::Node(index_t index) {
    this->index = index;
    parent = NULL;
    size = -1;
}

Node::~Node() {
    for (Node *child : children) 
        delete child;
}

void Node::new_states(index_t size) {
    this->size = size;
    singlet = new weight_t[size + 1];
    for (index_t i = 0; i < size + 1; i ++) 
        singlet[i] = 0;
    doublet = new std::vector<std::pair<index_t, weight_t>>();
    /*
    doublet = new weight_t*[size + 1];
    for (index_t i = 0; i < size + 1; i ++) {
        doublet[i] = new weight_t[size + 1];
        for (index_t j = 0; j < size + 1; j ++) 
            doublet[i][j] = 0;
    }
    */
}

void Node::delete_states() {
    if (size >= 0) {
        delete [] singlet;
        delete doublet;
        /*
        for (index_t i = 0; i < size + 1; i ++) 
            delete [] doublet[i];
        delete [] doublet;
        */
        size = -1;
    }
}

weight_t Node::get_doublet(index_t a, index_t b) {
    index_t key = a * (size + 1) + b;
    index_t l = 0, r = doublet->size() - 1;
    while (l <= r) {
        index_t m = (l + r) >> 1;
        if (key == (*doublet)[m].first) return (*doublet)[m].second;
        if (key < (*doublet)[m].first) 
            r = m - 1;
        else 
            l = m + 1;
    }
    return 0;
}

void Node::add_doublet(index_t a, index_t b, weight_t c) {
    index_t key = a * (size + 1) + b;
    doublet->push_back(std::make_pair(key, c));
}

/*
weight_t Node::get_doublet(weight_t *singlet, weight_t s1, weight_t s2, index_t x, index_t y) {
    weight_t t0 = singlet[0], t1 = s1, t2 = s2;
    if (x != 0) {
        t1 -= singlet[x];
        t2 -= singlet[x] * singlet[x];
    }
    if (y != 0) {
        t1 -= singlet[y];
        t2 -= singlet[y] * singlet[y];
    }
    return (t1 * t1 - t2) / 2 + t1 * t0 + t0 * (t0 - 1) / 2;
}
*/

weight_t Node::get_doublet(weight_t *singlet, weight_t s1, weight_t s2, index_t x, index_t y) {
    if (x != 0) {
        weight_t sx = singlet[x];
        s1 -= sx;
        s2 -= sx * sx;
    }
    if (y != 0) {
        weight_t sy = singlet[y];
        s1 -= sy;
        s2 -= sy * sy;
    }
    weight_t s0 = singlet[0], s = s1 + s0;
    return s * s - s2 - s0;
}
